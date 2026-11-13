#!/usr/bin/env python3
"""Report LINgroup clade counts, direct assignments, and assigned sequence bases."""

import argparse
import ast
import csv
from collections import defaultdict

import pandas as pd


def lin_tokens(value):
    return tuple(token.strip().strip("'\"") for token in str(value).strip().strip('[]').split(','))


def taxid_tokens(value):
    text = str(value).strip()
    if text.startswith('['):
        return tuple(int(x) for x in ast.literal_eval(text))
    return tuple(int(x.strip()) for x in text.split(','))


def build_taxonomy(data):
    """Index complete LIN components and each custom taxid's ancestor path."""
    prefixes = {}
    paths = {}
    for lin, ids in zip(data['LIN'], data['taxid_LIN']):
        tokens = lin_tokens(lin)
        taxids = taxid_tokens(ids)
        if len(tokens) != len(taxids):
            raise ValueError('LIN and taxid_LIN lengths differ in taxonomy/data.txt')
        for depth, taxid in enumerate(taxids, 1):
            prefix = tokens[:depth]
            path = taxids[:depth]
            if prefix in prefixes and prefixes[prefix] != path:
                raise ValueError('Inconsistent taxonomy for LIN prefix ' + ','.join(prefix))
            if taxid in paths and paths[taxid] != path:
                raise ValueError('Inconsistent ancestor paths for taxid ' + str(taxid))
            prefixes[prefix] = path
            paths[taxid] = path
    return prefixes, paths


def read_counts(report_file):
    clade, direct = {}, {}
    with open(report_file) as handle:
        for row in csv.reader(handle, delimiter='\t'):
            if not row:
                continue
            if len(row) < 6:
                raise ValueError('Malformed Kraken report: expected six columns')
            taxid = int(row[4])
            if taxid in clade:
                raise ValueError('Duplicate taxid in Kraken report: ' + str(taxid))
            clade[taxid], direct[taxid] = int(row[1]), int(row[2])
    return clade, direct


def assigned_lengths(output_file, paths, group_taxids):
    """Stream Kraken output; include descendants and both mates' base counts."""
    lengths = defaultdict(int)
    total_length = 0
    with open(output_file) as handle:
        for row in csv.reader(handle, delimiter='\t'):
            if not row or row[0] != 'C':
                continue
            if len(row) < 5:
                raise ValueError('Malformed classified row in Kraken output')
            taxid = int(row[2])
            bases = sum(int(x) for x in row[3].split('|'))
            total_length += bases
            for ancestor in set(paths.get(taxid, (taxid,))):
                if ancestor in group_taxids:
                    lengths[ancestor] += bases
    return lengths, total_length


def make_report(lin_file, data_file, report_file, output_file):
    groups = pd.read_csv(lin_file, sep='\t', dtype=str, keep_default_na=False)
    data = pd.read_csv(data_file, sep='\t', dtype=str, keep_default_na=False)
    prefixes, paths = build_taxonomy(data)
    clade, direct = read_counts(report_file)
    # Sum direct counts once per taxon. Custom LINtax roots need not be NCBI taxa.
    total_reads = sum(count for taxid, count in direct.items() if taxid != 0)
    group_paths = [prefixes.get(lin_tokens(prefix), ()) for prefix in groups['LINgroup_prefix']]
    group_taxids = {path[-1] for path in group_paths if path}
    lengths, total_length = assigned_lengths(output_file, paths, group_taxids)

    rows = []
    for (_, group), path in zip(groups.iterrows(), group_paths):
        taxid = path[-1] if path else None
        reads = clade.get(taxid, 0)
        unique = direct.get(taxid, 0)
        row = group.to_dict()
        row.update({
            'Assigned_reads': reads,
            'Percentage_assigned_reads': 100.0 * reads / total_reads if total_reads else 0.0,
            # Legacy column name: this is a direct group-node count, not a
            # claim that these reads uniquely identify a strain.
            'Unique_Assigned_reads': unique,
            'Direct_Assigned_reads': unique,
            'Percentage_unique_assigned_reads': 100.0 * unique / total_reads if total_reads else 0.0,
            'Total_reads_length': lengths[taxid] if taxid is not None else 0,
            'LINtax_note': '' if path else 'group_prefix_missing_from_taxonomy',
        })
        rows.append(row)
    total_unique = sum(direct.get(taxid, 0) for taxid in group_taxids)
    total = {column: '' for column in groups.columns}
    total.update({
        'LINgroup_Name': 'Total_reads',
        'Assigned_reads': total_reads,
        'Percentage_assigned_reads': 100.0 if total_reads else 0.0,
        'Unique_Assigned_reads': total_unique,
        'Direct_Assigned_reads': total_unique,
        'Percentage_unique_assigned_reads': 100.0 * total_unique / total_reads if total_reads else 0.0,
        'Total_reads_length': total_length,
        'LINtax_note': 'classified_fragments; unique_count_is_direct_group_node_assignments',
    })
    rows.append(total)
    return pd.DataFrame(rows)


def filter_report(frame, min_assigned_reads=2):
    if min_assigned_reads < 1:
        raise ValueError('Minimum assigned reads must be at least 1')
    return frame[(frame['LINgroup_Name'] == 'Total_reads') |
                 (frame['Assigned_reads'] >= min_assigned_reads)].copy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lin_file', required=True)
    parser.add_argument('--data_file', required=True)
    parser.add_argument('--in_file_report', required=True)
    parser.add_argument('--in_file_output', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--raw-output', help='Optional unfiltered group-count table')
    parser.add_argument('--min-assigned-reads', type=int, default=2,
                        help='Minimum clade fragments for a reported group (default: 2)')
    args = parser.parse_args()
    frame = make_report(args.lin_file, args.data_file, args.in_file_report, args.in_file_output)
    filtered = filter_report(frame, args.min_assigned_reads)
    if args.raw_output:
        frame.to_csv(args.raw_output, index=False)
    filtered.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
