#!/usr/bin/env bash
set -euo pipefail
sample_id=$1
sample_hdf5=$2
database=$3

make_placeholders() {
    # Header-only output carries zero hits; never invent an NA strain hit.
    printf 'strain\tgcov\n' > "$sample_id.strains.tsv"
    printf 'metric\tvalue\nstatus\t%s\n' "$1" > "$sample_id.stats.tsv"
}

if [[ ! -s "$sample_hdf5" ]]; then
    make_placeholders missing_or_empty_sample_hdf5
    exit 0
fi

# StrainGST uses Path.with_suffix(). A fixed prefix avoids losing dotted IDs.
status=0
straingst run -O -o strainge_result "$database" "$sample_hdf5" \
    > "$sample_id.strainge.stdout.log" 2> "$sample_id.strainge.stderr.log" || status=$?
if [[ "$status" -ne 0 ]]; then
    cat "$sample_id.strainge.stderr.log" >&2
    make_placeholders "straingst_failed_exit_$status"
elif [[ -s strainge_result.strains.tsv && -s strainge_result.stats.tsv ]]; then
    mv strainge_result.strains.tsv "$sample_id.strains.tsv"
    mv strainge_result.stats.tsv "$sample_id.stats.tsv"
else
    make_placeholders missing_or_empty_straingst_outputs
fi
