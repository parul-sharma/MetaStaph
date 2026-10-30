#!/usr/bin/env python3
"""Summarize Kraken2 Staphylococcus aureus reads and preserve tool evidence.

Classified reads are read from the Kraken2 report's clade count for the target
taxid. Tool hits remain uncollapsed in the evidence TSV; no species-specific
strain consensus is inferred here.
"""

import argparse
import re
from pathlib import Path
from typing import Dict

import pandas as pd

TOOLS = ["StrainGE", "StrainScan", "Sourmash", "LINtax"]


def normalize_hit(s: str) -> str:
    if s is None:
        return ""
    x = str(s).strip().lower()
    x = x.replace("-", "_")
    x = x.replace("/", "_")
    x = re.sub(r"\s+", " ", x).strip()
    x = x.replace("__clade", "_clade")
    x = x.replace("-clade", "_clade")
    x = re.sub(r"_+", "_", x)
    return x


def parse_kraken_report(report_path: str, target_taxid: int) -> Dict[str, object]:
    out = {"Target_reads": 0, "Kraken_note": ""}
    if not report_path:
        out["Kraken_note"] = "kraken_report_missing"
        return out

    path = Path(report_path)
    if not path.exists() or path.stat().st_size == 0:
        out["Kraken_note"] = "kraken_report_empty"
        return out

    matches = []
    for line in path.read_text(errors="replace").splitlines():
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        try:
            taxid = int(fields[4])
            clade_reads = int(fields[1])
        except ValueError:
            continue
        if taxid == target_taxid:
            matches.append(clade_reads)

    if len(matches) == 1:
        out["Target_reads"] = matches[0]
    elif len(matches) > 1:
        out["Kraken_note"] = "target_taxid_duplicated_in_report"
    else:
        out["Kraken_note"] = "target_taxid_not_found"
    return out

def _row(tool: str, raw: str, hit_rank: int, evidence_value, evidence_metric: str) -> Dict[str, object]:
    norm = normalize_hit(raw)
    return {
        "Tool": tool,
        "Raw_hit": raw,
        "Normalized_hit": norm,
        "Hit_rank": hit_rank,
        "Evidence_value": evidence_value,
        "Evidence_metric": evidence_metric,
    }


def read_strainge(strainge_file: str) -> pd.DataFrame:
    rows = []
    tool = "StrainGE"
    if not strainge_file:
        return pd.DataFrame(rows)
    p = Path(strainge_file)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame(rows)
    try:
        df = pd.read_csv(p, sep="\t")
    except Exception:
        return pd.DataFrame(rows)
    name_col = next((c for c in ["strain", "reference", "ref", "name", "strain_name"] if c in df.columns), None)
    score_col = next((c for c in ["gcov", "cov", "coverage"] if c in df.columns), None)
    if name_col is None:
        return pd.DataFrame(rows)
    if score_col:
        df[score_col] = pd.to_numeric(df[score_col], errors="coerce").fillna(0.0)
        df = df.sort_values(score_col, ascending=False)
    for i, (_, rec) in enumerate(df.iterrows(), start=1):
        raw = str(rec[name_col])
        ev_val = float(rec[score_col]) if score_col else ""
        rows.append(_row(tool, raw, i, ev_val, score_col if score_col else ""))
    return pd.DataFrame(rows)


def read_strainscan(strainscan_dir: str) -> pd.DataFrame:
    rows = []
    tool = "StrainScan"
    if not strainscan_dir:
        return pd.DataFrame(rows)
    p = Path(strainscan_dir) / "final_report.txt"
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame(rows)
    try:
        df = pd.read_csv(p, sep="\t")
    except Exception:
        return pd.DataFrame(rows)
    if "Strain_Name" not in df.columns:
        return pd.DataFrame(rows)
    score_col = "Predicted_Depth (Enet)" if "Predicted_Depth (Enet)" in df.columns else None
    if score_col:
        df[score_col] = pd.to_numeric(df[score_col], errors="coerce").fillna(0.0)
        df = df.sort_values(score_col, ascending=False)
    for i, (_, rec) in enumerate(df.iterrows(), start=1):
        raw = str(rec["Strain_Name"])
        ev_val = float(rec[score_col]) if score_col else ""
        rows.append(_row(tool, raw, i, ev_val, score_col if score_col else ""))
    return pd.DataFrame(rows)


def read_sourmash(sourmash_file: str) -> pd.DataFrame:
    rows = []
    tool = "Sourmash"
    if not sourmash_file:
        return pd.DataFrame(rows)
    p = Path(sourmash_file)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame(rows)
    try:
        df = pd.read_csv(p, sep=",")
    except Exception:
        return pd.DataFrame(rows)
    if "name" not in df.columns:
        return pd.DataFrame(rows)
    score_col = "f_orig_query" if "f_orig_query" in df.columns else None
    if score_col:
        df[score_col] = pd.to_numeric(df[score_col], errors="coerce").fillna(0.0)
        df = df.sort_values(score_col, ascending=False)
    for i, (_, rec) in enumerate(df.iterrows(), start=1):
        raw = str(rec["name"])
        ev_val = float(rec[score_col]) if score_col else ""
        rows.append(_row(tool, raw, i, ev_val, score_col if score_col else ""))
    return pd.DataFrame(rows)


def read_lintax(linreport_file: str) -> pd.DataFrame:
    rows = []
    tool = "LINtax"
    if not linreport_file:
        return pd.DataFrame(rows)
    p = Path(linreport_file)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame(rows)
    try:
        df = pd.read_csv(p, sep=",")
    except Exception:
        return pd.DataFrame(rows)
    required = {"LINgroup_Name", "Assigned_reads"}
    if not required.issubset(set(df.columns)):
        return pd.DataFrame(rows)
    if "Percentage_assigned_reads" in df.columns:
        df["Percentage_assigned_reads"] = pd.to_numeric(df["Percentage_assigned_reads"], errors="coerce").fillna(0.0)
        df = df[(df["Percentage_assigned_reads"] > 0)]
    df["Assigned_reads"] = pd.to_numeric(df["Assigned_reads"], errors="coerce").fillna(0).astype(int)
    df = df.sort_values("Assigned_reads", ascending=False)
    for i, (_, rec) in enumerate(df.iterrows(), start=1):
        raw = str(rec["LINgroup_Name"])
        rows.append(_row(tool, raw, i, int(rec["Assigned_reads"]), "Assigned_reads"))
    return pd.DataFrame(rows)


def main(kraken_report: str, target_taxid: int, strainge_file: str, linreport_file: str,
         sourmash_file: str, strainscan_dir: str, sample_id: str,
         output_summary: str, output_evidence: str):
    kraken = parse_kraken_report(kraken_report, target_taxid)
    target_reads = int(kraken["Target_reads"])

    evidence_frames = [
        read_strainge(strainge_file),
        read_strainscan(strainscan_dir),
        read_sourmash(sourmash_file),
        read_lintax(linreport_file),
    ]
    frames = [frame for frame in evidence_frames if frame is not None and not frame.empty]
    if frames:
        evidence = pd.concat(frames, ignore_index=True)
    else:
        evidence = pd.DataFrame({
            "Tool": [], "Raw_hit": [], "Normalized_hit": [], "Hit_rank": [],
            "Evidence_value": [], "Evidence_metric": [],
        })
    evidence.insert(0, "Sample", sample_id)

    hit_counts = {tool: int((evidence["Tool"] == tool).sum()) for tool in TOOLS}
    tools_reporting = [tool for tool in TOOLS if hit_counts[tool] > 0]
    notes = []
    if kraken["Kraken_note"]:
        notes.append(kraken["Kraken_note"])
    if target_reads > 0 and not tools_reporting:
        notes.append("target_reads_present_but_no_strain_tool_hits")
    if target_reads == 0 and tools_reporting:
        notes.append("strain_tool_hits_present_but_no_target_reads")

    summary = pd.DataFrame([{
        "Sample": sample_id,
        "Target_species": "Staphylococcus aureus",
        "Target_taxid": target_taxid,
        "Target_species_reads": target_reads,
        "Target_species_detected": int(target_reads > 0),
        "Kraken_note": kraken["Kraken_note"] or "-",
        "Tools_reporting_n": f"{len(tools_reporting)}/{len(TOOLS)}",
        "Tools_reporting": ";".join(tools_reporting) if tools_reporting else "-",
        "Tool_hit_counts": ";".join(f"{tool}:{hit_counts[tool]}" for tool in TOOLS),
        "Total_tool_hits": int(sum(hit_counts.values())),
        "Notes": ";".join(notes) if notes else "-",
    }])

    Path(output_summary).parent.mkdir(parents=True, exist_ok=True)
    Path(output_evidence).parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_summary, index=False)
    evidence.to_csv(output_evidence, sep="\t", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="MetaStaph: summarize target reads from Kraken2 report and preserve tool evidence.")
    ap.add_argument("-kr", "--kraken_report", required=True, help="Kraken2 report generated with the LINtax database.")
    ap.add_argument("--target_taxid", type=int, default=1280, help="NCBI taxid for the target species (default: 1280).")
    ap.add_argument("-st", "--strainge", required=False, help="Path to StrainGE strains TSV.")
    ap.add_argument("-lr", "--linreport", required=False, help="Path to LINtax report CSV.")
    ap.add_argument("-sm", "--sourmash", required=False, help="Path to sourmash gather CSV.")
    ap.add_argument("-sc", "--strainscan", required=True, help="Path to StrainScan output directory (contains final_report.txt).")
    ap.add_argument("-s", "--sample", required=True, help="Sample ID.")
    ap.add_argument("-o", "--output", required=True, help="Summary CSV output path.")
    ap.add_argument("-e", "--evidence", required=True, help="Evidence TSV output path.")
    args = ap.parse_args()
    main(args.kraken_report, args.target_taxid, args.strainge, args.linreport,
         args.sourmash, args.strainscan, args.sample, args.output, args.evidence)
