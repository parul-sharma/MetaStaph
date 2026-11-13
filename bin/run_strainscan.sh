#!/usr/bin/env bash
# Run from a Nextflow task directory, keeping diagnostics with the results.
set -euo pipefail
sample_id=$1
read1=$2
read2=$3
database=$4
output="$PWD/strainscan_output"
mkdir -p "$output"
# Some installed launchers change directory. Keep all tool paths absolute.
[[ "$read1" = /* ]] || read1="$PWD/$read1"
[[ "$read2" = /* ]] || read2="$PWD/$read2"
status=0
strainscan -i "$read1" -j "$read2" -d "$database" -o "$output" \
    > "$output/strainscan.stdout.log" 2> "$output/strainscan.stderr.log" || status=$?

result=failed
note="exit_code=$status; see strainscan.stdout.log and strainscan.stderr.log"
report="$output/final_report.txt"
if [[ -s "$report" ]] && awk -F '\t' '
    NR == 1 { for (i=1; i<=NF; i++) if ($i == "Strain_Name") found=1 }
    END { exit !found }
' "$report"; then
    if [[ "$status" -ne 0 ]]; then
        result=report_with_error
        note="final_report.txt present but exit_code=$status; inspect logs for incomplete execution"
    elif [[ $(wc -l < "$report") -gt 1 ]]; then
        result=success
        note=completed
    else
        result=no_hits
        note="final_report.txt contains a header only"
    fi
elif [[ "$status" -eq 0 ]]; then
    result=no_report
    note="exit_code=0 but no valid final_report.txt; inspect logs for no-cluster messages"
fi
printf 'sample_id\ttool\tstatus\tnote\n%s\tStrainScan\t%s\t%s\n' \
    "$sample_id" "$result" "$note" > "$output/strainscan_status.tsv"
if [[ "$result" == failed || "$result" == report_with_error ]]; then
    printf 'StrainScan %s: %s\n' "$sample_id" "$note" >&2
    tail -n 30 "$output/strainscan.stderr.log" >&2
fi
