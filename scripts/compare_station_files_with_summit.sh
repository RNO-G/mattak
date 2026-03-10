#!/bin/bash
# compare_runs.sh
# Usage: ./compare_runs.sh <input_file> <dir_compare> [-q]
#
# Compares non-empty runs from input_file against dir_compare.
# A run is considered non-empty if its size exceeds EMPTY_THRESHOLD_KB.
# get  input_file by ssh into station then du -sh /data/logs/* > input_file
# Arguments:
#   input_file   - tab-separated file with columns: <size> <run_name>
#   dir_compare  - directory to check runs against
#   -q           - quiet mode, only print MISSING and EMPTY runs

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 <input_file> <dir_compare> [-q]"
    exit 1
fi

INPUT_FILE="$1"
DIR_COMPARE="$2"
QUIET=false
[[ "$3" == "-q" ]] && QUIET=true

EMPTY_THRESHOLD_KB=24
ok=0; missing=0; empty=0

mapfile -t lines < "$INPUT_FILE"

for line in "${lines[@]}"; do
    read -r size run_name <<< "$line"
    [[ "$size" == "4.0K" || "$size" == "4K" || "$size" == "24K" ]] && continue

    run_path="$DIR_COMPARE/$run_name"
    if [[ ! -d "$run_path" ]]; then
        echo "MISSING: $run_name"; ((missing++)); continue
    fi

    run_kb=$(du -sk "$run_path" 2>/dev/null | cut -f1)
    run_mb=$(awk "BEGIN {printf \"%.1f\", $run_kb/1024}")

    if (( run_kb <= EMPTY_THRESHOLD_KB )); then
        echo "EMPTY:   $run_name (${run_mb}M)"; ((empty++))
    else
        [[ "$QUIET" == false ]] && echo "OK:      $run_name (${run_mb}M)"; ((ok++))
    fi
done

echo -e "\nOK: $ok | Missing: $missing | Empty: $empty"
