#!/usr/bin/env bash
# batch_profile_plots.sh
# Recursively finds every sample_diameter_profile.json under a results directory
# and generates a sample_profile_I_Apm.html in the same folder.
#
# Usage:
#   ./batch_profile_plots.sh                          # uses ./results/
#   ./batch_profile_plots.sh /path/to/results         # explicit path
#   ./batch_profile_plots.sh --force /path/to/results # overwrite existing plots

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLOTTER="$SCRIPT_DIR/plot_runlog_compact.py"

FORCE=0
RESULTS_DIR=""

for arg in "$@"; do
    if [[ "$arg" == "--force" ]]; then
        FORCE=1
    else
        RESULTS_DIR="$arg"
    fi
done

if [[ -z "$RESULTS_DIR" ]]; then
    RESULTS_DIR="$SCRIPT_DIR/results"
fi

if [[ ! -d "$RESULTS_DIR" ]]; then
    echo "ERROR: results directory not found: $RESULTS_DIR"
    exit 1
fi

echo "Scanning: $RESULTS_DIR"
echo ""

FOUND=0
SKIPPED=0
DONE=0
FAILED=0

while IFS= read -r -d '' json_file; do
    dir="$(dirname "$json_file")"
    out="$dir/sample_profile_I_Apm.html"

    FOUND=$((FOUND + 1))

    if [[ -f "$out" && "$FORCE" -eq 0 ]]; then
        echo "[skip]  $(basename "$dir")"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    echo "[plot]  $(basename "$dir")"
    err=$(python "$PLOTTER" sample-profile \
        --json "$json_file" \
        --x-units mm \
        --y I_Apm \
        --out "$out" 2>&1)
    if [[ $? -eq 0 ]]; then
        DONE=$((DONE + 1))
    else
        echo "        ERROR: $err"
        FAILED=$((FAILED + 1))
    fi

done < <(find "$RESULTS_DIR" -name "sample_diameter_profile.json" -print0 | sort -z)

echo ""
echo "Done. found=$FOUND plotted=$DONE skipped=$SKIPPED failed=$FAILED"
echo "(use --force to overwrite existing plots)"
