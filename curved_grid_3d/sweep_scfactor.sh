#!/usr/bin/env bash
# Sweep drift space-charge compensation. Submits one SLURM job per SC_FACTOR
# value, each writing to its own results dir. The case file reads SC_FACTOR and
# RESULTS_DIR from the environment (${VAR:-default}), which SLURM exports.
#
#   bash sweep_scfactor.sh [case_file.sh]
set -euo pipefail
cd "$(dirname "$0")"
CASE="${1:-case_baseline.sh}"
for sc in 1.0 0.5 0.0005; do
    tag=$(echo "$sc" | sed 's/\./p/')
    echo "submit SC_FACTOR=$sc -> results_3d_sc$tag"
    SC_FACTOR="$sc" RESULTS_DIR="results_3d_sc$tag" sbatch run_grid_3d.slurm "$CASE"
done
echo "compare afterward: for d in results_3d_sc*; do echo \$d; grep -o '\"merged_half_angle_deg\": [0-9.]*' \$d/diagnostics.json; done"
