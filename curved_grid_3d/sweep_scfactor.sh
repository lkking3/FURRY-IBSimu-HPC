#!/usr/bin/env bash
# Sweep drift space-charge compensation SC_FACTOR.
# The FIRST job builds the shared geometry cache (geom_cache.dat); the rest run
# with --dependency=afterok so they reuse it instead of rebuilding the mesh.
#   bash sweep_scfactor.sh [case_file.sh]
set -euo pipefail
cd "$(dirname "$0")"
CASE="${1:-case_baseline.sh}"
VALS=(0.0005 1.0 0.5)        # first value builds the cache
first=""
for sc in "${VALS[@]}"; do
    tag=$(echo "$sc" | sed 's/\./p/')
    if [[ -z "$first" ]]; then
        out=$(SC_FACTOR="$sc" RESULTS_DIR="results_3d_sc$tag" sbatch run_grid_3d.slurm "$CASE")
        first=$(echo "$out" | awk '{print $NF}')
        echo "SC_FACTOR=$sc  builds geom cache  -> job $first  -> results_3d_sc$tag"
    else
        SC_FACTOR="$sc" RESULTS_DIR="results_3d_sc$tag" \
            sbatch --dependency=afterok:"$first" run_grid_3d.slurm "$CASE"
        echo "SC_FACTOR=$sc  after job $first  -> results_3d_sc$tag"
    fi
done
