#!/usr/bin/env bash
# =============================================================================
# run_grid_3d.slurm -- single 3-D curved-grid extraction run.
# Submit from curved_grid_3d/:   sbatch run_grid_3d.slurm [case_file.sh]
# IBSimu is single-threaded; request 1 CPU and generous memory for the
# h = 0.2 mm half-domain (~481 M nodes ~ 39 GB; 48 GB gives headroom).
# All compute runs here under SLURM -- never on the login node.
# =============================================================================
#SBATCH --job-name=grid3d
#SBATCH --output=logs/grid3d_%j.out
#SBATCH --error=logs/grid3d_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16          # IBSimu threads the mesh build + solver
#SBATCH --mem=96G                   # ILU0 solve at h=0.2mm needs ~65-95 GB
#SBATCH --time=48:00:00
## #SBATCH --partition=compute      # set to your large-memory partition
## For a coarse h=0.3mm run, 32G is plenty; drop --mem to save queue time.

export ION_THREADS=${SLURM_CPUS_PER_TASK:-16}   # match thread count to allocation

set -euo pipefail
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$SUBMIT_DIR"
mkdir -p logs

CASE_FILE="${1:-case_baseline.sh}"
echo "case file : $CASE_FILE"
# shellcheck disable=SC1090
source "$CASE_FILE"

BIN="$(realpath ../sim/grid_3d_curved)"
[[ -x "$BIN" ]] || { echo "ERROR: build first: (cd ../sim && make -f Makefile.3d)"; exit 1; }

echo "binary    : $BIN"
echo "results   : $RESULTS_DIR"
echo "started   : $(date)"
/usr/bin/time -v "$BIN"
echo "finished  : $(date)"
