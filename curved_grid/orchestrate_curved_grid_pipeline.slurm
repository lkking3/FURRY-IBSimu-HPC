#!/bin/bash
# ============================================================
# IBSimu curved-grid pipeline — SLURM job script
# ============================================================
# Submit from the curved_grid/ directory:
#   cd curved_grid && sbatch orchestrate_curved_grid_pipeline.slurm
#
# All defaults are resolved relative to that directory, so the
# correct binary, config, and pipeline script are picked up
# automatically — no env vars needed for a standard run.
#
# Override any #SBATCH line on the sbatch command line, e.g.:
#   sbatch --mem=8GB orchestrate_curved_grid_pipeline.slurm
#   sbatch --mem=16GB --time=02:00:00 orchestrate_curved_grid_pipeline.slurm
#   sbatch --partition=week --mem=8GB orchestrate_curved_grid_pipeline.slurm
#
# Override paths if needed:
#   CURVED_GRID_CONFIG=/path/to/my_case.py \
#     sbatch orchestrate_curved_grid_pipeline.slurm
# ============================================================

#SBATCH -J ibsimu_curved_grid
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH -t 06:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err

set -euo pipefail

# Prevent accidental OpenMP/MKL threading inside IBSimu (uses all cores
# if not constrained, without any actual speedup from our serial code).
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# IBSimu's GeomPlotter does not need an X display to write PNG files,
# but setting DISPLAY= avoids spurious Xlib errors on headless nodes.
export DISPLAY=

# ---- LMOD init (Omega) ----
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh || true
[ -f /usr/share/lmod/lmod/init/bash ] && . /usr/share/lmod/lmod/init/bash || true

# ---- Load required modules (edit to match your environment) ----
# module load ibsimu/1.0.6dev   # uncomment if IBSimu is a loadable module
# module load gcc/12             # if a specific compiler runtime is needed

# SLURM_SUBMIT_DIR is set to wherever sbatch was called from.
# Submit from curved_grid/ so all relative paths resolve correctly.
PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
SIM_DIR="$(realpath "$PIPELINE_DIR/../sim")"

# Defaults: curved-grid binary, case file, and pipeline — all within this sector.
CONFIG_PATH="${CURVED_GRID_CONFIG:-$PIPELINE_DIR/curved_grid_case_example.py}"
BINARY_PATH="${CURVED_GRID_BINARY:-$SIM_DIR/multi_grid_2d_curved}"
PIPELINE_SCRIPT="${CURVED_GRID_PIPELINE:-$PIPELINE_DIR/curved_grid_pipeline.py}"

if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "[curved-grid] missing config: $CONFIG_PATH" >&2
  exit 2
fi

if [[ ! -x "$BINARY_PATH" ]]; then
  echo "[curved-grid] missing solver binary: $BINARY_PATH" >&2
  echo "[curved-grid] build it first with: cd ../sim && make multi_grid_2d_curved" >&2
  exit 2
fi

if [[ ! -f "$PIPELINE_SCRIPT" ]]; then
  echo "[curved-grid] missing pipeline script: $PIPELINE_SCRIPT" >&2
  exit 2
fi

echo "[curved-grid] pipeline : $PIPELINE_SCRIPT"
echo "[curved-grid] config   : $CONFIG_PATH"
echo "[curved-grid] binary   : $BINARY_PATH"
echo "[curved-grid] CPUs     : ${SLURM_CPUS_PER_TASK:-1}  mem: ${SLURM_MEM_PER_NODE:-?} MB"

python3 "$PIPELINE_SCRIPT" \
  --config "$CONFIG_PATH" \
  --binary "$BINARY_PATH"
