#!/bin/bash
# ============================================================
# IBSimu multi-grid pipeline — SLURM job script
# ============================================================
# Submit from the Multi_Grid/ directory:
#   cd Multi_Grid && sbatch orchestrate_multi_grid_pipeline.slurm
#
# All defaults are resolved relative to that directory, so the
# correct binary, config, and pipeline script are picked up
# automatically — no env vars needed for a standard run.
#
# Override any #SBATCH line on the sbatch command line, e.g.:
#   sbatch --mem=4GB orchestrate_multi_grid_pipeline.slurm
#   sbatch --time=01:00:00 orchestrate_multi_grid_pipeline.slurm
#
# Override paths if needed:
#   MULTI_GRID_CONFIG=/path/to/my_case.py \
#     sbatch orchestrate_multi_grid_pipeline.slurm
# ============================================================

#SBATCH -J ibsimu_multi_grid
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH -t 00:30:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err

set -euo pipefail

# Prevent accidental OpenMP/MKL threading inside IBSimu.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

export DISPLAY=

# ---- LMOD init (Omega) ----
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh || true
[ -f /usr/share/lmod/lmod/init/bash ] && . /usr/share/lmod/lmod/init/bash || true

# ---- Load required modules (edit to match your environment) ----
# module load ibsimu/1.0.6dev
# module load gcc/12

# SLURM_SUBMIT_DIR is set to wherever sbatch was called from.
# Submit from Multi_Grid/ so all relative paths resolve correctly.
PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
SIM_DIR="$(realpath "$PIPELINE_DIR/../sim")"

# Defaults: flat multi-grid binary, case file, and pipeline — all within this sector.
CONFIG_PATH="${MULTI_GRID_CONFIG:-$PIPELINE_DIR/multi_grid_case_example.py}"
BINARY_PATH="${MULTI_GRID_BINARY:-$SIM_DIR/multi_grid_2d}"
PIPELINE_SCRIPT="${MULTI_GRID_PIPELINE:-$PIPELINE_DIR/multi_grid_pipeline.py}"

if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "[multi-grid] missing config: $CONFIG_PATH" >&2
  exit 2
fi

if [[ ! -x "$BINARY_PATH" ]]; then
  echo "[multi-grid] missing solver binary: $BINARY_PATH" >&2
  echo "[multi-grid] build it first with: cd ../sim && make multi_grid_2d" >&2
  exit 2
fi

if [[ ! -f "$PIPELINE_SCRIPT" ]]; then
  echo "[multi-grid] missing pipeline script: $PIPELINE_SCRIPT" >&2
  exit 2
fi

echo "[multi-grid] pipeline : $PIPELINE_SCRIPT"
echo "[multi-grid] config   : $CONFIG_PATH"
echo "[multi-grid] binary   : $BINARY_PATH"
echo "[multi-grid] CPUs     : ${SLURM_CPUS_PER_TASK:-1}  mem: ${SLURM_MEM_PER_NODE:-?} MB"

python3 "$PIPELINE_SCRIPT" \
  --config "$CONFIG_PATH" \
  --binary "$BINARY_PATH"
