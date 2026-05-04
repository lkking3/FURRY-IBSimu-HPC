#!/usr/bin/env bash
# =============================================================================
# orchestrate_curved_grid_sweep.slurm
#
# SLURM job array for sweeping over curved-grid parameters.
#
# Usage (submit from the curved_grid/ directory):
#   sbatch \
#     --array=0-<N-1> \
#     --export=SWEEP_MANIFEST=<path/to/sweep_manifest.json>,\
#              SWEEP_CONFIG=<path/to/curved_grid_case.py>,\
#              SWEEP_BINARY=<path/to/multi_grid_2d_curved> \
#     orchestrate_curved_grid_sweep.slurm
#
# The pipeline generates these values when you run:
#   python curved_grid_pipeline.py --sweep --config <case_file.py>
#
# It prints the exact sbatch command with the correct --array=0-<N-1> bound
# and --export values, so you can copy-paste it directly.
#
# Environment variables
# ----------------------
# SWEEP_MANIFEST   Path to the JSON manifest written by --sweep (required)
# SWEEP_CONFIG     Path to the case .py file that defines make_case + SWEEP (required)
# SWEEP_BINARY     Path to the compiled binary (default: ../sim/multi_grid_2d_curved)
# REDUCE_SCRIPT    Path to reduce_best.py (default: ../tools/reduce_best.py)
#
# Tunable SLURM directives — edit below or override with sbatch --<flag>=<val>
# =============================================================================
#SBATCH --job-name=cg_sweep
#SBATCH --output=logs/sweep_%A_%a.out
#SBATCH --error=logs/sweep_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
# Uncomment and edit for GPU/partition requirements:
##SBATCH --partition=compute
##SBATCH --gres=gpu:1

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
SIM_DIR="$(realpath "$PIPELINE_DIR/../sim")"
TOOLS_DIR="$(realpath "$PIPELINE_DIR/../tools")"

PIPELINE_SCRIPT="$PIPELINE_DIR/curved_grid_pipeline.py"
BINARY="${SWEEP_BINARY:-$SIM_DIR/multi_grid_2d_curved}"
REDUCE="${REDUCE_SCRIPT:-$TOOLS_DIR/reduce_best.py}"

# ---------------------------------------------------------------------------
# Validate required inputs
# ---------------------------------------------------------------------------
if [[ -z "$SWEEP_MANIFEST" ]]; then
    echo "ERROR: SWEEP_MANIFEST is not set." >&2
    exit 1
fi
if [[ -z "$SWEEP_CONFIG" ]]; then
    echo "ERROR: SWEEP_CONFIG is not set." >&2
    exit 1
fi
if [[ ! -f "$SWEEP_MANIFEST" ]]; then
    echo "ERROR: SWEEP_MANIFEST not found: $SWEEP_MANIFEST" >&2
    exit 1
fi
if [[ ! -f "$SWEEP_CONFIG" ]]; then
    echo "ERROR: SWEEP_CONFIG not found: $SWEEP_CONFIG" >&2
    exit 1
fi
if [[ ! -f "$BINARY" ]]; then
    echo "ERROR: solver binary not found: $BINARY" >&2
    exit 1
fi
if [[ ! -f "$PIPELINE_SCRIPT" ]]; then
    echo "ERROR: pipeline script not found: $PIPELINE_SCRIPT" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Create logs directory if needed
# ---------------------------------------------------------------------------
mkdir -p "$PIPELINE_DIR/logs"

# ---------------------------------------------------------------------------
# Run this task's simulation
# ---------------------------------------------------------------------------
echo "============================================================"
echo "  Sweep task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}"
echo "  Array job ${SLURM_ARRAY_JOB_ID}"
echo "  Manifest:  $SWEEP_MANIFEST"
echo "  Config:    $SWEEP_CONFIG"
echo "  Binary:    $BINARY"
echo "  Node:      $(hostname)"
echo "  Started:   $(date)"
echo "============================================================"

python "$PIPELINE_SCRIPT" \
    --sweep-task "$SLURM_ARRAY_TASK_ID" \
    --manifest   "$SWEEP_MANIFEST" \
    --config     "$SWEEP_CONFIG" \
    --binary     "$BINARY"

TASK_EXIT=$?

echo "Task ${SLURM_ARRAY_TASK_ID} finished at $(date) with exit code ${TASK_EXIT}"
exit $TASK_EXIT
