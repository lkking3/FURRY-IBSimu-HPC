#!/usr/bin/env bash
# =============================================================================
# orchestrate_curved_grid_sweep.slurm
#
# SLURM job array for sweeping over curved-grid parameters.
#
# Usage — just run --sweep from the curved_grid/ directory:
#   python3 curved_grid_pipeline.py --sweep --config <case_file.py>
#
# That command generates the manifest (for provenance) and submits the
# sbatch array job in one step.  No separate manifest generation needed.
#
# Environment variables
# ----------------------
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
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
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

# Resolve SWEEP_CONFIG to an absolute path if it isn't already
if [[ "$SWEEP_CONFIG" != /* ]]; then
    SWEEP_CONFIG="$PIPELINE_DIR/$SWEEP_CONFIG"
fi

# ---------------------------------------------------------------------------
# Validate required inputs
# ---------------------------------------------------------------------------
if [[ -z "$SWEEP_CONFIG" ]]; then
    echo "ERROR: SWEEP_CONFIG is not set." >&2
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
echo "  Config:    $SWEEP_CONFIG"
echo "  Binary:    $BINARY"
echo "  Node:      $(hostname)"
echo "  Started:   $(date)"
echo "============================================================"

python3 "$PIPELINE_SCRIPT" \
    --sweep-task "$SLURM_ARRAY_TASK_ID" \
    --config     "$SWEEP_CONFIG" \
    --binary     "$BINARY"

TASK_EXIT=$?

echo "Task ${SLURM_ARRAY_TASK_ID} finished at $(date) with exit code ${TASK_EXIT}"
exit $TASK_EXIT
