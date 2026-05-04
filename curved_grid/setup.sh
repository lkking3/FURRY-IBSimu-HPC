#!/bin/bash
# setup.sh
#
# One-command setup for the two_grid_curved simulation branch.
# Run this once after transferring the folder to the HPC to build the solver
# and verify the Python pipeline before submitting jobs.
#
# Usage:
#   bash setup.sh              # build + verify only
#   bash setup.sh --submit     # build + verify + submit the example job
#   bash setup.sh --dry-run    # print geometry table only, no build
#   DEBUG=1 bash setup.sh      # enable bash trace

set -euo pipefail
[ "${DEBUG:-0}" = "1" ] && set -x
trap 'echo "[setup] ERROR at line $LINENO: $BASH_COMMAND" >&2' ERR

# ---- resolve script directory (works when called from any cwd) ----
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ---- parse args ----
DO_SUBMIT=0
DRY_RUN_ONLY=0
for arg in "$@"; do
    case "$arg" in
        --submit)       DO_SUBMIT=1 ;;
        --dry-run)      DRY_RUN_ONLY=1 ;;
        *) echo "[setup] unknown argument: $arg" >&2; exit 1 ;;
    esac
done

# ---- LMOD init (Omega) ----
[ -f /etc/profile.d/modules.sh ]          && . /etc/profile.d/modules.sh          || true
[ -f /usr/share/lmod/lmod/init/bash ]     && . /usr/share/lmod/lmod/init/bash     || true

# ============================================================
# 1. Verify required source files are present
# ============================================================
echo "[setup] checking source files..."
REQUIRED_FILES=(
    multi_grid_2d.cpp
    multi_grid_pipeline.py
    curved_grid_pipeline.py
    curved_grid_case_example.py
    Makefile
    orchestrate_curved_grid_pipeline.slurm
)
MISSING=0
for f in "${REQUIRED_FILES[@]}"; do
    if [[ ! -f "$SCRIPT_DIR/$f" ]]; then
        echo "[setup] MISSING: $f" >&2
        MISSING=1
    fi
done
if [[ "$MISSING" -eq 1 ]]; then
    echo "[setup] one or more required files are missing — aborting." >&2
    exit 1
fi
echo "[setup] all source files present."

# ============================================================
# 2. Python pipeline dry-run (geometry verification, no binary needed)
# ============================================================
echo "[setup] running Python pipeline dry-run..."
python3 "$SCRIPT_DIR/curved_grid_pipeline.py" \
    --config "$SCRIPT_DIR/curved_grid_case_example.py" \
    --dry-run
echo "[setup] Python pipeline OK."

if [[ "$DRY_RUN_ONLY" -eq 1 ]]; then
    echo "[setup] --dry-run only, stopping here."
    exit 0
fi

# ============================================================
# 3. Build the C++ solver
# ============================================================
echo "[setup] building multi_grid_2d..."
make -C "$SCRIPT_DIR" multi_grid_2d

if [[ ! -x "$SCRIPT_DIR/multi_grid_2d" ]]; then
    echo "[setup] build failed — multi_grid_2d not found after make." >&2
    exit 1
fi
echo "[setup] build successful: $SCRIPT_DIR/multi_grid_2d"

# ============================================================
# 4. Smoke-test: binary responds to --help or exits cleanly on bad env
# ============================================================
echo "[setup] verifying binary is executable..."
# Run with no GRID_STACK so it exits with usage error (not a crash).
# We just confirm it at least launches and links correctly.
if GRID_STACK="" "$SCRIPT_DIR/multi_grid_2d" 2>&1 | grep -qiE 'error|usage|ibsimu|grid'; then
    echo "[setup] binary responds correctly."
else
    echo "[setup] WARNING: binary smoke-test produced unexpected output — check manually." >&2
fi

# ============================================================
# 5. Summary
# ============================================================
echo ""
echo "============================================"
echo " two_grid_curved setup complete"
echo "  directory : $SCRIPT_DIR"
echo "  binary    : $SCRIPT_DIR/multi_grid_2d"
echo "  config    : $SCRIPT_DIR/curved_grid_case_example.py"
echo "============================================"
echo ""
echo " To submit the example job:"
echo "   sbatch $SCRIPT_DIR/orchestrate_curved_grid_pipeline.slurm"
echo ""
echo " To use a custom config:"
echo "   CURVED_GRID_CONFIG=/path/to/my_case.py \\"
echo "     sbatch $SCRIPT_DIR/orchestrate_curved_grid_pipeline.slurm"
echo ""

# ============================================================
# 6. Optional immediate submission
# ============================================================
if [[ "$DO_SUBMIT" -eq 1 ]]; then
    echo "[setup] submitting job..."
    sbatch "$SCRIPT_DIR/orchestrate_curved_grid_pipeline.slurm"
fi
