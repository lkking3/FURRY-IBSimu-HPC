#!/bin/bash
#SBATCH -J ibsimu_multi_grid
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH -t 00:30:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err

set -euo pipefail

CONFIG_PATH="${MULTI_GRID_CONFIG:-$PWD/multi_grid_case_example.py}"
BINARY_PATH="${MULTI_GRID_BINARY:-$PWD/multi_grid_2d}"

if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "[multi-grid] missing config: $CONFIG_PATH" >&2
  exit 2
fi

if [[ ! -x "$BINARY_PATH" ]]; then
  echo "[multi-grid] missing solver binary: $BINARY_PATH" >&2
  echo "[multi-grid] build it first with: make multi_grid_2d" >&2
  exit 2
fi

python3 "$PWD/multi_grid_pipeline.py" --config "$CONFIG_PATH" --binary "$BINARY_PATH"
