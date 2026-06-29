#!/usr/bin/env bash
# Run the meniscus micro-model. Cheap (single bore, axisymmetric).
#   sbatch run_meniscus.slurm [case_meniscus_cell.sh]
#SBATCH --job-name=meniscus
#SBATCH --output=logs/meniscus_%j.out
#SBATCH --error=logs/meniscus_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"; mkdir -p logs
CASE="${1:-case_meniscus_cell.sh}"
# shellcheck disable=SC1090
source "$CASE"
export ION_THREADS=${SLURM_CPUS_PER_TASK:-4}
BIN="$(realpath ../sim/meniscus_cell)"
[[ -x "$BIN" ]] || { echo "build first: (cd ../sim && make -f Makefile.3d meniscus_cell)"; exit 1; }
echo "started $(date)"; /usr/bin/time -v "$BIN"; echo "finished $(date)"
