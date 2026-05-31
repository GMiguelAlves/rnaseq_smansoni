#!/bin/bash
#SBATCH --job-name=merge_samples
#SBATCH --output=logs/merge_samples/merge_samples_%A_%a.out
#SBATCH --error=logs/merge_samples/merge_samples_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Uso: sbatch --array=1-N $0 <QC_PLAN.csv>"
    exit 1
fi

PLAN=$1

if [ -z "${PROJECT_DIR:-}" ]; then
    if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "${SLURM_SUBMIT_DIR}/../config/pipeline_config.sh" ]; then
        PROJECT_DIR="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
    else
        PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
    fi
fi

source "${PROJECT_DIR}/config/pipeline_config.sh"
activate_python_env

SCRIPT_DIR="${SCRIPT_DIR:-${PROJECT_DIR}/030-qc-fastq}"

mkdir -p logs/merge_samples

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "[ERRO] Execute como job array: sbatch --array=1-N ..."
    exit 1
fi

SAMPLE_ID=$(
    python -c "import csv,sys; samples=sorted({r['sample_id'] for r in csv.DictReader(open(sys.argv[1], newline=''))}); print(samples[int(sys.argv[2])-1])" \
        "$PLAN" "$SLURM_ARRAY_TASK_ID"
)

python "${SCRIPT_DIR}/merge_sample_from_plan.py" \
    --plan "$PLAN" \
    --sample-id "$SAMPLE_ID"

