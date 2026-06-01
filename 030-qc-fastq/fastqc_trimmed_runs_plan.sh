#!/bin/bash
#SBATCH --job-name=fastqc_trim_runs
#SBATCH --output=logs/qc_trimmed_runs/fastqc_trimmed_runs_%A_%a.out
#SBATCH --error=logs/qc_trimmed_runs/fastqc_trimmed_runs_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Uso: sbatch --array=1-N $0 <QC_PLAN.csv> <OUTPUT_DIR>"
    exit 1
fi

PLAN=$1
OUT_DIR=$2

if [ -z "${PROJECT_DIR:-}" ]; then
    if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "${SLURM_SUBMIT_DIR}/../config/pipeline_config.sh" ]; then
        PROJECT_DIR="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
    else
        PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
    fi
fi

source "${PROJECT_DIR}/config/pipeline_config.sh"
activate_rna_tools

mkdir -p "$OUT_DIR" logs/qc_trimmed_runs

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "[ERRO] Execute como job array: sbatch --array=1-N ..."
    exit 1
fi

read -r TRIM_R1 TRIM_R2 SAMPLE_ID RUN_ACCESSION < <(
    python -c "import csv,sys; rows=list(csv.DictReader(open(sys.argv[1], newline=''))); r=rows[int(sys.argv[2])-1]; print(r['trimmed_run_r1'], r['trimmed_run_r2'], r['sample_id'], r['run_accession'], sep='\t')" \
        "$PLAN" "$SLURM_ARRAY_TASK_ID"
)

if [[ ! -f "$TRIM_R1" || ! -f "$TRIM_R2" ]]; then
    echo "[ERRO] FASTQs trimmados ausentes para ${RUN_ACCESSION}"
    exit 1
fi

echo "[INFO] FastQC run trimmado: ${SAMPLE_ID} ${RUN_ACCESSION}"
fastqc "$TRIM_R1" "$TRIM_R2" --outdir "$OUT_DIR" --threads "${SLURM_CPUS_PER_TASK}"
echo "[OK] FastQC run trimmado concluido: ${RUN_ACCESSION}"
