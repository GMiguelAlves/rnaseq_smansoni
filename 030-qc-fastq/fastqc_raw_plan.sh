#!/bin/bash
#SBATCH --job-name=fastqc_raw
#SBATCH --output=logs/qc_raw/fastqc_raw_%A_%a.out
#SBATCH --error=logs/qc_raw/fastqc_raw_%A_%a.err
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

mkdir -p "$OUT_DIR" logs/qc_raw

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "[ERRO] Execute como job array: sbatch --array=1-N ..."
    exit 1
fi

read -r RAW_R1 RAW_R2 SAMPLE_ID RUN_ACCESSION < <(
    python -c "import csv,sys; rows=list(csv.DictReader(open(sys.argv[1], newline=''))); r=rows[int(sys.argv[2])-1]; print(r['raw_r1'], r['raw_r2'], r['sample_id'], r['run_accession'], sep='\t')" \
        "$PLAN" "$SLURM_ARRAY_TASK_ID"
)

if [[ ! -f "$RAW_R1" || ! -f "$RAW_R2" ]]; then
    echo "[ERRO] FASTQs brutos ausentes para ${RUN_ACCESSION}:"
    echo "$RAW_R1"
    echo "$RAW_R2"
    exit 1
fi

echo "[INFO] FastQC bruto: ${SAMPLE_ID} ${RUN_ACCESSION}"
fastqc "$RAW_R1" "$RAW_R2" --outdir "$OUT_DIR" --threads "${SLURM_CPUS_PER_TASK}"
echo "[OK] FastQC bruto concluido: ${RUN_ACCESSION}"
