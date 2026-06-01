#!/bin/bash
#SBATCH --job-name=trim_runs
#SBATCH --output=logs/trim_runs/trim_runs_%A_%a.out
#SBATCH --error=logs/trim_runs/trim_runs_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=08:00:00

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
activate_rna_tools

mkdir -p logs/trim_runs

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "[ERRO] Execute como job array: sbatch --array=1-N ..."
    exit 1
fi

read -r RAW_R1 RAW_R2 TRIM_R1 TRIM_R2 SAMPLE_ID RUN_ACCESSION < <(
    python -c "import csv,sys; rows=list(csv.DictReader(open(sys.argv[1], newline=''))); r=rows[int(sys.argv[2])-1]; print(r['raw_r1'], r['raw_r2'], r['trimmed_run_r1'], r['trimmed_run_r2'], r['sample_id'], r['run_accession'], sep='\t')" \
        "$PLAN" "$SLURM_ARRAY_TASK_ID"
)

if [[ ! -f "$RAW_R1" || ! -f "$RAW_R2" ]]; then
    echo "[ERRO] FASTQs brutos ausentes para ${RUN_ACCESSION}"
    exit 1
fi

if [[ -s "$TRIM_R1" && -s "$TRIM_R2" ]]; then
    echo "[SKIP] Trimmed run ja existe: ${RUN_ACCESSION}"
    exit 0
fi

TRIM_DIR=$(dirname "$TRIM_R1")
mkdir -p "$TRIM_DIR"

echo "[INFO] Trim Galore run-level: ${SAMPLE_ID} ${RUN_ACCESSION}"
trim_galore --paired \
    --quality 20 \
    --length 20 \
    --cores "${SLURM_CPUS_PER_TASK}" \
    --output_dir "$TRIM_DIR" \
    "$RAW_R1" "$RAW_R2"

RAW_R1_BASE=$(basename "$RAW_R1")
RAW_R2_BASE=$(basename "$RAW_R2")
TG_R1="${TRIM_DIR}/${RAW_R1_BASE%.fastq.gz}_val_1.fq.gz"
TG_R2="${TRIM_DIR}/${RAW_R2_BASE%.fastq.gz}_val_2.fq.gz"

if [[ ! -f "$TG_R1" || ! -f "$TG_R2" ]]; then
    echo "[ERRO] Saidas esperadas do Trim Galore ausentes:"
    echo "$TG_R1"
    echo "$TG_R2"
    exit 1
fi

mv "$TG_R1" "$TRIM_R1"
mv "$TG_R2" "$TRIM_R2"

echo "[OK] Trimming run-level concluido: ${RUN_ACCESSION}"
