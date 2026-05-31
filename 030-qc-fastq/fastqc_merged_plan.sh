#!/bin/bash
#SBATCH --job-name=fastqc_merged
#SBATCH --output=logs/qc_merged/fastqc_merged_%A_%a.out
#SBATCH --error=logs/qc_merged/fastqc_merged_%A_%a.err
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

mkdir -p "$OUT_DIR" logs/qc_merged

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "[ERRO] Execute como job array: sbatch --array=1-N ..."
    exit 1
fi

read -r SAMPLE_ID MERGED_R1 MERGED_R2 < <(
    python -c "import csv,sys; rows=list(csv.DictReader(open(sys.argv[1], newline=''))); samples=sorted({r['sample_id'] for r in rows}); s=samples[int(sys.argv[2])-1]; r=next(x for x in rows if x['sample_id']==s); print(s, r['merged_sample_r1'], r['merged_sample_r2'], sep='\t')" \
        "$PLAN" "$SLURM_ARRAY_TASK_ID"
)

if [[ ! -f "$MERGED_R1" || ! -f "$MERGED_R2" ]]; then
    echo "[ERRO] FASTQs merged ausentes para ${SAMPLE_ID}"
    echo "$MERGED_R1"
    echo "$MERGED_R2"
    exit 1
fi

echo "[INFO] FastQC pos-merge: ${SAMPLE_ID}"
fastqc "$MERGED_R1" "$MERGED_R2" --outdir "$OUT_DIR" --threads "${SLURM_CPUS_PER_TASK}"
echo "[OK] FastQC pos-merge concluido: ${SAMPLE_ID}"
