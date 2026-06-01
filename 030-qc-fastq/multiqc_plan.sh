#!/bin/bash
#SBATCH --job-name=multiqc_030
#SBATCH --output=logs/multiqc/multiqc_%j.out
#SBATCH --error=logs/multiqc/multiqc_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Uso: sbatch $0 <PROJECT> <SCRATCH_ROOT>"
    echo "Exemplo: sbatch $0 PRJEB14695 /scratch/Schisto-epigenetics/gustavo"
    exit 1
fi

PROJECT=$1
SCRATCH_ROOT=$2
PROJECT_SCRATCH="${SCRATCH_ROOT}/${PROJECT}"
OUT_DIR="${PROJECT_SCRATCH}/multiqc_030"

if [ -z "${PROJECT_DIR:-}" ]; then
    if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "${SLURM_SUBMIT_DIR}/../config/pipeline_config.sh" ]; then
        PROJECT_DIR="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
    else
        PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
    fi
fi

source "${PROJECT_DIR}/config/pipeline_config.sh"
activate_rna_tools

mkdir -p "$OUT_DIR" logs/multiqc

multiqc \
    "${PROJECT_SCRATCH}/fastqc_raw" \
    "${PROJECT_SCRATCH}/fastqc_trimmed_runs" \
    "${PROJECT_SCRATCH}/fastqc_merged" \
    -o "$OUT_DIR" \
    -n "${PROJECT}_multiqc_030.html"

echo "[OK] MultiQC concluido: ${OUT_DIR}/${PROJECT}_multiqc_030.html"
