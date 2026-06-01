#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --output=logs/multi/multiqc_%j.out
#SBATCH --error=logs/multi/multiqc_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

set -euo pipefail


source /home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate
conda activate rna-tools


PROJECT_DIR="/home/${USER}@bio.ib.unicamp.br/rnaseq_smansoni"
SCRATCH_DIR="/scratch/Schisto-epigenetics/gustavo"

FILE_DIR="${PROJECT_DIR}/030-qc-fastq"

QC_PRE_DIR="${SCRATCH_DIR}/fastqc_pre"
QC_POST_DIR="${SCRATCH_DIR}/fastqc_post"

OUT_PRE="${FILE_DIR}/multiqc_report_pre"
OUT_POST="${FILE_DIR}/multiqc_report_post"
OUT_COMBINED="${FILE_DIR}/multiqc_report_all"

mkdir -p "$OUT_PRE" "$OUT_POST" "$OUT_COMBINED" logs/multi


echo "[INFO] Gerando MultiQC para fastqc_pre..."
multiqc "${QC_PRE_DIR}" -o "$OUT_PRE" --interactive

echo "[INFO] Gerando MultiQC para fastqc_post..."
multiqc "${QC_POST_DIR}" -o "$OUT_POST" --interactive

echo "[INFO] Gerando MultiQC combinado (pré + pós)..."
multiqc "${QC_PRE_DIR}" "${QC_POST_DIR}" -o "$OUT_COMBINED" --interactive

echo "[OK] Relatórios MultiQC gerados:"
echo " - Pré-trimming:   ${OUT_PRE}/multiqc_report.html"
echo " - Pós-trimming:   ${OUT_POST}/multiqc_report.html"
echo " - Combinado:      ${OUT_COMBINED}/multiqc_report.html"
