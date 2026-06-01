#!/bin/bash
#SBATCH --job-name=salmon_quant
#SBATCH --output=logs/salmon/salmon_%A_%a.out
#SBATCH --error=logs/salmon/salmon_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-75%5

set -euo pipefail

source "/home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate"
conda activate rna-tools

PROJECT_DIR="/home/${USER}@bio.ib.unicamp.br/rnaseq_smansoni"
TRIM_DIR="/scratch/Schisto-epigenetics/gustavo/trimmed"
REF_DIR="${PROJECT_DIR}/010-reference"
INDEX_DIR="${REF_DIR}/salmon_index"
META="${REF_DIR}/RNAseq_metadata.tsv"
OUT_DIR="${PROJECT_DIR}/040-alignment/quants"

mkdir -p "$OUT_DIR" logs/salmon

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
    echo "[ERRO] Este script deve ser executado como job array."
    exit 1
fi

line_num=$((SLURM_ARRAY_TASK_ID + 1))

sample_name=$(awk -F '\t' -v line="$line_num" 'NR==line {print $3}' "$META")

if [ -z "$sample_name" ]; then
    echo "[INFO] Nenhuma amostra encontrada na linha $line_num"
    exit 0
fi

r1="${TRIM_DIR}/${sample_name}_R1_trimmed.fastq.gz"
r2="${TRIM_DIR}/${sample_name}_R2_trimmed.fastq.gz"

if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    echo "[ERRO] FASTQs não encontrados para $sample_name"
    exit 1
fi

echo "[INFO] Rodando Salmon para $sample_name"
salmon quant \
    -i "$INDEX_DIR" \
    -l A \
    -1 "$r1" \
    -2 "$r2" \
    -p $SLURM_CPUS_PER_TASK \
    --validateMappings \
    -o "${OUT_DIR}/${sample_name}"

echo "[OK] Quantificação concluída para $sample_name"
