#!/bin/bash
#SBATCH --job-name=fastqc_post
#SBATCH --output=logs/qc_post/fastqc_post_%A_%a.out
#SBATCH --error=logs/qc_post/fastqc_post_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --array=1-75%5

set -euo pipefail

source /home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate
conda activate rna-tools

PROJECT_DIR="/home/${USER}@bio.ib.unicamp.br/rnaseq_smansoni"
SCRATCH_DIR="/scratch/Schisto-epigenetics/gustavo"
TRIM_DIR="${SCRATCH_DIR}/trimmed"
QC_POST_DIR="${SCRATCH_DIR}/fastqc_post"
META="${PROJECT_DIR}/010-reference/RNAseq_metadata.tsv"

mkdir -p "$QC_POST_DIR" logs/qc_post

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
    echo "[ERRO] Este script deve ser executado como um job array."
    exit 1
fi

line_num=$((SLURM_ARRAY_TASK_ID + 1))

sample_name=$(awk -F '\t' -v line="$line_num" 'NR==line {print $3}' "$META")

if [ -z "$sample_name" ]; then
    echo "[INFO] Nenhuma amostra encontrada para linha $line_num"
    exit 0
fi

sample_name_base=$(echo "$sample_name" | sed -E 's/_R[0-9]+$//')
rep=$(echo "$sample_name" | sed -E 's/.*_R([0-9]+)$/\1/')

R1="${TRIM_DIR}/${sample_name_base}_R${rep}_R1_trimmed.fastq.gz"
R2="${TRIM_DIR}/${sample_name_base}_R${rep}_R2_trimmed.fastq.gz"

if [[ -f "$R1" && -f "$R2" ]]; then
    echo "[INFO] Processando ${sample_name_base} (rep${rep})"
    echo "[INFO] R1: $R1"
    echo "[INFO] R2: $R2"

    fastqc "$R1" "$R2" --outdir "$QC_POST_DIR" --threads "${SLURM_CPUS_PER_TASK}"

    echo "[OK] FastQC concluído para ${sample_name_base} (rep${rep})"
else
    echo "[ERRO] Arquivos trimmados não encontrados para ${sample_name_base} (rep${rep})"
    exit 1
fi
