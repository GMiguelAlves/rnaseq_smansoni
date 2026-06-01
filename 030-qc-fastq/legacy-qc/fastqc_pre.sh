#!/bin/bash
#SBATCH --job-name=fastqc_pre
#SBATCH --output=logs/qc_pre/fastqc_pre_%A_%a.out
#SBATCH --error=logs/qc_pre/fastqc_pre_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --array=1-180%5

set -euo pipefail

# =====================
# Ativar conda
# =====================
source /home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate
conda activate rna-tools

# =====================
# Caminhos
# =====================
PROJECT_DIR="/home/${USER}@bio.ib.unicamp.br/rnaseq_smansoni"
SCRATCH_DIR="/scratch/Schisto-epigenetics/gustavo"
RAW_DIR="${SCRATCH_DIR}/fastq_ftp"
QC_PRE_DIR="${SCRATCH_DIR}/fastqc_pre"
META="${PROJECT_DIR}/010-reference/RNAseq_metadata.tsv"

mkdir -p "$QC_PRE_DIR" logs/qc_pre


if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
    echo "[ERRO] Este script deve ser executado como um job array."
    exit 1
fi

line_num=$((SLURM_ARRAY_TASK_ID + 1))


sample_name=$(awk -F '\t' -v line="$line_num" 'NR==line {print $3}' "$META")
runs=$(awk -F '\t' -v line="$line_num" 'NR==line {print $8","$13}' "$META")

if [ -z "$sample_name" ] || [ -z "$runs" ]; then
    echo "[INFO] Nenhuma amostra encontrada para linha $line_num"
    exit 0
fi

IFS=',' read -r -a run_ids <<< "$runs"


find_fastq() {
    local run_id="$1"
    local read_end="$2"
    local file="${RAW_DIR}/${run_id}_${read_end}.fastq.gz"
    [[ -f "$file" ]] && echo "$file" || echo ""
}

for run in "${run_ids[@]}"; do
    R1=$(find_fastq "$run" "1")
    R2=$(find_fastq "$run" "2")

    if [[ -n "$R1" && -n "$R2" ]]; then
        echo "[INFO] Processando ${sample_name} - ${run}"
        fastqc "$R1" "$R2" --outdir "$QC_PRE_DIR" --threads "${SLURM_CPUS_PER_TASK}"
    else
        echo "[ERRO] Arquivos não encontrados para run ${run} (${sample_name})"
        exit 1
    fi
done

echo "[OK] FastQC concluído para ${sample_name}"
