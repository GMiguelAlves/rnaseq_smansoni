#!/bin/bash
#SBATCH --job-name=trim_galore
#SBATCH --output=logs/trim/trimming_%A_%a.out
#SBATCH --error=logs/trim/trimming_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=18:00:00
#SBATCH --array=1-75%5

set -euo pipefail

# =====================
# Ativar conda
# =====================
source "/home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate"
conda activate rna-tools

# =====================
# Caminhos
# =====================
PROJECT_DIR="/home/${USER}@bio.ib.unicamp.br/rnaseq_smansoni"
RAW_DIR="/scratch/Schisto-epigenetics/gustavo/fastq_ftp"
TRIM_DIR="/scratch/Schisto-epigenetics/gustavo/trimmed"
META="${PROJECT_DIR}/010-reference/RNAseq_metadata.tsv"

mkdir -p "$TRIM_DIR" logs/trim

if [ -z "${SLURM_ARRAY_TASK_ID}" ]; then
    echo "[ERRO] Este script deve ser executado como um job array."
    exit 1
fi

line_num=$((SLURM_ARRAY_TASK_ID + 1))

# Extrair SampleName e runs
sample_name=$(awk -F '\t' -v line="$line_num" 'NR==line {print $3}' "$META")
runs=$(awk -F '\t' -v line="$line_num" 'NR==line {print $8","$13}' "$META")

if [ -z "$sample_name" ] || [ -z "$runs" ]; then
    echo "[INFO] Nenhuma amostra encontrada para linha $line_num"
    exit 0
fi

IFS=',' read -r -a run_ids <<< "$runs"
n_runs=${#run_ids[@]}
threads_per_run=$(( SLURM_CPUS_PER_TASK / n_runs ))
[[ $threads_per_run -lt 1 ]] && threads_per_run=1

echo "[INFO] Amostra: $sample_name"
echo "[INFO] Runs: ${run_ids[*]}"
echo "[INFO] Usando $threads_per_run threads por run"

# Arrays para guardar outputs
trimmed_r1_files=()
trimmed_r2_files=()

# Processar cada run individualmente
for run in "${run_ids[@]}"; do
    r1_file="${RAW_DIR}/${run}_1.fastq.gz"
    r2_file="${RAW_DIR}/${run}_2.fastq.gz"

    if [[ ! -f "$r1_file" || ! -f "$r2_file" ]]; then
        echo "[ERRO] Arquivos não encontrados para $run"
        exit 1
    fi

    echo "[INFO] Rodando Trim Galore para $run"
    trim_galore --paired \
        --quality 20 \
        --length 20 \
        --cores "$threads_per_run" \
        --output_dir "$TRIM_DIR" \
        "$r1_file" "$r2_file"

    # Renomear para evitar confusão
    mv "${TRIM_DIR}/${run}_1_val_1.fq.gz" "${TRIM_DIR}/${run}_R1_trimmed.fastq.gz"
    mv "${TRIM_DIR}/${run}_2_val_2.fq.gz" "${TRIM_DIR}/${run}_R2_trimmed.fastq.gz"

    trimmed_r1_files+=("${TRIM_DIR}/${run}_R1_trimmed.fastq.gz")
    trimmed_r2_files+=("${TRIM_DIR}/${run}_R2_trimmed.fastq.gz")
done

# Concatenar arquivos já trimmados
echo "[INFO] Concatenando runs para $sample_name"
zcat "${trimmed_r1_files[@]}" | gzip > "${TRIM_DIR}/${sample_name}_R1_trimmed.fastq.gz"
zcat "${trimmed_r2_files[@]}" | gzip > "${TRIM_DIR}/${sample_name}_R2_trimmed.fastq.gz"

# Remover arquivos intermediários
for file in "${trimmed_r1_files[@]}" "${trimmed_r2_files[@]}"; do
    rm -f "$file"
done

echo "[OK] Trimming concluído para $sample_name"
