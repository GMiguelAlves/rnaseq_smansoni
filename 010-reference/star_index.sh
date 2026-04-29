#!/bin/bash
#SBATCH --job-name=ref_index
#SBATCH --output=logs/indexstar_%j.out
#SBATCH --error=logs/indexstar_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00

set -euo pipefail

source "$(dirname "$0")/../config/pipeline_config.sh"
activate_conda_env "$RNA_TOOLS_ENV"
REF_URL="$GENOME_URL"
REF_DATA="$DATA_DIR"
REF_FA="$REF_GENOME_FA"

mkdir -p "$REF_DATA" "$STAR_INDEX_DIR"

if ! command -v STAR &> /dev/null; then
    echo "[ERRO] STAR não encontrado no ambiente Conda!"
    echo "Instale com: conda install -c bioconda star"
    exit 1
fi

echo "[INFO] STAR encontrado em: $(command -v STAR)"
echo "[INFO] Versão do STAR: $(STAR --version)"

########################################
# 1. Baixar genoma se necessário
########################################
if [[ ! -f "$REF_FA" ]]; then
    echo "[INFO] Baixando genoma de $REF_URL"
    wget -O "${REF_FA}.gz" "$REF_URL"
    echo "[INFO] Descompactando genoma"
    gunzip -f "${REF_FA}.gz"
else
    echo "[INFO] Genoma já presente em $REF_FA"
fi

########################################
# 2. Criar índice STAR
########################################
if [[ -z "$(ls -A "$STAR_INDEX_DIR")" ]]; then
    echo "[INFO] Criando índice para STAR..."
    STAR --runMode genomeGenerate \
         --genomeDir "$STAR_INDEX_DIR" \
         --genomeFastaFiles "$REF_FA" \
         --runThreadN "$SLURM_CPUS_PER_TASK" \
         --genomeSAindexNbases 12
else
    echo "[INFO] Índice STAR já existe em $STAR_INDEX_DIR"
fi

echo "[SUCESSO] Índices criados em:"
echo " - STAR: $STAR_INDEX_DIR"
