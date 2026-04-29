#!/bin/bash
#SBATCH --job-name=star_index_gtf
#SBATCH --output=logs/star_index_gtf_%j.out
#SBATCH --error=logs/star_index_gtf_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=08:00:00

set -euo pipefail

source "$(dirname "$0")/../config/pipeline_config.sh"
activate_conda_env "$RNA_TOOLS_ENV"

	REF_DATA="${DATA_DIR}"

FA="${REF_GENOME_FA}"
GFF3="${REF_GFF3}"
GTF="${REF_GTF}"


mkdir -p "${REF_DIR}/logs" "${STAR_INDEX_DIR}"

if [[ -f "$GFF3" && ! -f "$GTF" ]]; then
  echo "[INFO] Convertendo GFF3 para GTF..."
  gffread "$GFF3" -T -o "$GTF"
fi

echo "[INFO] Criando índice STAR com GTF..."
STAR --runMode genomeGenerate \
  --runThreadN ${SLURM_CPUS_PER_TASK} \
  --genomeDir "${STAR_INDEX_DIR}" \
  --genomeFastaFiles "${FA}" \
  --sjdbGTFfile "${GTF}" \
  --genomeSAindexNbases 10 \
  --limitGenomeGenerateRAM 170000000000

echo "[OK] Índice STAR criado em ${STAR_INDEX_DIR}"
