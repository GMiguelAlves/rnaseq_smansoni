#!/bin/bash
#SBATCH --job-name=fastq_download
#SBATCH --output=logs/fastq_download_%j.log
#SBATCH --error=logs/fastq_download_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=5:00:00

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Uso: $0 <PROJECT_ID>"
    exit 1
fi

PROJECT_ID=$1
BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"

DATASET_DIR="${BASE_DIR}/datasets/${PROJECT_ID}"
CONFIG="${DATASET_DIR}/config.yaml"
OUTDIR=""

mkdir -p "$OUTDIR" "${BASE_DIR}/logs"

if [ ! -f "$CONFIG" ]; then
    echo "Config não encontrado: $CONFIG"
    exit 1
fi

THREADS=$(grep "threads:" "$CONFIG" | awk '{print $2}')
LINK_FILE=$(grep "link_file:" "$CONFIG" | awk '{print $2}')

LINK_PATH="${DATASET_DIR}/${LINK_FILE}"

if [ ! -f "$LINK_PATH" ]; then
    echo "Arquivo de links não encontrado: $LINK_PATH"
    exit 1
fi

cd "$OUTDIR"

TMP_LINKS=$(mktemp)
grep -Eo '(ftp|https)://[^ ]+' "$LINK_PATH" > "$TMP_LINKS"

xargs -n 1 -P "$THREADS" wget -c < "$TMP_LINKS"

rm "$TMP_LINKS"

echo "Download concluído: $PROJECT_ID"
