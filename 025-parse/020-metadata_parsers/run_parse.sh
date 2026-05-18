#!/bin/bash
#SBATCH --job-name=025_parse_02
#SBATCH --output=logs/parse/%x_%A.out
#SBATCH --error=logs/parse/%x_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00

set -euo pipefail

if [ $# -eq 0 ]; then
    echo "Uso: sbatch $0 <PROJECT>"
    exit 1
fi

PROJ=$1

source "/home/${USER}@bio.ib.unicamp.br/miniconda3/etc/profile.d/conda.sh"
conda activate python-list

INTERMEDIATE_DIR="../015-intermediate_folder"

BASE_FILE="${INTERMEDIATE_DIR}/${PROJ}_base.csv"
ENRICHED_FILE="${INTERMEDIATE_DIR}/${PROJ}_enriched.csv"

PARSE_CONFIG="${PROJ}/configs/${PROJ}.yaml"
ENRICH_CONFIG="${PROJ}/configs/${PROJ}_enrich.yaml"

AUTHOR_METADATA="${PROJ}/author_metadata.tsv"

OUTPUT_FILE="Allprojects/${PROJ}_parsed.csv"

mkdir -p "$(dirname "$OUTPUT_FILE")"

if [ ! -f "$BASE_FILE" ]; then
    echo "Erro: base metadata não encontrado:"
    echo "$BASE_FILE"
    exit 1
fi

if [ ! -f "$PARSE_CONFIG" ]; then
    echo "Erro: parse config não encontrado:"
    echo "$PARSE_CONFIG"
    exit 1
fi


INPUT_FILE="$BASE_FILE"

if [ -f "$ENRICH_CONFIG" ]; then
    echo "Enrich config detectado."

    if [ ! -f "$AUTHOR_METADATA" ]; then
        echo "Erro: enrich requerido mas author metadata ausente:"
        echo "$AUTHOR_METADATA"
        exit 1
    fi

    echo "Running enrich..."

    metaqc enrich \
        "$BASE_FILE" \
        "$AUTHOR_METADATA" \
        "$ENRICH_CONFIG" \
        --output "$ENRICHED_FILE"

    INPUT_FILE="$ENRICHED_FILE"
else
    echo "No enrich config detected. Skipping enrich."
fi


echo "Running parse..."

metaqc parse \
    "$INPUT_FILE" \
    --config "$PARSE_CONFIG" \
    --output "$OUTPUT_FILE"

echo "Done:"
echo "$OUTPUT_FILE"
