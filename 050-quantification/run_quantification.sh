#!/bin/bash
#
# Coordinate step 050: import Salmon quant.sf files into gene-level matrices.
#
# Usage:
#   bash run_quantification.sh PRJEB32839
#   bash run_quantification.sh --all
#   bash run_quantification.sh PRJEB32839 --allow-missing --dry-run
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
source "${PROJECT_DIR}/config/pipeline_config.sh"
cd "$SCRIPT_DIR"

DEFAULT_METADATA="$METADATA_FINAL"
if [ ! -f "$DEFAULT_METADATA" ]; then
    DEFAULT_METADATA="$METADATA_FINAL"
fi

PROJECT=""
METADATA="$DEFAULT_METADATA"
QUANT_ROOT="$QUANT_DIR"
GTF="$REF_GTF"
OUTPUT_DIR="$QUANTIFICATION_DIR"
COUNTS_NAME=""
TPM_NAME=""
SAMPLE_TABLE_NAME=""
TX2GENE_OUT=""
ALLOW_MISSING=0
DRY_RUN=0

usage() {
    echo "Uso: $0 [PROJECT|--all] [opcoes]"
    echo ""
    echo "Opcoes:"
    echo "  --metadata PATH          Default: METADATA_FINAL, com fallback para METADATA_FINAL"
    echo "  --quant-root PATH        Default: QUANT_DIR do config/pipeline_config.sh"
    echo "  --gtf PATH               Default: REF_GTF do config/pipeline_config.sh"
    echo "  --output-dir PATH        Default: QUANTIFICATION_DIR do config/pipeline_config.sh"
    echo "  --counts-name NAME       Nome do arquivo de counts"
    echo "  --tpm-name NAME          Nome do arquivo de TPM"
    echo "  --sample-table-name NAME Nome da tabela de amostras importadas"
    echo "  --tx2gene-out PATH       Default: <output-dir>/tx2gene.tsv"
    echo "  --allow-missing          Importa somente quant.sf existentes"
    echo "  --dry-run                Mostra o comando sem executar"
    exit 1
}

if [ $# -eq 0 ]; then
    usage
fi

while [ $# -gt 0 ]; do
    case "$1" in
        --all)
            PROJECT=""
            shift
            ;;
        --metadata)
            METADATA=$2
            shift 2
            ;;
        --quant-root)
            QUANT_ROOT=$2
            shift 2
            ;;
        --gtf)
            GTF=$2
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR=$2
            shift 2
            ;;
        --counts-name)
            COUNTS_NAME=$2
            shift 2
            ;;
        --tpm-name)
            TPM_NAME=$2
            shift 2
            ;;
        --sample-table-name)
            SAMPLE_TABLE_NAME=$2
            shift 2
            ;;
        --tx2gene-out)
            TX2GENE_OUT=$2
            shift 2
            ;;
        --allow-missing)
            ALLOW_MISSING=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            usage
            ;;
        --*)
            echo "[ERRO] Opcao desconhecida: $1"
            exit 1
            ;;
        *)
            if [ -n "$PROJECT" ]; then
                echo "[ERRO] Projeto informado mais de uma vez: $PROJECT e $1"
                exit 1
            fi
            PROJECT=$1
            shift
            ;;
    esac
done

if [ ! -f "$METADATA" ]; then
    echo "[ERRO] Metadata nao encontrado: $METADATA"
    exit 1
fi

if [ ! -d "$QUANT_ROOT" ]; then
    echo "[ERRO] Quant root nao encontrado: $QUANT_ROOT"
    exit 1
fi

if [ ! -f "$GTF" ]; then
    echo "[ERRO] GTF nao encontrado: $GTF"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

CMD=(
    Rscript txtimport_quant.R
    --metadata "$METADATA"
    --quant-root "$QUANT_ROOT"
    --gtf "$GTF"
    --output-dir "$OUTPUT_DIR"
)

if [ -n "$PROJECT" ]; then
    CMD+=(--project "$PROJECT")
fi

if [ -n "$COUNTS_NAME" ]; then
    CMD+=(--counts-name "$COUNTS_NAME")
fi

if [ -n "$TPM_NAME" ]; then
    CMD+=(--tpm-name "$TPM_NAME")
fi

if [ -n "$SAMPLE_TABLE_NAME" ]; then
    CMD+=(--sample-table-name "$SAMPLE_TABLE_NAME")
fi

if [ -n "$TX2GENE_OUT" ]; then
    CMD+=(--tx2gene-out "$TX2GENE_OUT")
fi

if [ "$ALLOW_MISSING" -eq 1 ]; then
    CMD+=(--allow-missing)
fi

echo "[INFO] Projeto: ${PROJECT:-TODOS}"
echo "[INFO] Metadata: $METADATA"
echo "[INFO] Quant root: $QUANT_ROOT"
echo "[INFO] GTF: $GTF"
echo "[INFO] Output dir: $OUTPUT_DIR"

if [ "$DRY_RUN" -eq 1 ]; then
    echo "[DRY-RUN] Comando:"
    printf ' %q' "${CMD[@]}"
    echo
    exit 0
fi

activate_r_analysis
check_command Rscript

echo "+ ${CMD[*]}"
"${CMD[@]}"
