#!/bin/bash
#
# Coordinate step 040 for one project using the step-030 qc_plan.
#
# Usage:
#   bash run_alignment_project.sh PRJEB32839 ../030-qc-fastq/work/PRJEB32839_qc_plan.csv
#   bash run_alignment_project.sh PRJEB14695 ../030-qc-fastq/work/PRJEB14695_qc_plan.csv --dry-run
#

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Uso: $0 <PROJECT> <QC_PLAN.csv> [opcoes]"
    echo ""
    echo "Opcoes:"
    echo "  --index PATH          Default: SALMON_INDEX_DIR do config/pipeline_config.sh"
    echo "  --output-root PATH    Default: QUANT_DIR do config/pipeline_config.sh"
    echo "  --plan PATH           Default: work/<PROJECT>_salmon_plan.csv"
    echo "  --concurrency N       Default: 10"
    echo "  --allow-missing       Permite gerar plano mesmo com FASTQs merged ausentes"
    echo "  --dry-run             Mostra comandos sem submeter jobs"
    exit 1
fi

PROJECT=$1
QC_PLAN=$2
shift 2

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
source "${PROJECT_DIR}/config/pipeline_config.sh"
cd "$SCRIPT_DIR"

INDEX_DIR="$SALMON_INDEX_DIR"
OUTPUT_ROOT="$QUANT_DIR"
PLAN=""
CONCURRENCY=10
ALLOW_MISSING=0
DRY_RUN=0

while [ $# -gt 0 ]; do
    case "$1" in
        --index)
            INDEX_DIR=$2
            shift 2
            ;;
        --output-root)
            OUTPUT_ROOT=$2
            shift 2
            ;;
        --plan)
            PLAN=$2
            shift 2
            ;;
        --concurrency)
            CONCURRENCY=$2
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
        *)
            echo "[ERRO] Opcao desconhecida: $1"
            exit 1
            ;;
    esac
done

if [ -z "$PLAN" ]; then
    PLAN="work/${PROJECT}_salmon_plan.csv"
fi

mkdir -p work logs

run_cmd() {
    echo "+ $*"
    if [ "$DRY_RUN" -eq 0 ]; then
        "$@"
    fi
}

submit_job() {
    local output
    echo "+ $*" >&2
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "DRYRUN_JOB"
        return 0
    fi
    output=$("$@")
    echo "$output" >&2
    echo "$output" | tail -n 1 | cut -d';' -f1
}

echo "[INFO] Projeto: $PROJECT"
echo "[INFO] QC plan: $QC_PLAN"
echo "[INFO] Salmon index: $INDEX_DIR"
echo "[INFO] Output root: $OUTPUT_ROOT"
echo "[INFO] Salmon plan: $PLAN"

GEN_PLAN_CMD=(
    python generate_salmon_plan.py
    --qc-plan "$QC_PLAN"
    --project "$PROJECT"
    --output-root "$OUTPUT_ROOT"
    --output "$PLAN"
)

if [ "$ALLOW_MISSING" -eq 1 ]; then
    GEN_PLAN_CMD+=(--allow-missing)
fi

activate_python_env
run_cmd "${GEN_PLAN_CMD[@]}"

if [ "$DRY_RUN" -eq 1 ]; then
    echo "[DRY-RUN] Exemplo de submissao:"
    echo "  sbatch --parsable --array=1-SAMPLES%${CONCURRENCY} salmon_quant_plan.sh $PLAN $INDEX_DIR"
    exit 0
fi

COUNTS=$(python salmon_plan_counts.py "$PLAN")
SAMPLES=$(echo "$COUNTS" | awk -F= '$1=="samples" {print $2}')

if [ -z "$SAMPLES" ]; then
    echo "[ERRO] Nao foi possivel ler samples de $PLAN"
    echo "$COUNTS"
    exit 1
fi

echo "[INFO] Amostras para Salmon: $SAMPLES"

SALMON_JOB=$(
    submit_job sbatch --parsable \
        --export="ALL,PROJECT_DIR=${PROJECT_DIR},PIPELINE_CONFIG=${PROJECT_DIR}/config/pipeline_config.sh" \
        --array="1-${SAMPLES}%${CONCURRENCY}" \
        "${SCRIPT_DIR}/salmon_quant_plan.sh" "$PLAN" "$INDEX_DIR"
)

echo "[OK] Job Salmon submetido: $SALMON_JOB"
