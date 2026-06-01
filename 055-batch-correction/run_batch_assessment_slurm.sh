#!/bin/bash
#
# Submit step 055 batch assessment to SLURM.
#
# Usage:
#   bash run_batch_assessment_slurm.sh --all --batch-column dataset --covariates life_stage,tissue,sex
#   bash run_batch_assessment_slurm.sh --all --sbatch-dry-run
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
source "${PROJECT_DIR}/config/pipeline_config.sh"
cd "$SCRIPT_DIR"

SBATCH_DRY_RUN=0
ARGS=()

while [ $# -gt 0 ]; do
    case "$1" in
        --sbatch-dry-run)
            SBATCH_DRY_RUN=1
            shift
            ;;
        *)
            ARGS+=("$1")
            shift
            ;;
    esac
done

if [ "${#ARGS[@]}" -eq 0 ]; then
    echo "Uso: $0 [PROJECT|--all] [opcoes de run_batch_assessment.sh] [--sbatch-dry-run]"
    exit 1
fi

mkdir -p logs/batch

CMD=(
    sbatch --parsable
    --export="ALL,PROJECT_DIR=${PROJECT_DIR},PIPELINE_CONFIG=${PROJECT_DIR}/config/pipeline_config.sh,STEP_DIR=${SCRIPT_DIR}"
    "${SCRIPT_DIR}/batch_assessment_job.sh"
    "${ARGS[@]}"
)

echo "[INFO] Submetendo avaliacao de batch via SLURM"
echo "[INFO] STEP_DIR: $SCRIPT_DIR"
echo "[INFO] PROJECT_DIR: $PROJECT_DIR"

if [ "$SBATCH_DRY_RUN" -eq 1 ]; then
    echo "[DRY-RUN] Comando:"
    printf ' %q' "${CMD[@]}"
    echo
    exit 0
fi

JOB_ID=$("${CMD[@]}" | tail -n 1 | cut -d';' -f1)
echo "[OK] Job batch assessment submetido: $JOB_ID"
