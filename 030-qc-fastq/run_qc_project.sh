#!/bin/bash
#
# Coordinate step 030 for one project using the metadata-driven qc_plan.
#
# Usage:
#   bash run_qc_project.sh PRJEB32839
#   bash run_qc_project.sh PRJEB14695 --dry-run
#

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Uso: $0 <PROJECT> [opcoes]"
    echo ""
    echo "Opcoes:"
    echo "  --metadata PATH          Default: ../025-parse/030-metadata_final/AllProjects_metadata_new.csv"
    echo "  --scratch-root PATH      Default: /scratch/Schisto-epigenetics/gustavo"
    echo "  --plan PATH              Default: work/<PROJECT>_qc_plan.csv"
    echo "  --run-concurrency N      Default: 10"
    echo "  --sample-concurrency N   Default: 10"
    echo "  --allow-missing          Permite gerar plano mesmo com FASTQs ausentes"
    echo "  --dry-run                Mostra comandos sem submeter jobs"
    exit 1
fi

PROJECT=$1
shift

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
source "${PROJECT_DIR}/config/pipeline_config.sh"
cd "$SCRIPT_DIR"

DEFAULT_METADATA="$METADATA_FINAL_NEW"
if [ ! -f "$DEFAULT_METADATA" ]; then
    DEFAULT_METADATA="$METADATA_FINAL"
fi

METADATA="$DEFAULT_METADATA"
SCRATCH_ROOT="$SCRATCH_ROOT"
PLAN=""
RUN_CONCURRENCY=10
SAMPLE_CONCURRENCY=10
ALLOW_MISSING=0
DRY_RUN=0

while [ $# -gt 0 ]; do
    case "$1" in
        --metadata)
            METADATA=$2
            shift 2
            ;;
        --scratch-root)
            SCRATCH_ROOT=$2
            shift 2
            ;;
        --plan)
            PLAN=$2
            shift 2
            ;;
        --run-concurrency)
            RUN_CONCURRENCY=$2
            shift 2
            ;;
        --sample-concurrency)
            SAMPLE_CONCURRENCY=$2
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
    PLAN="work/${PROJECT}_qc_plan.csv"
fi

PROJECT_SCRATCH="${SCRATCH_ROOT}/${PROJECT}"
FASTQC_RAW_DIR="${PROJECT_SCRATCH}/fastqc_raw"
FASTQC_TRIMMED_RUNS_DIR="${PROJECT_SCRATCH}/fastqc_trimmed_runs"
FASTQC_MERGED_DIR="${PROJECT_SCRATCH}/fastqc_merged"

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
echo "[INFO] Metadata: $METADATA"
echo "[INFO] Scratch root: $SCRATCH_ROOT"
echo "[INFO] Plano: $PLAN"

GEN_PLAN_CMD=(
    python generate_qc_plan.py
    --metadata "$METADATA"
    --project "$PROJECT"
    --scratch-root "$SCRATCH_ROOT"
    --output "$PLAN"
)

if [ "$ALLOW_MISSING" -eq 1 ]; then
    GEN_PLAN_CMD+=(--allow-missing)
fi

activate_python_env

run_cmd "${GEN_PLAN_CMD[@]}"

if [ "$DRY_RUN" -eq 1 ]; then
    echo "[DRY-RUN] Pulando leitura de contagens do plano."
    echo "[DRY-RUN] Exemplo de submissao seria:"
    echo "  sbatch --parsable --array=1-RUNS%${RUN_CONCURRENCY} fastqc_raw_plan.sh $PLAN $FASTQC_RAW_DIR"
    echo "  sbatch --parsable --array=1-RUNS%${RUN_CONCURRENCY} trim_runs_plan.sh $PLAN"
    echo "  sbatch --parsable --dependency=afterok:<trim_job> --array=1-RUNS%${RUN_CONCURRENCY} fastqc_trimmed_runs_plan.sh $PLAN $FASTQC_TRIMMED_RUNS_DIR"
    echo "  sbatch --parsable --dependency=afterok:<trim_job> --array=1-SAMPLES%${SAMPLE_CONCURRENCY} merge_samples_plan.sh $PLAN"
    echo "  sbatch --parsable --dependency=afterok:<merge_job> --array=1-SAMPLES%${SAMPLE_CONCURRENCY} fastqc_merged_plan.sh $PLAN $FASTQC_MERGED_DIR"
    echo "  sbatch --parsable --dependency=afterok:<raw_job>:<trim_qc_job>:<merged_qc_job> multiqc_plan.sh $PROJECT $SCRATCH_ROOT"
    exit 0
fi

COUNTS=$(python plan_counts.py "$PLAN")
RUNS=$(echo "$COUNTS" | awk -F= '$1=="runs" {print $2}')
SAMPLES=$(echo "$COUNTS" | awk -F= '$1=="samples" {print $2}')

if [ -z "$RUNS" ] || [ -z "$SAMPLES" ]; then
    echo "[ERRO] Nao foi possivel ler runs/samples de $PLAN"
    echo "$COUNTS"
    exit 1
fi

echo "[INFO] Runs tecnicos: $RUNS"
echo "[INFO] Amostras biologicas: $SAMPLES"

RAW_QC_JOB=$(
    submit_job sbatch --parsable \
        --export=ALL,PROJECT_DIR="$PROJECT_DIR",SCRIPT_DIR="$SCRIPT_DIR" \
        --array="1-${RUNS}%${RUN_CONCURRENCY}" \
        "${SCRIPT_DIR}/fastqc_raw_plan.sh" "$PLAN" "$FASTQC_RAW_DIR"
)

TRIM_JOB=$(
    submit_job sbatch --parsable \
        --export=ALL,PROJECT_DIR="$PROJECT_DIR",SCRIPT_DIR="$SCRIPT_DIR" \
        --array="1-${RUNS}%${RUN_CONCURRENCY}" \
        "${SCRIPT_DIR}/trim_runs_plan.sh" "$PLAN"
)

TRIM_QC_JOB=$(
    submit_job sbatch --parsable \
        --export=ALL,PROJECT_DIR="$PROJECT_DIR",SCRIPT_DIR="$SCRIPT_DIR" \
        --dependency="afterok:${TRIM_JOB}" \
        --array="1-${RUNS}%${RUN_CONCURRENCY}" \
        "${SCRIPT_DIR}/fastqc_trimmed_runs_plan.sh" "$PLAN" "$FASTQC_TRIMMED_RUNS_DIR"
)

MERGE_JOB=$(
    submit_job sbatch --parsable \
        --export=ALL,PROJECT_DIR="$PROJECT_DIR",SCRIPT_DIR="$SCRIPT_DIR" \
        --dependency="afterok:${TRIM_JOB}" \
        --array="1-${SAMPLES}%${SAMPLE_CONCURRENCY}" \
        "${SCRIPT_DIR}/merge_samples_plan.sh" "$PLAN"
)

MERGED_QC_JOB=$(
    submit_job sbatch --parsable \
        --export=ALL,PROJECT_DIR="$PROJECT_DIR",SCRIPT_DIR="$SCRIPT_DIR" \
        --dependency="afterok:${MERGE_JOB}" \
        --array="1-${SAMPLES}%${SAMPLE_CONCURRENCY}" \
        "${SCRIPT_DIR}/fastqc_merged_plan.sh" "$PLAN" "$FASTQC_MERGED_DIR"
)

MULTIQC_JOB=$(
    submit_job sbatch --parsable \
        --export=ALL,PROJECT_DIR="$PROJECT_DIR",SCRIPT_DIR="$SCRIPT_DIR" \
        --dependency="afterok:${RAW_QC_JOB}:${TRIM_QC_JOB}:${MERGED_QC_JOB}" \
        "${SCRIPT_DIR}/multiqc_plan.sh" "$PROJECT" "$SCRATCH_ROOT"
)

echo "[OK] Jobs submetidos:"
echo "  FastQC raw:           $RAW_QC_JOB"
echo "  Trim runs:            $TRIM_JOB"
echo "  FastQC trimmed runs:  $TRIM_QC_JOB"
echo "  Merge samples:        $MERGE_JOB"
echo "  FastQC merged:        $MERGED_QC_JOB"
echo "  MultiQC:              $MULTIQC_JOB"
