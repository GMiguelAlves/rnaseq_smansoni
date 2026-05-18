#!/bin/bash
#SBATCH --job-name=merge_meta
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

source "/home/${USER}@bio.ib.unicamp.br/miniconda3/etc/profile.d/conda.sh"
conda activate python-list

INPUT_DIR="../020-metadata_parsers/Allprojects"
OUTPUT_FILE="AllProjects_metadata.csv"

mkdir -p logs/

echo "Searching parsed metadata in: $INPUT_DIR"
echo "Output: $OUTPUT_FILE"

metaqc merge "$INPUT_DIR" \
    --output "$OUTPUT_FILE"


echo "Merge completed successfully."
echo "Generated:"
echo "$OUTPUT_FILE"
