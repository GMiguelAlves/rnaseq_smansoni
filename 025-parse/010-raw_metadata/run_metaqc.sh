#!/bin/bash
#SBATCH --job-name=metaqc
#SBATCH --output=logs/metaqc_%A.out
#SBATCH --error=logs/err/metaqc_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00

set -euo pipefail

source "/home/${USER}@bio.ib.unicamp.br/miniconda3/bin/activate"
conda activate python-list

wget -O PRJEB32839.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB32839&result=read_run&fields=study_accession,sample_accession,run_accession,tax_id,scientific_name,library_name,center_name,study_title,fastq_ftp,submitted_ftp,sample_alias&format=tsv&download=true&limit=0"

metaqc validate PRJEB32839.tsv \
--keep sample_id,run_accession,study_accession  \
--output ../015-intermediate_folder/PRJEB32839_base.csv \
--schema-file /home/${USER}@bio.ib.unicamp.br/metaQC/schemas/rnaseq.yaml

wget -O PRJEB14695.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB14695&result=read_run&fields=study_accession,sample_accession,run_accession,tax_id,scientific_name,library_name,center_name,study_title,fastq_ftp,submitted_ftp,sample_alias&format=tsv&download=true&limit=0"

metaqc validate PRJEB14695.tsv  \
--keep sample_id,run_accession,study_accession  \
--output ../015-intermediate_folder/PRJEB14695_base.csv \
--schema-file /home/${USER}@bio.ib.unicamp.br/metaQC/schemas/rnaseq.yaml

rm validation_report.txt
