
# 030-qc-fastq

Esta etapa faz controle de qualidade e trimming dos FASTQs usando o metadata
final da etapa `025-parse` como fonte.

A unidade tecnica e o `run_accession`. A unidade biologica e o `sample_id`.
Por isso, o pipeline processa primeiro cada run tecnico e so depois junta os
runs que pertencem a uma mesma amostra biologica.

## Visao Geral

```text
FASTQs brutos ou renomeados
  -> FastQC bruto por run
  -> Trim Galore por run
  -> FastQC dos runs trimmados
  -> merge por sample_id
  -> FastQC das amostras merged
  -> MultiQC por projeto
```

## Entradas

O metadata final deve estar em:

```text
../025-parse/030-metadata_final/AllProjects_metadata.csv
```

Os FASTQs devem estar no scratch do projeto:

```text
/scratch/<PROJECT>/fastq_ftp
```

O plano aceita FASTQs ainda no nome original do ENA:

```text
ERRxxxx_1.fastq.gz
ERRxxxx_2.fastq.gz
```

ou FASTQs ja renomeados pela etapa `020`:

```text
<file_prefix>_<run_accession>_R1.fastq.gz
<file_prefix>_<run_accession>_R2.fastq.gz
```

Nao e necessario colocar diretorio completo nos exemplos de reads. O caminho
real e resolvido pelo `qc_plan.csv`.

## Configuracao

Os caminhos e ambientes padrao estao centralizados em:

```text
../config/pipeline_config.sh
```

Variaveis importantes:

```bash
SCRATCH_ROOT=/scratch
RNA_TOOLS_ENV=rna-tools
PYTHON_ENV=python-list
METADATA_FINAL=../025-parse/030-metadata_final/AllProjects_metadata.csv
```

Scripts de FastQC, Trim Galore e MultiQC usam `activate_rna_tools`.
Scripts auxiliares em Python usam `activate_python_env`.

## Execucao Recomendada

Use o orquestrador:

```bash
bash run_qc_project.sh PRJEB14695
```

Para outro projeto:

```bash
bash run_qc_project.sh PRJEB32839
```

Antes de submeter jobs, confira os comandos:

```bash
bash run_qc_project.sh PRJEB14695 --dry-run
```

Com opcoes explicitas:

```bash
bash run_qc_project.sh PRJEB32839 \
  --metadata ../025-parse/030-metadata_final/AllProjects_metadata.csv \
  --scratch-root /scratch \
  --run-concurrency 10 \
  --sample-concurrency 10
```

O orquestrador:

```text
1. gera work/<PROJECT>_qc_plan.csv
2. calcula o numero de runs e amostras
3. submete FastQC bruto por run
4. submete trimming por run
5. submete FastQC dos runs trimmados apos o trimming
6. submete merge por sample_id apos o trimming
7. submete FastQC das amostras merged apos o merge
8. submete MultiQC apos os FastQCs
```

O FastQC bruto pode rodar em paralelo com o trimming. As etapas que dependem do
trimming ficam em estado `Dependency` no Slurm ate a etapa anterior terminar
com sucesso.

## Execucao Manual

Use esta forma apenas quando quiser rodar etapa por etapa.

### 1. Gerar o plano

```bash
python generate_qc_plan.py \
  --metadata ../025-parse/030-metadata_final/AllProjects_metadata.csv \
  --project PRJEB14695 \
  --scratch-root /scratch \
  --output work/PRJEB14695_qc_plan.csv
```

O plano tem uma linha por run tecnico e contem os caminhos usados pelos scripts:

```text
dataset
sample_id
file_prefix
run_accession
raw_r1
raw_r2
trimmed_run_r1
trimmed_run_r2
merged_sample_r1
merged_sample_r2
```

### 2. Conferir tamanho dos arrays

```bash
python plan_counts.py work/PRJEB14695_qc_plan.csv
```

Exemplo:

```text
runs=138
samples=23
```

Use `runs` nos jobs por run e `samples` nos jobs por amostra biologica.

### 3. FastQC bruto por run

```bash
sbatch --array=1-N%10 fastqc_raw_plan.sh \
  work/PRJEB14695_qc_plan.csv \
  /scratch/PRJEB14695/fastqc_raw
```

### 4. Trimming por run

```bash
sbatch --array=1-N%10 trim_runs_plan.sh \
  work/PRJEB14695_qc_plan.csv
```

### 5. FastQC dos runs trimmados

```bash
sbatch --array=1-N%10 fastqc_trimmed_runs_plan.sh \
  work/PRJEB14695_qc_plan.csv \
  /scratch/PRJEB14695/fastqc_trimmed_runs
```

### 6. Merge por amostra biologica

```bash
sbatch --array=1-M%10 merge_samples_plan.sh \
  work/PRJEB14695_qc_plan.csv
```

### 7. FastQC das amostras merged

```bash
sbatch --array=1-M%10 fastqc_merged_plan.sh \
  work/PRJEB14695_qc_plan.csv \
  /scratch/PRJEB14695/fastqc_merged
```

### 8. MultiQC

```bash
sbatch multiqc_plan.sh PRJEB14695 /scratch
```

## Saidas

Para cada projeto, as saidas ficam organizadas em:

```text
/scratch/<PROJECT>/fastqc_raw
/scratch/<PROJECT>/trimmed_runs
/scratch/<PROJECT>/fastqc_trimmed_runs
/scratch/<PROJECT>/trimmed_merged
/scratch/<PROJECT>/fastqc_merged
/scratch/<PROJECT>/multiqc_030
```

Os FASTQs finais para quantificacao ficam em:

```text
/scratch/<PROJECT>/trimmed_merged
```

com nomes:

```text
<sample_id>_R1_trimmed.fastq.gz
<sample_id>_R2_trimmed.fastq.gz
```
