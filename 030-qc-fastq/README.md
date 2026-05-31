

---

## Fluxo orientado por `qc_plan.csv`

O fluxo novo usa o metadata final da etapa `025-parse` como fonte unica de
verdade. A unidade tecnica e o `run_accession`; a unidade biologica e o
`sample_id`.

Ordem recomendada:

```text
FASTQs renomeados
  -> FastQC bruto por run
  -> Trim Galore por run
  -> FastQC dos runs trimmados
  -> merge/concatenacao por sample_id
  -> FastQC das amostras biologicas merged
  -> MultiQC por projeto
```

### 1. Gerar o plano

```bash
cd 

python generate_qc_plan.py \
  --metadata ../025-parse/030-metadata_final/AllProjects_metadata.csv \
  --project PRJEB14695 \
  --scratch-root /scratch \
  --output work/PRJEB14695_qc_plan.csv
```

O plano contem uma linha por run tecnico:

```text
dataset,sample_id,run_accession,raw_r1,raw_r2,trimmed_run_r1,trimmed_run_r2,merged_sample_r1,merged_sample_r2
```

Conferir tamanho dos arrays:

```bash
python plan_counts.py work/PRJEB14695_qc_plan.csv
```

### 2. FastQC bruto por run

Use `runs=N` retornado por `plan_counts.py`.

```bash
sbatch --array=1-N%10 fastqc_raw_plan.sh \
  work/PRJEB14695_qc_plan.csv \
  /scratch/PRJEB14695/fastqc_raw
```

### 3. Trimming por run

```bash
sbatch --array=1-N%10 trim_runs_plan.sh \
  work/PRJEB14695_qc_plan.csv
```

### 4. FastQC dos runs trimmados

```bash
sbatch --array=1-N%10 fastqc_trimmed_runs_plan.sh \
  work/PRJEB14695_qc_plan.csv \
  /scratch/PRJEB14695/fastqc_trimmed_runs
```

### 5. Merge por amostra biologica

Use `samples=M` retornado por `plan_counts.py`.

```bash
sbatch --array=1-M%10 merge_samples_plan.sh \
  work/PRJEB14695_qc_plan.csv
```

### 6. FastQC das amostras merged

```bash
sbatch --array=1-M%10 fastqc_merged_plan.sh \
  work/PRJEB14695_qc_plan.csv \
  /scratch/PRJEB14695/fastqc_merged
```

### 7. MultiQC

```bash
sbatch multiqc_plan.sh \
  PRJEB14695 \
  /scratch
```

Repita os mesmos comandos trocando `PRJEB14695` por `PRJEB32839`.

Observacao: se os FASTQs ainda nao foram renomeados, `generate_qc_plan.py`
tenta usar os nomes ENA originais (`ERRxxxx_1.fastq.gz`). Se ja foram
renomeados pela etapa `020`, ele usa os nomes novos com `file_prefix`.

### Execucao coordenada em um unico comando

Para submeter todas as etapas acima com dependencias Slurm:

```bash
bash run_qc_project.sh PRJEB14695
```

O script:

```text
1. gera work/PRJEB14695_qc_plan.csv
2. calcula runs e samples
3. submete FastQC bruto por run
4. submete trimming por run
5. submete FastQC dos runs trimmados apos o trimming
6. submete merge por sample_id apos o trimming
7. submete FastQC das amostras merged apos o merge
8. submete MultiQC apos os FastQCs
```

Os ambientes e caminhos padrao vêm de `../config/pipeline_config.sh`:

```bash
SCRATCH_ROOT=/scratch
RNA_TOOLS_ENV=rna-tools
PYTHON_ENV=python-list
METADATA_FINAL_NEW=../025-parse/030-metadata_final/AllProjects_metadata.csv
```

Os scripts `*_plan.sh` carregam essa config automaticamente. Scripts de
FastQC/Trim/MultiQC usam `activate_rna_tools`; o merge por amostra usa
`activate_python_env`.

Antes de submeter, confira os comandos:

```bash
bash run_qc_project.sh PRJEB14695 --dry-run
```

Opcoes uteis:

```bash
bash run_qc_project.sh PRJEB32839 \
  --metadata ../025-parse/030-metadata_final/AllProjects_metadata.csv \
  --scratch-root /scratch \
  --run-concurrency 10 \
  --sample-concurrency 10
```

---
---
