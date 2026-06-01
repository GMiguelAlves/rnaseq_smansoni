# 040-alignment

Esta etapa quantifica as amostras biologicas com Salmon usando os FASTQs merged
gerados na etapa `030-qc-fastq`.

O fluxo novo nao usa `010-reference/RNAseq_metadata.tsv`. A entrada e o
`qc_plan.csv` criado no passo 030.

## Entradas

```text
030-qc-fastq/work/<PROJECT>_qc_plan.csv
010-reference/salmon_index/
```

O `qc_plan.csv` deve conter:

```text
sample_id
merged_sample_r1
merged_sample_r2
```

## Gerar plano Salmon

```bash
source ../config/pipeline_config.sh
cd "$ALIGN_DIR"

python generate_salmon_plan.py \
  --qc-plan ../030-qc-fastq/work/PRJEB32839_qc_plan.csv \
  --project PRJEB32839 \
  --output-root quants \
  --output work/PRJEB32839_salmon_plan.csv
```

Conferir contagem:

```bash
python salmon_plan_counts.py work/PRJEB32839_salmon_plan.csv
```

## Rodar Salmon

Use `samples=N` retornado por `salmon_plan_counts.py`.

```bash
sbatch --array=1-N%10 salmon_quant_plan.sh \
  work/PRJEB32839_salmon_plan.csv \
  "$SALMON_INDEX_DIR"
```

As saidas ficam em:

```text
040-alignment/quants/<PROJECT>/<sample_id>/quant.sf
```

## Execucao coordenada

Para gerar o plano e submeter o array Salmon em um unico comando:

```bash
bash run_alignment_project.sh \
  PRJEB32839 \
  ../030-qc-fastq/work/PRJEB32839_qc_plan.csv
```

Dry run:

```bash
bash run_alignment_project.sh \
  PRJEB32839 \
  ../030-qc-fastq/work/PRJEB32839_qc_plan.csv \
  --dry-run
```

Opcoes uteis:

```bash
bash run_alignment_project.sh \
  PRJEB32839 \
  ../030-qc-fastq/work/PRJEB32839_qc_plan.csv \
  --index "$SALMON_INDEX_DIR" \
  --output-root "$QUANT_DIR" \
  --concurrency 10
```
