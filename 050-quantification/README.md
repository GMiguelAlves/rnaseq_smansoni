# Etapa 050 - Quantificacao

Esta etapa importa os arquivos `quant.sf` gerados pelo Salmon na etapa 040 e produz matrizes gene-level de counts e TPM usando `tximport`.

## Entradas

- Metadata final: `METADATA_FINAL_NEW` definido em `config/pipeline_config.sh`, com fallback para `METADATA_FINAL`.
- Saidas do Salmon: `QUANT_DIR/<PROJECT>/<sample_id>/quant.sf`.
- Anotacao GTF: `REF_GTF`.

Todos os defaults vem de `config/pipeline_config.sh`. Sobrescreva por flags apenas quando necessario.

## Execucao

Execucao preferencial via SLURM, para importar um projeto:

```bash
bash run_quantification_slurm.sh PRJEB32839
```

Execucao preferencial via SLURM, para importar todos os projetos juntos:

```bash
bash run_quantification_slurm.sh --all
```

Simular apenas a submissao SLURM:

```bash
bash run_quantification_slurm.sh PRJEB32839 --sbatch-dry-run
```

Execucao direta, util para debug fora do SLURM:

```bash
bash run_quantification.sh PRJEB32839
```

Execucao direta para todos os projetos juntos:

```bash
bash run_quantification.sh --all
```

Simular sem executar:

```bash
bash run_quantification.sh PRJEB32839 --dry-run
```

Se parte dos arquivos `quant.sf` estiver intencionalmente ausente:

```bash
bash run_quantification_slurm.sh PRJEB32839 --allow-missing
```

## Saidas

Para um projeto:

- `<PROJECT>_counts_matrix.tsv`
- `<PROJECT>_tpm_matrix.tsv`
- `<PROJECT>_quant_samples.tsv`
- `tx2gene.tsv`

Para `--all`:

- `counts_matrix.tsv`
- `tpm_matrix.tsv`
- `quant_samples.tsv`
- `tx2gene.tsv`

Ao importar todos os projetos juntos, os nomes das colunas das matrizes usam `dataset__sample_id` para garantir IDs unicos entre projetos. Em importacoes de um unico projeto, as colunas mantem apenas `sample_id`.
