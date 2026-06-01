# Etapa 055 - Correcao de batch effect

Esta etapa aplica correcao de batch effect sobre a matriz gene-level de counts produzida na etapa 050.

A ferramenta usada e `pycombat_seq`, implementacao Python de ComBat-Seq no pacote `inmoose`. Ela foi escolhida porque o paper  descreve o `pyComBat`/`inmoose` como uma implementacao Python de ComBat e ComBat-Seq com o mesmo arcabouco matematico das versoes originais. Para RNA-seq, a funcao adequada e `pycombat_seq`, pois modela counts com distribuicao binomial negativa.

## Quando usar

Use principalmente em analises combinando projetos/lotes diferentes, por exemplo a matriz `counts_matrix.tsv` gerada por:

```bash
bash ../050-quantification/run_quantification_slurm.sh --all
```

Nao aplique cegamente quando o batch estiver confundido com a biologia. Por exemplo: se cada projeto tiver apenas um estagio, tecido ou sexo exclusivo, a correcao pode remover sinal biologico real. O script tenta detectar esse problema quando covariaveis biologicas sao informadas.

## Ambiente

Crie o ambiente uma vez:

```bash
conda env create -f ../envs/batch-correction.yml
```

O nome do ambiente pode ser sobrescrito por `BATCH_CORRECTION_ENV` em `config/pipeline_config.sh`.

## Execucao

Execucao preferencial via SLURM para avaliar se ha evidencia de batch effect antes de corrigir:

```bash
bash run_batch_assessment_slurm.sh --all --batch-column dataset --covariates life_stage,tissue,sex
```

Execucao direta, util para debug fora do SLURM:

```bash
bash run_batch_assessment.sh --all --batch-column dataset --covariates life_stage,tissue,sex
```

Correcao entre projetos via SLURM, usando `dataset` como batch:

```bash
bash run_batch_correction_slurm.sh --all --batch-column dataset
```

Correcao direta, util para debug:

```bash
bash run_batch_correction.sh --all --batch-column dataset
```

Preservando covariaveis biologicas categoricas via SLURM:

```bash
bash run_batch_correction_slurm.sh --all --batch-column dataset --covariates life_stage,tissue,sex
```

Equivalente direto:

```bash
bash run_batch_correction.sh --all --batch-column dataset --covariates life_stage,tissue,sex
```

Avaliar a eficiencia da correcao depois de gerar `counts_batch_corrected.tsv`, via SLURM:

```bash
bash run_batch_assessment_slurm.sh --all --batch-column dataset --covariates life_stage,tissue,sex
```

Equivalente direto:

```bash
bash run_batch_assessment.sh --all --batch-column dataset --covariates life_stage,tissue,sex
```

Quando a matriz corrigida existe no diretorio default da 055, o wrapper a detecta automaticamente. Tambem e possivel informar explicitamente:

```bash
bash run_batch_assessment_slurm.sh --all --batch-column dataset --corrected-counts all_projects/counts_batch_corrected.tsv
```

Simular sem executar:

```bash
bash run_batch_correction_slurm.sh --all --sbatch-dry-run
```

Para um projeto individual, normalmente e necessario indicar uma coluna tecnica real, como `sequencing_batch`, se existir:

```bash
bash run_batch_correction_slurm.sh PRJEB32839 --batch-column sequencing_batch --skip-if-single-batch
```

## Saidas

Por default, as saidas ficam em:

- `055-batch-correction/all_projects/` para `--all`
- `055-batch-correction/<PROJECT>/` para um projeto especifico

Arquivos gerados:

- `counts_batch_corrected.tsv`: counts corrigidos por batch.
- `batch_correction_samples.tsv`: tabela de amostras alinhada a matriz corrigida.
- `batch_correction_report.json`: parametros, avisos e resumo da execucao.
- `batch_pca_before_after.png`: PCA log2(counts + 1) antes/depois, colorido por batch.
- `assessment/batch_effect_metrics.tsv`: metricas de associacao entre batch/covariaveis e a expressao.
- `assessment/batch_effect_assessment.json`: interpretacao resumida da presenca de batch e da eficiencia da correcao.
- `assessment/batch_effect_pc_scores.tsv`: coordenadas de PCA usadas na avaliacao.
- `assessment/batch_effect_pca.png`: PCA antes/depois colorido por batch, quando `matplotlib` estiver disponivel.

## Como interpretar a avaliacao

O avaliador usa a matriz de counts e a tabela de amostras da etapa 050. A matriz e convertida para `logCPM`, os genes de baixa contagem sao filtrados, e os genes mais variaveis da matriz bruta sao usados tanto antes quanto depois da correcao.

As principais metricas sao:

- `weighted_pc_eta2`: quanto da variacao capturada pelos PCs esta associada ao batch.
- `multivariate_r2`: fracao da variacao global de expressao explicada pelo batch.
- `pvalue`: p-valor por permutacao, quando `--permutations` e maior que zero.

Em geral, batch effect forte aparece como `multivariate_r2` alto e/ou p-valor baixo para a coluna de batch. A correcao e considerada eficiente quando o `multivariate_r2` do batch cai substancialmente sem exigir ignorar confundimento batch-biologia.

## Observacao importante

Para analise diferencial com DESeq2, a alternativa estatisticamente mais conservadora costuma ser incluir o batch no design do modelo quando possivel, em vez de testar diretamente counts corrigidos. Esta etapa existe para produzir uma matriz harmonizada e diagnosticos comparativos; a decisao de usar `counts_batch_corrected.tsv` na etapa 060 deve considerar o desenho biologico de cada conjunto de dados.
