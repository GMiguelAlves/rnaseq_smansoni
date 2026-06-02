# Etapa 090 - Relatorio exploratorio de genes

Esta etapa consolida expressao, anotacao e resultados de DEG para genes ou grupos de genes de interesse. A etapa nao ranqueia nem classifica genes; ela gera visualizacoes para explorar cada gene e cada grupo em diferentes projetos, batches e contextos biologicos.

## Entrada principal

Edite `genes.txt` no formato:

```text
Modificadores de histonas: Smp_000000, Smp_111111
Grupo 2: Smp_222222, nome_anotado
```

As entradas podem ser IDs `Smp_...` ou nomes anotados no GFF3. Quando houver anotacao disponivel, os graficos e tabelas usam rotulos no formato `nome_do_gene | Smp_xxxx`.

O ciclo de vida e ordenado de forma biologica nos graficos, mas o subestagio original nao e colapsado:

```text
eggs -> miracidium -> sporocyst -> cercariae -> schistosomules -> adult
```

Subestagios com tempo nao somem: `sporocyst_5d`, `sporocyst_10d`, `adult_7d` e `adult_21d` continuam aparecendo separadamente. A etapa cria uma classe ampla apenas para ordenar (`stage_class`) e preserva o valor detalhado em `stage`.

## Execucao

Preferencialmente via SLURM:

```bash
bash run_gene_report_slurm.sh --genes genes.txt --title "Genes epigeneticos"
```

Simular submissao:

```bash
bash run_gene_report_slurm.sh --genes genes.txt --sbatch-dry-run
```

Por default, a etapa usa:

- `050-quantification/tpm_matrix.tsv`
- `050-quantification/quant_samples.tsv`
- `060-deg-analysis/`
- `REF_GFF3`
- `METADATA_FINAL_NEW`

Todos os defaults vem de `config/pipeline_config.sh`.

## Saidas

As saidas ficam em `090-search-gene/results/`:

- `gene_set_report.html`: relatorio HTML navegavel.
- `tables/gene_catalog.tsv`: genes do `genes.txt`, anotacao, `gene_display_label` e status na matriz TPM.
- `tables/gene_expression_summary.tsv`: resumo descritivo por gene, sem ranking.
- `tables/expression_long.tsv`: expressao TPM por gene e amostra.
- `tables/expression_summary_by_context.tsv`: resumo por projeto, batch, condicao, estagio detalhado, tecido e sexo.
- `tables/deg_hits.tsv`: resultados de DEG da etapa 060 filtrados para os genes de interesse.
- `plots/all_groups_*.png`: visualizacoes integradas de todos os grupos.
- `groups/<grupo>/`: heatmaps, dotplots e graficos de batch/DEG por grupo.
- `genes/<grupo>/<gene>/`: graficos aprofundados por gene.

As visualizacoes globais incluem:

- `all_groups_expression_heatmap.png`: expressao media por contexto.
- `all_groups_expression_dotplot.png`: expressao media e fracao expressa.
- `all_groups_sample_heatmap_annotated.png`: heatmap gene x amostra com anotacoes de projeto, batch e biologia.
- `all_groups_gene_correlation.png`: correlacao de expressao entre genes.
- `all_groups_sample_pca.png`: PCA das amostras usando apenas os genes do `genes.txt`.
- `all_groups_sample_mds.png`: MDS das amostras usando apenas os genes do `genes.txt`.
- `all_groups_tissue_sex_heatmap.png`: padroes por tecido e sexo.
- `all_groups_ovary_testis_panel.png`: comparacao ovario versus testiculo.
- `all_groups_aggregate_profile.png`: perfil agregado medio por grupo.
- `all_groups_batch_project_boxplot.png`: distribuicao de expressao por batch/projeto.
- `all_groups_deg_log2fc_heatmap.png`: log2FC nos contrastes DEG.
- `all_groups_deg_context_tile.png`: sinais DEG por contraste/projeto.
- `all_groups_deg_direction_summary.png`: direcao DEG, separando up, down e nao significativo.

Cada pasta em `groups/<grupo>/` replica as principais visualizacoes para o grupo especifico.

## Visualizacoes por gene

Para cada gene sao gerados, quando os dados permitem:

- expressao por tecido, sexo, condicao, estagio detalhado e projeto;
- expressao por batch e projeto;
- perfil medio por estagio detalhado, tecido, batch e condicao;
- expressao em cada amostra individual;
- grafico de efeitos DEG por contraste;
- dispersao `log2FC` versus `-log10(padj)`.

## Scripts antigos

Os scripts antigos foram preservados em `legacy/` para referencia.

