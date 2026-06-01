# Etapa 090 - Relatorio de genes candidatos

Esta etapa consolida expressao, anotacao e resultados de DEG para genes ou grupos de genes de interesse.

## Entrada principal

Edite `genes.txt` no formato:

```text
Modificadores de histonas: gene1, genex, genez
Grupo 2: gene3, geney
```

Os nomes podem ser IDs de gene ou nomes anotados no GFF3. Comentarios com `#` e linhas vazias sao ignorados.

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

- `gene_set_report.html`: relatorio final navegavel.
- `tables/gene_catalog.tsv`: grupos, queries, genes encontrados e anotacao.
- `tables/expression_long.tsv`: expressao TPM por gene e amostra.
- `tables/expression_summary.tsv`: resumo por projeto, estagio, tecido, sexo e condicao.
- `tables/deg_hits.tsv`: resultados de DEG da etapa 060 filtrados para os genes de interesse.
- `tables/gene_candidate_scores.tsv`: ranking integrado por gene, combinando expressao, DEG e sinal ovario/testiculo.
- `plots/candidate_genes_expression_heatmap.png`: heatmap integrativo de expressao.
- `plots/candidate_genes_expression_dotplot.png`: dotplot de expressao media e fracao expressa.
- `plots/candidate_genes_deg_log2fc_heatmap.png`: heatmap de log2FC em contrastes DEG.
- `plots/candidate_genes_deg_significant_barplot.png`: numero de contrastes significativos por gene.
- `plots/candidate_genes_integrated_score.png`: ranking visual dos principais candidatos.
- `plots/candidate_genes_expression_vs_deg.png`: comparacao entre expressao TPM clara e recorrencia DEG.
- `plots/candidate_genes_deg_consistency_tile.png`: consistencia gene x contraste/projeto/modo de analise.
- `plots/candidate_genes_tissue_sex_heatmap.png`: padroes de expressao por tecido e sexo.
- `plots/candidate_genes_ovary_vs_testis.png`: comparacao direta entre ovario e testiculo.
- `genes/<gene>/`: boxplots e perfis de expressao individuais.

## Como interpretar o relatorio

O HTML prioriza quatro perguntas biologicas:

- Quais genes do `genes.txt` sao melhores candidatos quando expressao, DEG e especificidade gonadal sao combinados?
- Quais genes tem expressao clara em TPM, mesmo quando nao sao os mais significativos estatisticamente?
- Quais genes aparecem de forma consistente em mais de um projeto, contraste ou modo de analise?
- Quais genes mostram padroes especificos de tecido/sexo, especialmente ovario e testiculo?

O `candidate_score` nao substitui a interpretacao biologica. Ele e uma triagem para ordenar genes candidatos e destacar casos que combinam sinal estatistico com expressao biologicamente interpretavel.

## Scripts antigos

Os scripts antigos foram preservados em `legacy/` para referencia.
