#!/usr/bin/env Rscript

source("prem_utils.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Uso: Rscript prem_graphs.R <TPM_file> <gene>")

tpm_file <- args[1]
gene <- args[2]
out_dir <- file.path("../results", gene)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Carregar dados
gene_data <- load_tpm_data(tpm_file, gene)
write_csv(gene_data, file.path(out_dir, "TPM_data_long.csv"))

# Estatísticas
stats <- compute_stats(gene_data)
write_csv(stats, file.path(out_dir, "stats_by_category.csv"))

# Gráficos
plot_boxplot(gene_data, gene, file.path(out_dir, "boxplot.png"))
plot_barplot(stats, gene, file.path(out_dir, "barplot.png"))
plot_lineplot(gene_data, gene, file.path(out_dir, "lineplot.png"))
plot_violin(gene_data, gene, file.path(out_dir, "violin.png"))
plot_box_sex(gene_data, gene, file.path(out_dir, "boxplot_sex.png"))
plot_heatmap(gene_data, gene, file.path(out_dir, "heatmap.png"))

cat("prem_graphs.R concluído para", gene, "\n")
