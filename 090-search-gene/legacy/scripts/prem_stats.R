#!/usr/bin/env Rscript

source("prem_utils.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Uso: Rscript prem_stats.R <TPM_long_file.csv> <out_dir>")

tpm_file <- args[1]
out_dir  <- args[2]

data <- read_csv(tpm_file)

# Descritivas
stats_desc <- data %>%
  group_by(stage_category) %>%
  summarise(
    n = n(),
    mean_tpm = mean(TPM),
    median_tpm = median(TPM),
    sd_tpm = sd(TPM),
    .groups = "drop"
  )

# Kruskal-Wallis
kruskal_res <- kruskal_test(data, TPM ~ stage_category)

# Pós-hoc Dunn
dunn_res <- dunn_test(data, TPM ~ stage_category, p.adjust.method = "BH")

# Salvar
write_csv(dunn_res, file.path(out_dir, "comparisons.csv"))

sink(file.path(out_dir, "stats_report.txt"))
cat("Relatório estatístico -", basename(out_dir), "\n\n")
print(stats_desc)
print(kruskal_res)
print(dunn_res)
sink()

cat("prem_stats.R concluído para", out_dir, "\n")
