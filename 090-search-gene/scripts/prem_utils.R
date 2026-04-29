
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(rstatix)
  library(seqinr)
  library(rtracklayer)
})

# ============================================================
# Funções de utilidade
# ============================================================

load_tpm_data <- function(tpm_file, gene_target) {
  tpm_matrix <- readr::read_tsv(tpm_file, show_col_types = FALSE)

  gene_data <- tpm_matrix %>%
    dplyr::filter(stringr::str_detect(gene, stringr::fixed(gene_target))) %>%
    tidyr::pivot_longer(cols = -gene, names_to = "sample_full", values_to = "TPM") %>%
    dplyr::mutate(sample = sample_full) %>%
    tidyr::separate(sample, into = c("sex", "stage", "rep"), sep = "_", extra = "merge", fill = "left") %>%
    dplyr::mutate(
      sex = ifelse(is.na(sex), "U", sex),
      stage = dplyr::case_when(
        stage == "Mira" ~ "Miracidia",
        stage == "1d" ~ "Sporocysts_1d",
        stage == "5d" ~ "Sporocysts_5d",
        stage == "32d" ~ "Sporocysts_32d",
        stage == "Cerc" ~ "Cercariae",
        stage == "2d" ~ "Somules_2d",
        stage == "26d" ~ "Juveniles",
        TRUE ~ stage
      ),
      rep = stringr::str_remove(rep, "^R"),
      logTPM = log2(TPM + 1),
      z_score = as.numeric(scale(logTPM)),
      med_norm = logTPM - median(logTPM, na.rm = TRUE)
    )

  gene_data <- map_stages(gene_data)
  return(gene_data)
}

map_stages <- function(df) {
  stage_map <- c(
    "Eggs" = "Ovos",
    "Miracidia" = "Miracídio",
    "Sporocysts_1d" = "Esporocisto_1d",
    "Sporocysts_5d" = "Esporocisto_5d",
    "Sporocysts_32d" = "Esporocisto_32d",
    "Cercariae" = "Cercária",
    "Somules_2d" = "Esquistossômulo_2d",
    "Juveniles" = "Juvenis_26d"
  )
  df %>%
    dplyr::mutate(stage_category = dplyr::recode(stage, !!!stage_map, .default = stage)) %>%
    dplyr::mutate(stage_category = factor(stage_category, levels = unname(stage_map)))
}

compute_stats <- function(df) {
  df %>%
    dplyr::group_by(stage_category) %>%
    dplyr::summarise(
      mean_tpm = mean(TPM, na.rm = TRUE),
      sd_tpm = sd(TPM, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
}

# Graphs

plot_boxplot <- function(df, gene, outfile) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = stage_category, y = logTPM, fill = stage_category)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.6) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = paste("Expressão de", gene),
         x = "Estágio", y = "log2(TPM+1)")
  ggplot2::ggsave(outfile, p, width = 10, height = 7, dpi = 300)
}

plot_barplot <- function(stats, gene, outfile) {
  p <- ggplot2::ggplot(stats, ggplot2::aes(x = stage_category, y = mean_tpm, fill = stage_category)) +
    ggplot2::geom_col(color = "black") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_tpm - sd_tpm, ymax = mean_tpm + sd_tpm), width = 0.2) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = paste("TPM médio de", gene),
         x = "Estágio", y = "TPM médio")
  ggplot2::ggsave(outfile, p, width = 10, height = 7, dpi = 300)
}

plot_lineplot <- function(df, gene, outfile) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = stage_category, y = logTPM, group = sex, color = sex)) +
    ggplot2::stat_summary(fun = mean, geom = "line", linewidth = 1) +
    ggplot2::stat_summary(fun = mean, geom = "point", size = 3) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = paste("Expressão média de", gene, "ao longo do ciclo de vida"),
         x = "Estágio", y = "log2(TPM+1)", color = "Sexo")
  ggplot2::ggsave(outfile, p, width = 10, height = 7, dpi = 300)
}

plot_box_sex <- function(df, gene, outfile) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = stage_category, y = logTPM, fill = sex)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = paste("TPM por sexo e estágio -", gene),
         x = "Estágio", y = "log2(TPM+1)", fill = "Sexo")
  ggplot2::ggsave(outfile, p, width = 10, height = 7, dpi = 300)
}

# Violins

plot_violin_log2 <- function(df, gene, outfile) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = stage_category, y = logTPM, fill = stage_category)) +
    ggplot2::geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.5) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = paste("Distribuição log2(TPM+1) -", gene),
         x = "Estágio", y = "log2(TPM+1)")
  ggplot2::ggsave(outfile, p, width = 10, height = 7, dpi = 300)
}

plot_violin_raw <- function(df, gene, outfile) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = stage_category, y = TPM, fill = stage_category)) +
    ggplot2::geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.5) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = paste("Distribuição TPM Bruto -", gene),
         x = "Estágio", y = "TPM")
  ggplot2::ggsave(outfile, p, width = 10, height = 7, dpi = 300)
}

plot_violin_zscore <- function(df, gene, outfile) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = stage_category, y = z_score, fill = stage_category)) +
    ggplot2::geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.5) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = paste("Distribuição Z-Score -", gene),
         x = "Estágio", y = "Z-Score (log2TPM)")
  ggplot2::ggsave(outfile, p, width = 10, height = 7, dpi = 300)
}

plot_violin_median <- function(df, gene, outfile) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = stage_category, y = med_norm, fill = stage_category)) +
    ggplot2::geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.5) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = paste("Distribuição Normalizada pela Mediana -", gene),
         x = "Estágio", y = "log2TPM - Mediana")
  ggplot2::ggsave(outfile, p, width = 10, height = 7, dpi = 300)
}

plot_violin_fixed <- function(df, gene, outfile) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = stage_category, y = logTPM, fill = stage_category)) +
    ggplot2::geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.5) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = paste("Distribuição (Escala Fixa) -", gene),
         x = "Estágio", y = "log2(TPM+1)")
  ggplot2::ggsave(outfile, p, width = 10, height = 7, dpi = 300)
}

# Heatmaps

plot_heatmap_sample <- function(df, gene, outfile) {

  # Matriz logTPM usando sample_full direto (garantia de unicidade)
  heat_data <- df %>%
    dplyr::select(sample_full, logTPM) %>%
    tibble::column_to_rownames("sample_full") %>%
    t()

  # Annotation data
  annotation_col <- df %>%
    dplyr::select(sample_full, sex, stage_category) %>%
    distinct() %>%
    tibble::column_to_rownames("sample_full")

  color_breaks <- seq(0, 10, length.out = 100)
  colors <- colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(99)

  pheatmap::pheatmap(
    heat_data,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    annotation_col = annotation_col,
    color = colors,
    breaks = color_breaks,
    border_color = NA,
    main = paste("Expressão por amostra\n", gene, "(log2 TPM+1)"),
    filename = outfile,
    width = 12,
    height = 4
  )
}


plot_heatmap_stage_mean <- function(df, gene, outfile) {
  stage_means <- df %>%
    dplyr::group_by(stage_category) %>%
    dplyr::summarise(mean_logTPM = mean(logTPM, na.rm = TRUE), .groups = "drop")

  z_values <- scale(stage_means$mean_logTPM)

  heat_matrix <- matrix(
    z_values,
    ncol = 1
  )

  rownames(heat_matrix) <- stage_means$stage_category
  colnames(heat_matrix) <- gene

  max_abs <- max(abs(z_values), na.rm = TRUE)
  breaks <- seq(-max_abs, max_abs, length.out = 100)
  colors <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(99)

  pheatmap::pheatmap(
    heat_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colors,
    breaks = breaks,
    border_color = NA,
    main = paste("Expressão relativa por estágio\n", gene, "(Z-score)"),
    filename = outfile,
    width = 4,
    height = 8
  )
}
