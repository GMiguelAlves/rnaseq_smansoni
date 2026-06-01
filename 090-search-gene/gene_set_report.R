#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
  library(rtracklayer)
})

filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
summarise <- dplyr::summarise
distinct <- dplyr::distinct
left_join <- dplyr::left_join
bind_rows <- dplyr::bind_rows
bind_cols <- dplyr::bind_cols
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
arrange <- dplyr::arrange
rowwise <- dplyr::rowwise
case_when <- dplyr::case_when
n_distinct <- dplyr::n_distinct
rename <- dplyr::rename
pull <- dplyr::pull
coalesce <- dplyr::coalesce

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = "") {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

log_info <- function(msg) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")

genes_file <- get_arg("--genes", "genes.txt")
tpm_file <- get_arg("--tpm", file.path(Sys.getenv("QUANTIFICATION_DIR", unset = "../050-quantification"), "tpm_matrix.tsv"))
samples_file <- get_arg("--samples", file.path(Sys.getenv("QUANTIFICATION_DIR", unset = "../050-quantification"), "quant_samples.tsv"))
metadata_file <- get_arg("--metadata", Sys.getenv("METADATA_FINAL_NEW", unset = Sys.getenv("METADATA_FINAL", unset = "")))
deg_root <- get_arg("--deg-root", Sys.getenv("DEG_DIR", unset = "../060-deg-analysis"))
gff_file <- get_arg("--gff", Sys.getenv("REF_GFF3", unset = ""))
out_dir <- get_arg("--output-dir", file.path(Sys.getenv("GENE_REPORT_DIR", unset = "."), "results"))
report_title <- get_arg("--title", "Relatorio de genes candidatos")

if (!file.exists(genes_file)) stop("[ERRO] genes.txt nao encontrado: ", genes_file)
if (!file.exists(tpm_file)) stop("[ERRO] Matriz TPM nao encontrada: ", tpm_file)
if (!file.exists(samples_file)) stop("[ERRO] Tabela de amostras nao encontrada: ", samples_file)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "genes"), recursive = TRUE, showWarnings = FALSE)

sanitize <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "unknown"
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[x == ""] <- "unknown"
  x
}

write_tsv2 <- function(df, path) {
  readr::write_tsv(df, path, na = "")
}

parse_gene_groups <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  rows <- list()
  for (line in lines) {
    line <- trimws(line)
    if (line == "" || startsWith(line, "#")) next
    if (!grepl(":", line, fixed = TRUE)) {
      warning("Linha ignorada em genes.txt sem ':': ", line)
      next
    }
    parts <- strsplit(line, ":", fixed = TRUE)[[1]]
    group <- trimws(parts[1])
    genes <- trimws(unlist(strsplit(paste(parts[-1], collapse = ":"), "[,;]")))
    genes <- genes[genes != ""]
    if (length(genes) == 0) next
    rows[[length(rows) + 1]] <- data.frame(group = group, query = genes, stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) stop("[ERRO] Nenhum gene encontrado em ", path)
  bind_rows(rows) %>% distinct(group, query, .keep_all = TRUE)
}

read_matrix <- function(path) {
  df <- readr::read_tsv(path, show_col_types = FALSE, col_types = cols(.default = col_character()))
  if (ncol(df) < 2) stop("[ERRO] Matriz invalida: ", path)
  gene_col <- colnames(df)[1]
  colnames(df)[1] <- "gene_id"
  df %>%
    mutate(across(-gene_id, ~ as.numeric(.x)))
}

read_samples <- function(path, sample_names) {
  samples <- readr::read_tsv(path, show_col_types = FALSE, col_types = cols(.default = col_character()))
  if (!"import_id" %in% colnames(samples)) {
    if (all(c("dataset", "sample_id") %in% colnames(samples))) {
      combined <- paste(samples$dataset, samples$sample_id, sep = "__")
      if (all(sample_names %in% combined)) {
        samples$import_id <- combined
      } else {
        samples$import_id <- samples$sample_id
      }
    } else if ("sample_id" %in% colnames(samples)) {
      samples$import_id <- samples$sample_id
    } else {
      stop("[ERRO] Tabela de amostras precisa de import_id ou sample_id.")
    }
  }
  samples <- samples %>% distinct(import_id, .keep_all = TRUE)
  missing <- setdiff(sample_names, samples$import_id)
  if (length(missing) > 0) stop("[ERRO] Amostras sem metadata: ", paste(head(missing, 20), collapse = ", "))
  samples[match(sample_names, samples$import_id), , drop = FALSE]
}

load_annotations <- function(gff_file) {
  if (gff_file == "" || !file.exists(gff_file)) {
    return(tibble(gene_id = character(), gene_name = character(), biotype = character(), description = character()))
  }
  log_info("Lendo anotacao GFF3...")
  gff <- rtracklayer::import(gff_file)
  genes <- as.data.frame(gff[gff$type == "gene"])
  if (nrow(genes) == 0) {
    return(tibble(gene_id = character(), gene_name = character(), biotype = character(), description = character()))
  }
  pick_col <- function(df, names, default = NA_character_) {
    found <- intersect(names, colnames(df))
    if (length(found) == 0) return(rep(default, nrow(df)))
    as.character(df[[found[1]]])
  }
  gene_id <- pick_col(genes, c("ID", "gene_id"))
  gene_id <- gsub("^gene:", "", gene_id)
  gene_id <- gsub("\\.[0-9]+$", "", gene_id)
  gene_name <- pick_col(genes, c("Name", "gene_name", "symbol"))
  gene_name[is.na(gene_name) | gene_name == ""] <- gene_id[is.na(gene_name) | gene_name == ""]
  biotype <- pick_col(genes, c("biotype", "gene_biotype", "type"), "Unknown")
  biotype[is.na(biotype) | biotype == ""] <- "Unknown"
  description <- pick_col(genes, c("description", "product", "Note", "note"), "")
  tibble(gene_id = gene_id, gene_name = gene_name, biotype = biotype, description = description) %>%
    distinct(gene_id, .keep_all = TRUE)
}

build_gene_catalog <- function(gene_groups, tpm, annotations) {
  tpm_genes <- tpm$gene_id
  ann <- annotations
  rows <- gene_groups %>%
    rowwise() %>%
    mutate(
      matched_gene_id = case_when(
        query %in% tpm_genes ~ query,
        query %in% ann$gene_id ~ query,
        query %in% ann$gene_name ~ ann$gene_id[match(query, ann$gene_name)],
        TRUE ~ query
      ),
      match_type = case_when(
        query %in% tpm_genes ~ "gene_id",
        query %in% ann$gene_id ~ "annotation_gene_id",
        query %in% ann$gene_name ~ "gene_name",
        TRUE ~ "unmatched"
      )
    ) %>%
    ungroup() %>%
    left_join(ann, by = c("matched_gene_id" = "gene_id")) %>%
    mutate(
      gene_name = ifelse(is.na(gene_name) | gene_name == "", matched_gene_id, gene_name),
      biotype = ifelse(is.na(biotype) | biotype == "", "Unknown", biotype),
      description = ifelse(is.na(description), "", description),
      found_in_tpm = matched_gene_id %in% tpm_genes
    ) %>%
    distinct(group, query, matched_gene_id, .keep_all = TRUE)
  rows
}

load_deg_hits <- function(deg_root, gene_catalog) {
  files <- list.files(deg_root, pattern = "DEGs_all_results.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(tibble())
  rows <- lapply(files, function(path) {
    df <- tryCatch(readr::read_tsv(path, show_col_types = FALSE, col_types = cols(.default = col_character())), error = function(e) NULL)
    if (is.null(df) || !"gene_id" %in% colnames(df)) return(NULL)
    rel <- gsub("\\\\", "/", sub(paste0("^", normalizePath(deg_root, winslash = "/", mustWork = FALSE), "/?"), "", normalizePath(path, winslash = "/", mustWork = FALSE)))
    df %>%
      filter(gene_id %in% gene_catalog$matched_gene_id) %>%
      mutate(
        source_file = rel,
        result_dir = dirname(rel),
        padj_num = suppressWarnings(as.numeric(padj)),
        log2FoldChange_num = suppressWarnings(as.numeric(log2FoldChange)),
        significant = !is.na(padj_num) & padj_num < 0.05 & abs(log2FoldChange_num) >= 1
      )
  })
  bind_rows(rows)
}

plot_or_skip <- function(expr, plot_fun) {
  tryCatch(plot_fun(), error = function(e) {
    warning(expr, ": ", e$message)
    FALSE
  })
}

make_expression_long <- function(tpm, samples, gene_catalog) {
  selected <- tpm %>% filter(gene_id %in% gene_catalog$matched_gene_id)
  selected %>%
    pivot_longer(-gene_id, names_to = "import_id", values_to = "TPM") %>%
    left_join(samples, by = "import_id") %>%
    left_join(gene_catalog %>% select(group, query, matched_gene_id, gene_name, biotype), by = c("gene_id" = "matched_gene_id")) %>%
    mutate(
      TPM = as.numeric(TPM),
      log2TPM = log2(TPM + 1),
      dataset = ifelse(is.na(dataset) | dataset == "", "unknown", dataset),
      stage = ifelse(is.na(stage) | stage == "", "unknown", stage),
      tissue = ifelse(is.na(tissue) | tissue == "", "unknown", tissue),
      sex = ifelse(is.na(sex) | sex == "", "unknown", sex),
      condition = ifelse(is.na(condition) | condition == "", "unknown", condition)
    ) %>%
    group_by(gene_id) %>%
    mutate(z_log2TPM = as.numeric(scale(log2TPM))) %>%
    ungroup()
}

plot_group_heatmap <- function(summary_df, outfile) {
  mat <- summary_df %>%
    unite("context", dataset, stage, tissue, sex, remove = FALSE) %>%
    select(gene_id, gene_name, context, mean_log2TPM) %>%
    mutate(label = paste(gene_id, gene_name, sep = " | ")) %>%
    select(label, context, mean_log2TPM) %>%
    distinct() %>%
    pivot_wider(names_from = context, values_from = mean_log2TPM, values_fill = 0)
  if (nrow(mat) == 0 || ncol(mat) < 2) return(FALSE)
  matrix_data <- as.matrix(mat[, -1, drop = FALSE])
  rownames(matrix_data) <- mat$label
  pheatmap::pheatmap(matrix_data, scale = "row", border_color = NA,
                     fontsize_row = 7, fontsize_col = 7,
                     main = "Genes candidatos - log2(TPM+1) medio",
                     filename = outfile, width = 12, height = max(5, min(14, nrow(matrix_data) * 0.28 + 3)))
  TRUE
}

plot_dotplot <- function(summary_df, outfile) {
  df <- summary_df %>%
    mutate(context = paste(dataset, stage, sep = " | "),
           gene_label = paste(group, gene_name, sep = " - "))
  if (nrow(df) == 0) return(FALSE)
  p <- ggplot(df, aes(x = context, y = gene_label)) +
    geom_point(aes(size = fraction_expressed, color = mean_log2TPM), alpha = 0.85) +
    scale_color_viridis_c(option = "C") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Projeto | estagio", y = "Grupo - gene", color = "Media log2TPM", size = "Frac. TPM>1")
  ggsave(outfile, p, width = 12, height = max(5, min(16, length(unique(df$gene_label)) * 0.28 + 3)), dpi = 300)
  TRUE
}

plot_gene_profiles <- function(expr_long, out_dir) {
  genes <- unique(expr_long$gene_id)
  for (gene in genes) {
    df <- expr_long %>% filter(gene_id == gene)
    gene_dir <- file.path(out_dir, "genes", sanitize(gene))
    dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)
    label <- unique(paste(df$gene_id, df$gene_name, sep = " | "))[1]
    p1 <- ggplot(df, aes(x = stage, y = log2TPM, fill = dataset)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.65) +
      geom_jitter(aes(color = dataset), width = 0.18, alpha = 0.55, size = 1.6) +
      facet_grid(tissue ~ sex, scales = "free_x", space = "free_x") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = label, x = "Estagio", y = "log2(TPM+1)")
    ggsave(file.path(gene_dir, "expression_boxplot_stage_tissue_sex.png"), p1, width = 13, height = 8, dpi = 300)

    p2 <- df %>%
      group_by(dataset, stage, tissue, sex) %>%
      summarise(mean_log2TPM = mean(log2TPM, na.rm = TRUE), .groups = "drop") %>%
      ggplot(aes(x = stage, y = mean_log2TPM, color = sex, group = interaction(dataset, tissue, sex))) +
      geom_line(alpha = 0.75) +
      geom_point(size = 2) +
      facet_grid(dataset ~ tissue, scales = "free_x", space = "free_x") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste("Perfil medio:", label), x = "Estagio", y = "Media log2(TPM+1)")
    ggsave(file.path(gene_dir, "expression_profile_line.png"), p2, width = 13, height = 8, dpi = 300)
  }
}

plot_deg_heatmap <- function(deg_hits, gene_catalog, outfile) {
  if (nrow(deg_hits) == 0) return(FALSE)
  df <- deg_hits %>%
    mutate(contrast_label = paste(result_dir, contrast, sep = " | ")) %>%
    select(gene_id, contrast_label, log2FoldChange_num) %>%
    distinct() %>%
    pivot_wider(names_from = contrast_label, values_from = log2FoldChange_num, values_fill = 0) %>%
    left_join(gene_catalog %>% select(matched_gene_id, gene_name, group), by = c("gene_id" = "matched_gene_id")) %>%
    mutate(label = paste(group, gene_name, gene_id, sep = " | "))
  if (nrow(df) == 0 || ncol(df) <= 4) return(FALSE)
  mat <- as.matrix(df[, setdiff(colnames(df), c("gene_id", "gene_name", "group", "label")), drop = FALSE])
  rownames(mat) <- df$label
  pheatmap::pheatmap(mat, color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(101),
                     border_color = NA, fontsize_row = 7, fontsize_col = 6,
                     main = "log2FC dos genes candidatos em contrastes DEG",
                     filename = outfile, width = 14, height = max(5, min(14, nrow(mat) * 0.3 + 3)))
  TRUE
}

plot_deg_bar <- function(deg_hits, gene_catalog, outfile) {
  if (nrow(deg_hits) == 0) return(FALSE)
  df <- deg_hits %>%
    group_by(gene_id) %>%
    summarise(n_contrasts = n_distinct(contrast), n_significant = sum(significant, na.rm = TRUE), .groups = "drop") %>%
    left_join(gene_catalog %>% select(matched_gene_id, gene_name, group), by = c("gene_id" = "matched_gene_id")) %>%
    mutate(label = paste(group, gene_name, sep = " - "))
  p <- ggplot(df, aes(x = reorder(label, n_significant), y = n_significant, fill = group)) +
    geom_col() +
    coord_flip() +
    theme_bw(base_size = 10) +
    labs(x = "Gene", y = "Contrastes DEG significativos", fill = "Grupo")
  ggsave(outfile, p, width = 10, height = max(5, min(14, nrow(df) * 0.28 + 3)), dpi = 300)
  TRUE
}

parse_deg_context <- function(deg_hits) {
  if (nrow(deg_hits) == 0) return(deg_hits)
  deg_hits %>%
    mutate(
      result_dir = ifelse(is.na(result_dir) | result_dir == ".", "unknown", result_dir),
      deg_project = sub("/.*$", "", result_dir),
      deg_mode = ifelse(grepl("/", result_dir), sub("^.*/", "", result_dir), "unknown"),
      contrast_label = paste(result_dir, contrast, sep = " | ")
    )
}

top_context_value <- function(df, value_col, context_col) {
  if (nrow(df) == 0) return(tibble(gene_id = character(), value = numeric(), context = character()))
  df %>%
    group_by(gene_id, .data[[context_col]]) %>%
    summarise(value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    arrange(gene_id, desc(value)) %>%
    group_by(gene_id) %>%
    summarise(value = dplyr::first(value), context = dplyr::first(.data[[context_col]]), .groups = "drop")
}

compute_gene_scores <- function(expr_long, expr_summary, deg_hits, gene_catalog) {
  deg_hits <- parse_deg_context(deg_hits)

  expr_gene <- expr_long %>%
    group_by(group, gene_id, gene_name) %>%
    summarise(
      mean_TPM_all = mean(TPM, na.rm = TRUE),
      mean_log2TPM_all = mean(log2TPM, na.rm = TRUE),
      max_TPM_sample = max(TPM, na.rm = TRUE),
      fraction_samples_TPM_gt1 = mean(TPM > 1, na.rm = TRUE),
      n_samples_TPM_gt1 = sum(TPM > 1, na.rm = TRUE),
      .groups = "drop"
    )

  context_expr <- expr_summary %>%
    group_by(gene_id) %>%
    summarise(
      max_mean_log2TPM_context = max(mean_log2TPM, na.rm = TRUE),
      n_contexts_TPM_gt1 = sum(mean_TPM > 1, na.rm = TRUE),
      n_contexts = dplyr::n(),
      .groups = "drop"
    )

  dominant_stage <- top_context_value(expr_summary, "mean_log2TPM", "stage") %>%
    rename(dominant_stage_mean_log2TPM = value, dominant_stage = context)
  dominant_tissue <- top_context_value(expr_summary, "mean_log2TPM", "tissue") %>%
    rename(dominant_tissue_mean_log2TPM = value, dominant_tissue = context)
  dominant_sex <- top_context_value(expr_summary, "mean_log2TPM", "sex") %>%
    rename(dominant_sex_mean_log2TPM = value, dominant_sex = context)

  tissue_means <- expr_long %>%
    group_by(gene_id, tissue) %>%
    summarise(tissue_mean_log2TPM = mean(log2TPM, na.rm = TRUE), .groups = "drop")
  tissue_scores <- tissue_means %>%
    group_by(gene_id) %>%
    summarise(
      ovary_mean_log2TPM = ifelse(any(tissue == "ovary"), max(tissue_mean_log2TPM[tissue == "ovary"], na.rm = TRUE), NA_real_),
      testis_mean_log2TPM = ifelse(any(tissue == "testis"), max(tissue_mean_log2TPM[tissue == "testis"], na.rm = TRUE), NA_real_),
      non_gonad_mean_log2TPM = ifelse(any(!tissue %in% c("ovary", "testis")), mean(tissue_mean_log2TPM[!tissue %in% c("ovary", "testis")], na.rm = TRUE), NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      ovary_vs_testis_log2TPM = ovary_mean_log2TPM - testis_mean_log2TPM,
      gonad_specificity_log2TPM = pmax(ovary_mean_log2TPM, testis_mean_log2TPM, na.rm = TRUE) - non_gonad_mean_log2TPM
    )

  sex_means <- expr_long %>%
    group_by(gene_id, sex) %>%
    summarise(sex_mean_log2TPM = mean(log2TPM, na.rm = TRUE), .groups = "drop")
  sex_scores <- sex_means %>%
    group_by(gene_id) %>%
    summarise(
      female_mean_log2TPM = ifelse(any(sex == "female"), max(sex_mean_log2TPM[sex == "female"], na.rm = TRUE), NA_real_),
      male_mean_log2TPM = ifelse(any(sex == "male"), max(sex_mean_log2TPM[sex == "male"], na.rm = TRUE), NA_real_),
      .groups = "drop"
    ) %>%
    mutate(female_vs_male_log2TPM = female_mean_log2TPM - male_mean_log2TPM)

  if (nrow(deg_hits) > 0) {
    deg_summary <- deg_hits %>%
      group_by(gene_id) %>%
      summarise(
        n_deg_hits = dplyr::n(),
        n_deg_contrasts = n_distinct(contrast_label),
        n_significant_contrasts = sum(significant, na.rm = TRUE),
        n_projects_significant = n_distinct(deg_project[significant]),
        n_modes_significant = n_distinct(deg_mode[significant]),
        max_abs_log2FC = suppressWarnings(max(abs(log2FoldChange_num), na.rm = TRUE)),
        best_padj = suppressWarnings(min(padj_num, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      mutate(
        max_abs_log2FC = ifelse(is.infinite(max_abs_log2FC), NA_real_, max_abs_log2FC),
        best_padj = ifelse(is.infinite(best_padj), NA_real_, best_padj)
      )
  } else {
    deg_summary <- tibble(gene_id = gene_catalog$matched_gene_id) %>%
      mutate(
        n_deg_hits = 0L,
        n_deg_contrasts = 0L,
        n_significant_contrasts = 0L,
        n_projects_significant = 0L,
        n_modes_significant = 0L,
        max_abs_log2FC = NA_real_,
        best_padj = NA_real_
      )
  }

  gene_catalog %>%
    select(group, query, matched_gene_id, gene_name, biotype, found_in_tpm) %>%
    rename(gene_id = matched_gene_id) %>%
    left_join(expr_gene, by = c("group", "gene_id", "gene_name")) %>%
    left_join(context_expr, by = "gene_id") %>%
    left_join(dominant_stage, by = "gene_id") %>%
    left_join(dominant_tissue, by = "gene_id") %>%
    left_join(dominant_sex, by = "gene_id") %>%
    left_join(tissue_scores, by = "gene_id") %>%
    left_join(sex_scores, by = "gene_id") %>%
    left_join(deg_summary, by = "gene_id") %>%
    mutate(
      across(c(n_deg_hits, n_deg_contrasts, n_significant_contrasts, n_projects_significant, n_modes_significant), ~ ifelse(is.na(.x), 0, .x)),
      expression_score = coalesce(max_mean_log2TPM_context, 0) + coalesce(fraction_samples_TPM_gt1, 0) + log1p(coalesce(n_contexts_TPM_gt1, 0)),
      deg_consistency_score = coalesce(n_significant_contrasts, 0) + coalesce(n_projects_significant, 0) + coalesce(n_modes_significant, 0) + coalesce(max_abs_log2FC, 0),
      gonad_score = abs(coalesce(ovary_vs_testis_log2TPM, 0)) + pmax(coalesce(gonad_specificity_log2TPM, 0), 0),
      candidate_score = expression_score + deg_consistency_score + gonad_score
    ) %>%
    arrange(desc(candidate_score), group, gene_name)
}

plot_tissue_sex_heatmap <- function(expr_summary, outfile) {
  df <- expr_summary %>%
    mutate(context = paste(tissue, sex, sep = " | "),
           label = paste(group, gene_name, gene_id, sep = " | ")) %>%
    group_by(label, context) %>%
    summarise(mean_log2TPM = mean(mean_log2TPM, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = context, values_from = mean_log2TPM, values_fill = 0)
  if (nrow(df) == 0 || ncol(df) < 2) return(FALSE)
  mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- df$label
  pheatmap::pheatmap(mat, scale = "row", border_color = NA,
                     fontsize_row = 7, fontsize_col = 8,
                     main = "Padroes por tecido e sexo",
                     filename = outfile, width = 11, height = max(5, min(14, nrow(mat) * 0.3 + 3)))
  TRUE
}

plot_ovary_testis_scatter <- function(gene_scores, outfile) {
  df <- gene_scores %>%
    filter(!is.na(ovary_mean_log2TPM) | !is.na(testis_mean_log2TPM)) %>%
    mutate(
      ovary_mean_log2TPM = coalesce(ovary_mean_log2TPM, 0),
      testis_mean_log2TPM = coalesce(testis_mean_log2TPM, 0),
      label = paste(gene_name, gene_id, sep = " | ")
    )
  if (nrow(df) == 0) return(FALSE)
  p <- ggplot(df, aes(x = testis_mean_log2TPM, y = ovary_mean_log2TPM, color = group, size = fraction_samples_TPM_gt1)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray55") +
    geom_point(alpha = 0.85) +
    geom_text(aes(label = gene_name), check_overlap = TRUE, vjust = -0.8, size = 3) +
    theme_bw(base_size = 11) +
    labs(x = "Testiculo: media log2(TPM+1)", y = "Ovario: media log2(TPM+1)", color = "Grupo", size = "Frac. TPM>1")
  ggsave(outfile, p, width = 8, height = 7, dpi = 300)
  TRUE
}

plot_consistency_tile <- function(deg_hits, gene_catalog, outfile) {
  if (nrow(deg_hits) == 0) return(FALSE)
  df <- parse_deg_context(deg_hits) %>%
    left_join(gene_catalog %>% select(matched_gene_id, gene_name, group), by = c("gene_id" = "matched_gene_id")) %>%
    mutate(
      gene_label = paste(group, gene_name, gene_id, sep = " | "),
      contrast_short = paste(result_dir, contrast, sep = " | "),
      sig_label = ifelse(significant, "significativo", "nao_significativo")
    )
  if (nrow(df) == 0) return(FALSE)
  if (length(unique(df$contrast_short)) > 80) {
    keep <- df %>%
      group_by(contrast_short) %>%
      summarise(best = min(padj_num, na.rm = TRUE), .groups = "drop") %>%
      arrange(best) %>%
      head(80) %>%
      pull(contrast_short)
    df <- df %>% filter(contrast_short %in% keep)
  }
  p <- ggplot(df, aes(x = contrast_short, y = gene_label, fill = log2FoldChange_num, alpha = sig_label)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", na.value = "gray90") +
    scale_alpha_manual(values = c("significativo" = 1, "nao_significativo" = 0.35)) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank()) +
    labs(x = "Projeto/modo | contraste", y = "Gene", fill = "log2FC", alpha = "")
  ggsave(outfile, p, width = 15, height = max(5, min(16, length(unique(df$gene_label)) * 0.35 + 3)), dpi = 300)
  TRUE
}

plot_expression_vs_deg <- function(gene_scores, outfile) {
  df <- gene_scores %>%
    mutate(label = paste(gene_name, gene_id, sep = " | "))
  if (nrow(df) == 0) return(FALSE)
  p <- ggplot(df, aes(x = max_mean_log2TPM_context, y = n_significant_contrasts, color = group, size = n_projects_significant)) +
    geom_point(alpha = 0.85) +
    geom_text(aes(label = gene_name), check_overlap = TRUE, vjust = -0.7, size = 3) +
    theme_bw(base_size = 11) +
    labs(x = "Maior media log2(TPM+1) em um contexto", y = "Contrastes DEG significativos", color = "Grupo", size = "Projetos")
  ggsave(outfile, p, width = 9, height = 7, dpi = 300)
  TRUE
}

plot_candidate_score_bar <- function(gene_scores, outfile) {
  df <- gene_scores %>%
    arrange(desc(candidate_score)) %>%
    head(40) %>%
    mutate(label = paste(group, gene_name, sep = " - "))
  if (nrow(df) == 0) return(FALSE)
  p <- ggplot(df, aes(x = reorder(label, candidate_score), y = candidate_score, fill = group)) +
    geom_col() +
    coord_flip() +
    theme_bw(base_size = 10) +
    labs(x = "Gene", y = "Score candidato", fill = "Grupo")
  ggsave(outfile, p, width = 10, height = max(5, min(14, nrow(df) * 0.28 + 3)), dpi = 300)
  TRUE
}

html_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

table_to_html <- function(df, max_rows = 20) {
  if (is.null(df) || nrow(df) == 0) return("<p><em>Nenhum registro.</em></p>")
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(df) > max_rows) df <- df[seq_len(max_rows), , drop = FALSE]
  header <- paste0("<tr>", paste0("<th>", html_escape(colnames(df)), "</th>", collapse = ""), "</tr>")
  rows <- apply(df, 1, function(row) {
    paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>")
  })
  paste0("<table>", header, paste(rows, collapse = "\n"), "</table>")
}

write_html_report <- function(path, title, catalog, summary_df, deg_hits, gene_scores, plot_files) {
  img_tag <- function(file, caption) {
    if (!file.exists(file)) return("")
    rel <- gsub("\\\\", "/", file.path("plots", basename(file)))
    paste0("<figure><img src='", rel, "' style='max-width:100%;border:1px solid #ddd'><figcaption>", caption, "</figcaption></figure>")
  }
  section_img <- function(file, caption) {
    x <- img_tag(file, caption)
    if (x == "") paste0("<p><em>", caption, " nao foi gerado.</em></p>") else x
  }
  top_candidates <- gene_scores %>%
    select(group, gene_id, gene_name, candidate_score, expression_score, deg_consistency_score, gonad_score,
           n_significant_contrasts, n_projects_significant, max_mean_log2TPM_context,
           fraction_samples_TPM_gt1, dominant_stage, dominant_tissue, dominant_sex) %>%
    arrange(desc(candidate_score))
  top_deg <- gene_scores %>%
    select(group, gene_id, gene_name, n_significant_contrasts, n_projects_significant, n_modes_significant,
           max_abs_log2FC, best_padj, candidate_score) %>%
    arrange(desc(n_significant_contrasts), desc(n_projects_significant), best_padj)
  top_gonad <- gene_scores %>%
    select(group, gene_id, gene_name, ovary_mean_log2TPM, testis_mean_log2TPM, ovary_vs_testis_log2TPM,
           gonad_specificity_log2TPM, dominant_tissue, female_vs_male_log2TPM) %>%
    arrange(desc(abs(ovary_vs_testis_log2TPM)), desc(gonad_specificity_log2TPM))
  top_expression <- gene_scores %>%
    select(group, gene_id, gene_name, max_mean_log2TPM_context, mean_log2TPM_all, fraction_samples_TPM_gt1,
           n_contexts_TPM_gt1, dominant_stage, dominant_tissue, dominant_sex) %>%
    arrange(desc(max_mean_log2TPM_context), desc(fraction_samples_TPM_gt1))

  gene_sections <- paste(vapply(unique(catalog$matched_gene_id), function(gene) {
    gene_dir <- file.path("genes", sanitize(gene))
    gene_row <- catalog[catalog$matched_gene_id == gene, , drop = FALSE][1, , drop = FALSE]
    score_row <- gene_scores[gene_scores$gene_id == gene, , drop = FALSE]
    score_text <- if (nrow(score_row) > 0) {
      paste0(
        "<p><b>Score candidato:</b> ", round(score_row$candidate_score[1], 3),
        " | <b>DEG significativos:</b> ", score_row$n_significant_contrasts[1],
        " | <b>Projetos:</b> ", score_row$n_projects_significant[1],
        " | <b>Expressao max:</b> ", round(score_row$max_mean_log2TPM_context[1], 3),
        "</p>"
      )
    } else {
      ""
    }
    paste0(
      "<h3>", gene_row$gene_name, " <code>", gene, "</code></h3>",
      "<p><b>Grupo:</b> ", gene_row$group, " | <b>Biotipo:</b> ", gene_row$biotype, "</p>",
      "<p>", gene_row$description, "</p>",
      score_text,
      "<img src='", gene_dir, "/expression_boxplot_stage_tissue_sex.png' style='max-width:100%;border:1px solid #ddd'>",
      "<img src='", gene_dir, "/expression_profile_line.png' style='max-width:100%;border:1px solid #ddd'>"
    )
  }, character(1)), collapse = "\n")
  html <- c(
    "<!doctype html><html><head><meta charset='utf-8'>",
    paste0("<title>", title, "</title>"),
    "<style>
      body{font-family:Arial,sans-serif;max-width:1280px;margin:32px auto;line-height:1.45;color:#222}
      nav{position:sticky;top:0;background:#fff;border-bottom:1px solid #ddd;padding:10px 0;margin-bottom:24px;z-index:2}
      nav a{margin-right:16px;color:#1d4e89;text-decoration:none;font-weight:600}
      h1,h2{color:#17324d} h2{border-top:2px solid #e6e6e6;padding-top:22px;margin-top:34px}
      .cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:12px;margin:18px 0}
      .card{background:#f7f9fb;border:1px solid #dde5ed;border-radius:6px;padding:14px}
      .card .num{font-size:28px;font-weight:700;color:#17324d}
      table{border-collapse:collapse;width:100%;font-size:12px;margin:10px 0 22px 0}
      th,td{border:1px solid #ddd;padding:5px;vertical-align:top} th{background:#f3f3f3}
      code{background:#f6f6f6;padding:2px 4px}
      figure{margin:18px 0 28px 0} figcaption{font-size:13px;color:#555;margin-top:6px}
      img{margin:8px 0 18px 0}
    </style>",
    "</head><body>",
    paste0("<h1>", title, "</h1>"),
    "<nav><a href='#ranking'>Ranking</a><a href='#expressao'>Expressao</a><a href='#deg'>DEG</a><a href='#gonadas'>Tecido/Sexo</a><a href='#genes'>Genes</a></nav>",
    "<div class='cards'>",
    paste0("<div class='card'><div class='num'>", nrow(catalog), "</div><div>genes consultados</div></div>"),
    paste0("<div class='card'><div class='num'>", sum(catalog$found_in_tpm), "</div><div>genes encontrados em TPM</div></div>"),
    paste0("<div class='card'><div class='num'>", sum(gene_scores$n_significant_contrasts > 0, na.rm = TRUE), "</div><div>genes com DEG significativo</div></div>"),
    paste0("<div class='card'><div class='num'>", sum(!is.na(gene_scores$ovary_vs_testis_log2TPM) & abs(gene_scores$ovary_vs_testis_log2TPM) > 1, na.rm = TRUE), "</div><div>genes com vies ovario/testiculo forte</div></div>"),
    "</div>",
    "<h2 id='ranking'>Ranking integrado</h2>",
    "<p>O score candidato combina expressao TPM robusta, consistencia DEG entre contrastes/projetos e sinais de especificidade gonadal.</p>",
    section_img(plot_files$candidate_score, "Ranking integrado de genes candidatos"),
    table_to_html(top_candidates, 30),
    "<h2>Resumo dos genes</h2>",
    paste0("<p>Tabelas completas em <code>tables/gene_catalog.tsv</code>, <code>tables/expression_summary.tsv</code> e <code>tables/deg_hits.tsv</code>.</p>"),
    paste0("<pre>", paste(capture.output(print(catalog %>% select(group, query, matched_gene_id, gene_name, biotype, found_in_tpm), n = Inf)), collapse = "\n"), "</pre>"),
    "<h2 id='expressao'>Expressao TPM clara</h2>",
    "<p>Esta secao prioriza genes com expressao clara nos TPM plots, independentemente de significancia estatistica.</p>",
    section_img(plot_files$heatmap, "Heatmap de expressao media por contexto"),
    section_img(plot_files$dotplot, "Dotplot de expressao media e fracao expressa"),
    section_img(plot_files$expression_vs_deg, "Comparacao entre expressao maxima e numero de contrastes DEG significativos"),
    table_to_html(top_expression, 30),
    "<h2 id='deg'>Top genes DEG e consistencia</h2>",
    "<p>Genes que aparecem repetidamente em mais de um contraste, projeto ou modo de analise tendem a ser candidatos mais robustos.</p>",
    section_img(plot_files$deg_heatmap, "Heatmap de log2FC em contrastes DEG"),
    section_img(plot_files$deg_bar, "Numero de contrastes DEG significativos por gene"),
    section_img(plot_files$consistency_tile, "Mapa de consistencia gene x contraste/projeto"),
    table_to_html(top_deg, 30),
    "<h2 id='gonadas'>Padroes de tecido e sexo</h2>",
    "<p>Foco especial em sinais de ovario/testiculo e vies por sexo, mantendo a informacao biologica preservada no metadata.</p>",
    section_img(plot_files$tissue_sex_heatmap, "Heatmap por tecido e sexo"),
    section_img(plot_files$ovary_testis, "Comparacao direta entre expressao em testiculo e ovario"),
    table_to_html(top_gonad, 30),
    "<h2 id='genes'>Genes individuais</h2>",
    gene_sections,
    "</body></html>"
  )
  writeLines(html, path, useBytes = TRUE)
}

gene_groups <- parse_gene_groups(genes_file)
tpm <- read_matrix(tpm_file)
samples <- read_samples(samples_file, setdiff(colnames(tpm), "gene_id"))
annotations <- load_annotations(gff_file)
gene_catalog <- build_gene_catalog(gene_groups, tpm, annotations)

if (metadata_file != "" && file.exists(metadata_file)) {
  metadata <- readr::read_csv(metadata_file, show_col_types = FALSE, col_types = cols(.default = col_character()))
  if (all(c("dataset", "sample_id") %in% colnames(metadata))) {
    metadata$import_id_combined <- paste(metadata$dataset, metadata$sample_id, sep = "__")
    key <- if (all(samples$import_id %in% metadata$import_id_combined)) "import_id_combined" else "sample_id"
    extra <- metadata[match(samples$import_id, metadata[[key]]), , drop = FALSE]
    add_cols <- setdiff(colnames(extra), colnames(samples))
    samples <- bind_cols(samples, extra[, add_cols, drop = FALSE])
  }
}

expr_long <- make_expression_long(tpm, samples, gene_catalog)
expr_summary <- expr_long %>%
  group_by(group, gene_id, gene_name, dataset, stage, tissue, sex, condition) %>%
  summarise(
    n = n(),
    mean_TPM = mean(TPM, na.rm = TRUE),
    median_TPM = median(TPM, na.rm = TRUE),
    mean_log2TPM = mean(log2TPM, na.rm = TRUE),
    fraction_expressed = mean(TPM > 1, na.rm = TRUE),
    .groups = "drop"
  )

deg_hits <- load_deg_hits(deg_root, gene_catalog)
deg_hits <- parse_deg_context(deg_hits)
gene_scores <- compute_gene_scores(expr_long, expr_summary, deg_hits, gene_catalog)

write_tsv2(gene_catalog, file.path(out_dir, "tables", "gene_catalog.tsv"))
write_tsv2(expr_long, file.path(out_dir, "tables", "expression_long.tsv"))
write_tsv2(expr_summary, file.path(out_dir, "tables", "expression_summary.tsv"))
write_tsv2(deg_hits, file.path(out_dir, "tables", "deg_hits.tsv"))
write_tsv2(gene_scores, file.path(out_dir, "tables", "gene_candidate_scores.tsv"))

plot_files <- list(
  heatmap = file.path(out_dir, "plots", "candidate_genes_expression_heatmap.png"),
  dotplot = file.path(out_dir, "plots", "candidate_genes_expression_dotplot.png"),
  deg_heatmap = file.path(out_dir, "plots", "candidate_genes_deg_log2fc_heatmap.png"),
  deg_bar = file.path(out_dir, "plots", "candidate_genes_deg_significant_barplot.png"),
  tissue_sex_heatmap = file.path(out_dir, "plots", "candidate_genes_tissue_sex_heatmap.png"),
  ovary_testis = file.path(out_dir, "plots", "candidate_genes_ovary_vs_testis.png"),
  consistency_tile = file.path(out_dir, "plots", "candidate_genes_deg_consistency_tile.png"),
  expression_vs_deg = file.path(out_dir, "plots", "candidate_genes_expression_vs_deg.png"),
  candidate_score = file.path(out_dir, "plots", "candidate_genes_integrated_score.png")
)

plot_or_skip("expression heatmap", function() plot_group_heatmap(expr_summary, plot_files$heatmap))
plot_or_skip("expression dotplot", function() plot_dotplot(expr_summary, plot_files$dotplot))
plot_or_skip("gene profiles", function() { plot_gene_profiles(expr_long, out_dir); TRUE })
plot_or_skip("DEG heatmap", function() plot_deg_heatmap(deg_hits, gene_catalog, plot_files$deg_heatmap))
plot_or_skip("DEG barplot", function() plot_deg_bar(deg_hits, gene_catalog, plot_files$deg_bar))
plot_or_skip("tissue/sex heatmap", function() plot_tissue_sex_heatmap(expr_summary, plot_files$tissue_sex_heatmap))
plot_or_skip("ovary/testis scatter", function() plot_ovary_testis_scatter(gene_scores, plot_files$ovary_testis))
plot_or_skip("DEG consistency tile", function() plot_consistency_tile(deg_hits, gene_catalog, plot_files$consistency_tile))
plot_or_skip("expression vs DEG", function() plot_expression_vs_deg(gene_scores, plot_files$expression_vs_deg))
plot_or_skip("candidate score", function() plot_candidate_score_bar(gene_scores, plot_files$candidate_score))

write_html_report(file.path(out_dir, "gene_set_report.html"), report_title, gene_catalog, expr_summary, deg_hits, gene_scores, plot_files)

log_info(paste("[OK] Relatorio 090 concluido:", file.path(out_dir, "gene_set_report.html")))
