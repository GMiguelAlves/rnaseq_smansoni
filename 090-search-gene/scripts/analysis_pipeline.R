#!/usr/bin/env Rscript

cat("=== INICIANDO PIPELINE DE ANÁLISE COMPLETO ===\n")

cat("Carregando pacotes...\n")

library(methods)
library(utils)
library(stats)
library(graphics)
library(grDevices)

#tidyr aux
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
})

#aux
suppressPackageStartupMessages({
  library(pheatmap)
  library(rstatix)
  library(jsonlite)
})

#bioconductor
suppressPackageStartupMessages({
  library(seqinr)
  library(rtracklayer)
  library(GO.db)
})

suppressPackageStartupMessages({
  library(Peptides)
})

cat("✓ Todos os pacotes carregados\n")

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
count <- dplyr::count

cat("✓ Conflitos de funções resolvidos\n")

# =====================
# Configurar caminhos
# =====================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Uso: Rscript analysis_pipeline.R <TPM_file> <gene> <project_dir>")
}

tpm_file <- args[1]
gene <- args[2]
project_dir <- args[3]

# Definir caminhos absolutos
scripts_dir <- file.path(project_dir, "090-search-gene", "scripts")
results_dir <- file.path(project_dir, "090-search-gene", "results", gene)
logs_dir <- file.path(project_dir, "090-search-gene", "logs")
main_results_dir <- file.path(project_dir, "090-search-gene", "results")

# Criar diretórios
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

cat("Carregando TODOS os módulos...\n")

module_files <- c(
  "config_analysis.R",
  "prem_utils.R", 
  "qc_utils.R",
  "logging_utils.R",
  "reporting_utils.R",
  "monitoring_utils.R"
)

for (module_file in module_files) {
  module_path <- file.path(scripts_dir, module_file)
  if (file.exists(module_path)) {
    source(module_path)
    cat("✓ Módulo carregado:", module_file, "\n")
  } else {
    cat("✗ Módulo não encontrado:", module_file, "\n")
  }
}

cat("✓ Todos os módulos carregados\n")

run_complete_gene_analysis <- function(tpm_file, gene, config) {
  
  # Setup
  out_dir <- file.path(config$results_dir, gene)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Iniciar sistema de logging
  logger <- create_logger(gene, config$log_dir)
  logger$info("=== INICIANDO ANÁLISE COMPLETA PARA: ", gene, " ===")
  
  results <- list()
  
  tryCatch({
    # FASE 1: CONTROLE DE QUALIDADE
    logger$info("FASE 1: Controle de qualidade dos dados")
    gene_data <- load_tpm_data(tpm_file, gene)
    results$qc <- perform_data_qc(gene_data, gene, config)
    generate_qc_report(results$qc, out_dir)
    
    if (results$qc$checks$gene_presence$status == "FAIL") {
      logger$warn("Gene não encontrado - pulando análises downstream")
      return(list(success = FALSE, reason = "Gene não encontrado"))
    }
    
    # Salvar dados brutos
    readr::write_csv(gene_data, file.path(out_dir, "TPM_data_long.csv"))
    
    # FASE 2: ANÁLISE DESCRITIVA
    logger$info("FASE 2: Análise descritiva")
    results$stats <- compute_stats(gene_data)
    readr::write_csv(results$stats, file.path(out_dir, "stats_by_category.csv"))
    
    # FASE 3: VISUALIZAÇÃO COMPLETA
    logger$info("FASE 3: Geração de visualizações completas")
    
    # Gerar TODOS os gráficos
    plot_files <- list(
    boxplot = plot_boxplot(gene_data, gene, file.path(out_dir, "boxplot.png")),
    barplot = plot_barplot(results$stats, gene, file.path(out_dir, "barplot.png")),
    lineplot = plot_lineplot(gene_data, gene, file.path(out_dir, "lineplot.png")),
    
    violin_log2 = plot_violin_log2(gene_data, gene, file.path(out_dir, "violin_log2.png")),
    violin_raw = plot_violin_raw(gene_data, gene, file.path(out_dir, "violin_raw.png")),
    violin_zscore = plot_violin_zscore(gene_data, gene, file.path(out_dir, "violin_zscore.png")),
    violin_median = plot_violin_median(gene_data, gene, file.path(out_dir, "violin_median.png")),
    violin_fixed = plot_violin_fixed(gene_data, gene, file.path(out_dir, "violin_fixed.png")),
    
    boxplot_sex = plot_box_sex(gene_data, gene, file.path(out_dir, "boxplot_sex.png")),
    heatmap_sample = plot_heatmap_sample(gene_data, gene, file.path(out_dir, "heatmap_sample.png")),
    heatmap_stage = plot_heatmap_stage_mean(gene_data, gene, file.path(out_dir, "heatmap_stage_mean.png"))
  )
    
    logger$info(paste("Gráficos gerados:", length(plot_files)))
    
    # FASE 4: ANOTAÇÃO FUNCIONAL
    logger$info("FASE 4: Anotação funcional")
    annotation_script <- file.path(scripts_dir, "prem_annot.R")
    if (file.exists(annotation_script)) {
      system2("Rscript", args = c(annotation_script, shQuote(tpm_file), shQuote(gene)), wait = TRUE)
      logger$info("Anotação concluída")
    } else {
      logger$warn("Script de anotação não encontrado")
    }
    
    # FASE 5: ANÁLISE ESTATÍSTICA AVANÇADA
    logger$info("FASE 5: Análise estatística avançada")
    stats_script <- file.path(scripts_dir, "prem_stats.R")
    tpm_data_file <- file.path(out_dir, "TPM_data_long.csv")
    if (file.exists(tpm_data_file) && file.exists(stats_script)) {
      system2("Rscript", args = c(stats_script, shQuote(tpm_data_file), shQuote(out_dir)), wait = TRUE)
      logger$info("Análise estatística concluída")
    } else {
      logger$warn("Dados ou script de estatística não encontrados")
    }
    
    # FASE 6: RELATÓRIO FINAL COMPLETO
    logger$info("FASE 6: Gerando relatório final completo")
    generate_final_report(gene, out_dir, config, results)
    
    # FASE 7: ATUALIZAR MONITORAMENTO
    logger$info("FASE 7: Atualizando sistema de monitoramento")
    if (exists("create_analysis_dashboard")) {
      tryCatch({
        create_analysis_dashboard(main_results_dir)
      }, error = function(e) {
        logger$warn("Erro no dashboard: ", e$message)
      })
    }
    
    logger$info("=== ANÁLISE COMPLETA CONCLUÍDA COM SUCESSO ===")
    return(list(success = TRUE, output_dir = out_dir, results = results))
    
  }, error = function(e) {
    logger$error("FALHA NO PIPELINE: ", e$message)
    return(list(success = FALSE, error = e$message))
  }, finally = {
    logger$close()
  })
}

# =====================
# Execução Principal
# =====================
main <- function() {
  cat("=== EXECUTANDO PIPELINE COMPLETO ===\n")
  cat("Gene:", gene, "\n")
  cat("Arquivo TPM:", tpm_file, "\n")
  cat("Diretório do projeto:", project_dir, "\n")
  cat("Diretório de resultados:", results_dir, "\n")
  
  # Atualizar configuração com caminhos absolutos
  ANALYSIS_CONFIG$input_dir <- file.path(project_dir, "inputs")
  ANALYSIS_CONFIG$results_dir <- file.path(project_dir, "090-search-gene", "results")
  ANALYSIS_CONFIG$log_dir <- file.path(project_dir, "090-search-gene", "logs")
  ANALYSIS_CONFIG$tpm_file <- tpm_file
  
  # Validar configuração
  validate_config(ANALYSIS_CONFIG)
  
  # Executar análise COMPLETA
  result <- run_complete_gene_analysis(tpm_file, gene, ANALYSIS_CONFIG)
  
  # Status final
  if (result$success) {
    cat("✓ PIPELINE COMPLETO EXECUTADO COM SUCESSO PARA:", gene, "\n")
    cat("✓ RESULTADOS EM:", result$output_dir, "\n")
    cat("✓ ARQUIVOS GERADOS:", length(list.files(result$output_dir, recursive = TRUE)), "\n")
  } else {
    cat("✗ PIPELINE FALHOU PARA:", gene, "\n")
    cat("✗ ERRO:", result$error, "\n")
    quit(status = 1)
  }
}

# Executar
main()
