# reporting_utils.R
generate_final_report <- function(gene, out_dir, config, results) {
  # Summary JSON
  summary <- list(
    gene = gene,
    analysis_date = Sys.time(),
    status = "COMPLETED",
    output_files = list.files(out_dir, recursive = TRUE),
    summary_stats = list(
      n_samples = if (!is.null(results$qc$stage_coverage)) nrow(results$qc$stage_coverage) else 0,
      mean_expression = if (!is.null(results$stats)) mean(results$stats$mean_tpm, na.rm = TRUE) else NA,
      expression_range = if (!is.null(results$stats)) range(results$stats$mean_tpm, na.rm = TRUE) else c(NA, NA)
    ),
    qc_checks = if (!is.null(results$qc$checks)) results$qc$checks else NULL
  )
  
  writeLines(
    jsonlite::toJSON(summary, pretty = TRUE),
    file.path(out_dir, "analysis_summary.json")
  )
  
  # Relatório detalhado em texto
  sink(file.path(out_dir, "analysis_report.txt"))
  cat("RELATÓRIO COMPLETO DE ANÁLISE -", gene, "\n")
  cat("==========================================\n\n")
  cat("Data da análise:", as.character(Sys.time()), "\n")
  cat("Diretório de saída:", out_dir, "\n\n")
  
  if (!is.null(results$qc)) {
    cat("CONTROLE DE QUALIDADE:\n")
    cat("----------------------\n")
    for (check_name in names(results$qc$checks)) {
      check <- results$qc$checks[[check_name]]
      cat(sprintf("- %s: %s (%s)\n", check_name, check$status, check$message))
    }
    cat("\n")
    
    if (!is.null(results$qc$stage_coverage)) {
      cat("COBERTURA POR ESTÁGIO:\n")
      print(results$qc$stage_coverage)
      cat("\n")
    }
  }
  
  if (!is.null(results$stats)) {
    cat("ESTATÍSTICAS DESCRITIVAS:\n")
    cat("-------------------------\n")
    print(results$stats)
    cat("\n")
  }
  
  cat("ARQUIVOS GERADOS:\n")
  cat("-----------------\n")
  files <- list.files(out_dir, recursive = TRUE)
  for (file in files) {
    cat("-", file, "\n")
  }
  
  sink()
  
  cat("✓ Relatório gerado:", file.path(out_dir, "analysis_report.txt"), "\n")
}
