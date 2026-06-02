# qc_utils.R
perform_data_qc <- function(tpm_data, gene, config) {
  qc_report <- list(
    gene = gene,
    timestamp = Sys.time(),
    checks = list()
  )
  
  # Check 1: Presença do gene
  if (nrow(tpm_data) == 0) {
    qc_report$checks$gene_presence <- list(
      status = "FAIL",
      message = "Gene não encontrado nos dados TPM"
    )
    return(qc_report)
  }
  
  qc_report$checks$gene_presence <- list(
    status = "PASS", 
    message = paste("Gene encontrado com", nrow(tpm_data), "amostras")
  )
  
  # Check 2: Valores TPM válidos
  tpm_values <- tpm_data$TPM
  valid_tpm <- sum(is.finite(tpm_values) & tpm_values >= 0)
  valid_percent <- valid_tpm / length(tpm_values) * 100
  
  qc_report$checks$tpm_validity <- list(
    status = ifelse(valid_percent >= 95, "PASS", "WARN"),
    message = paste0(round(valid_percent, 1), "% dos valores TPM são válidos")
  )
  
  # Check 3: Cobertura por estágio
  stage_coverage <- tpm_data %>%
    group_by(stage_category) %>%
    summarise(
      n_samples = n(),
      mean_tpm = mean(TPM, na.rm = TRUE),
      detected = sum(TPM > config$min_tpm_threshold, na.rm = TRUE)
    )
  
  qc_report$stage_coverage <- stage_coverage
  
  # Check 4: Outliers
  outlier_threshold <- median(tpm_values, na.rm = TRUE) + 
    config$quality_control$outlier_sd_threshold * sd(tpm_values, na.rm = TRUE)
  
  outliers <- sum(tpm_values > outlier_threshold, na.rm = TRUE)
  qc_report$checks$outliers <- list(
    status = ifelse(outliers == 0, "PASS", "WARN"),
    message = paste(outliers, "possíveis outliers detectados")
  )
  
  return(qc_report)
}

generate_qc_report <- function(qc_report, out_dir) {
  qc_file <- file.path(out_dir, "quality_control_report.json")
  writeLines(jsonlite::toJSON(qc_report, pretty = TRUE), qc_file)
  
  # Relatório resumido para log
  cat("=== RELATÓRIO DE QUALIDADE ===\n")
  for (check_name in names(qc_report$checks)) {
    check <- qc_report$checks[[check_name]]
    cat(sprintf("%s: %s - %s\n", check_name, check$status, check$message))
  }
  cat("==============================\n")
}
