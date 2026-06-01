# config_analysis.R

# Esta configuração será sobrescrita pelo analysis_pipeline.R
# com os caminhos absolutos corretos

ANALYSIS_CONFIG <- list(
  # Diretórios (serão sobrescritos)
  input_dir = "../inputs",
  results_dir = "../results",
  log_dir = "../logs",
  
  # Arquivos (serão sobrescritos)
  tpm_file = "../../050-quantification/tpm_matrix.tsv",
  gff_file = "../inputs/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3",
  protein_fasta = "../inputs/schistosoma_mansoni.PRJEA36577.WBPS19.protein.fa",
  cds_fasta = "../inputs/schistosoma_mansoni.PRJEA36577.WBPS19.CDS_transcripts.fa",
  
  # Parâmetros de análise
  min_tpm_threshold = 0.1,
  quality_control = list(
    min_samples_per_stage = 2,
    max_missing_rate = 0.3,
    outlier_sd_threshold = 3
  ),
  
  # Visualização
  plot_formats = c("png", "pdf"),
  plot_dpi = 300,
  plot_width = 10,
  plot_height = 7,
  
  # Estatísticas
  p_value_threshold = 0.05,
  fdr_threshold = 0.1
)

# Função para validar configuração
validate_config <- function(config) {
  cat("Validando configuração...\n")
  
  # Verificar diretórios
  dirs_to_check <- c(config$input_dir, config$results_dir, config$log_dir)
  for (dir in dirs_to_check) {
    if (!dir.exists(dir)) {
      cat("AVISO: Diretório não existe:", dir, "\n")
    }
  }
  
  # Verificar arquivos essenciais
  required_files <- c(config$tpm_file, config$gff_file, 
                     config$protein_fasta, config$cds_fasta)
  
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    cat("AVISO: Arquivos de configuração ausentes: ", paste(missing_files, collapse = ", "), "\n")
  }
  
  cat("✓ Configuração validada\n")
  return(TRUE)
}
