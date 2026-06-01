# monitoring_utils.R
create_analysis_dashboard <- function(results_dir) {
  analysis_summaries <- list.files(
    results_dir, 
    pattern = "analysis_summary.json", 
    recursive = TRUE, 
    full.names = TRUE
  )
  
  dashboard_data <- map_df(analysis_summaries, function(summary_file) {
    summary <- jsonlite::fromJSON(summary_file)
    tibble(
      gene = summary$gene,
      status = summary$status,
      date = as.POSIXct(summary$analysis_date),
      n_samples = summary$summary_stats$n_samples,
      mean_expression = summary$summary_stats$mean_expression
    )
  })
  
  # Gerar dashboard
  dashboard_plot <- ggplot(dashboard_data, aes(x = date, y = mean_expression, color = status)) +
    geom_point(alpha = 0.6) +
    theme_bw() +
    labs(title = "Dashboard de Análises Concluídas",
         x = "Data de Análise", y = "Expressão Média (TPM)")
  
  ggsave(file.path(results_dir, "analysis_dashboard.png"), dashboard_plot)
  
  return(dashboard_data)
}
