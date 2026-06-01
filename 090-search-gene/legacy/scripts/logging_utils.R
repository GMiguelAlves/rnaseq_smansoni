# logging_utils.R
create_logger <- function(gene, log_dir) {
  log_file <- file.path(log_dir, paste0("prem_analysis_", gene, ".log"))
  con <- file(log_file, open = "a")
  
  list(
    info = function(...) {
      msg <- paste("[INFO]", Sys.time(), "-", ...)
      writeLines(msg, con)
      flush(con)
    },
    warn = function(...) {
      msg <- paste("[WARN]", Sys.time(), "-", ...)
      writeLines(msg, con)
      warning(msg, call. = FALSE)
    },
    error = function(...) {
      msg <- paste("[ERROR]", Sys.time(), "-", ...)
      writeLines(msg, con)
      stop(msg, call. = FALSE)
    },
    close = function() close(con)
  )
}

