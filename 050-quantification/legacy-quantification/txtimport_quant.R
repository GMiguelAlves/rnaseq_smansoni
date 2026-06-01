#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(readr)
  library(dplyr)
  library(rtracklayer)
  library(tibble)
})

project_dir <- "../"
meta_file   <- file.path(project_dir, "010-reference", "RNAseq_metadata.tsv")
gtf_file    <- file.path(project_dir, "010-reference", "data", "schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3.gtf")
quant_dir   <- file.path(project_dir, "040-alignment", "quants")
out_dir     <- file.path(project_dir, "050-quantification")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- read.delim(meta_file, header = TRUE, sep = "\t")
samples <- meta[[3]]

cat("[INFO] Extraindo relação transcript-gene do GTF...\n")

gtf_data <- rtracklayer::import(gtf_file)
tx2gene <- as.data.frame(gtf_data) %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, gene_id) %>%
  dplyr::mutate(
    
    transcript_id = gsub("^transcript:", "", transcript_id),
    transcript_id = gsub("\\.[0-9]+$", "", transcript_id),
    gene_id = ifelse(grepl("^transcript:", gene_id),
                     gsub("^transcript:", "gene:", gene_id),
                     gene_id),
    
    gene_id = gsub("^gene:", "", gene_id),
    gene_id = gsub("\\.[0-9]+$", "", gene_id)
  ) %>%
  distinct()

files <- file.path(quant_dir, samples, "quant.sf")
names(files) <- samples

missing <- files[!file.exists(files)]
if (length(missing) > 0) {
  stop("[ERRO] Arquivos quant.sf não encontrados para: ", paste(names(missing), collapse = ", "))
}

cat("[INFO] Rodando tximport...\n")
txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                countsFromAbundance = "no",
                ignoreTxVersion = TRUE,
                ignoreAfterBar = TRUE)

cat("[INFO] Salvando matrizes...\n")

# Contagens brutas
write_tsv(as.data.frame(txi$counts) %>% rownames_to_column("gene_id"),
          file.path(out_dir, "counts_matrix.tsv"))

# TPM
write_tsv(as.data.frame(txi$abundance) %>% rownames_to_column("gene_id"),
          file.path(out_dir, "tpm_matrix.tsv"))

cat("[OK] Matrizes salvas em 050-quantification/\n")
