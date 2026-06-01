#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(seqinr)
  library(rtracklayer)
  library(GO.db)
  library(Peptides)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Uso: Rscript prem_annot.R <TPM_file> <gene>")

tpm_file <- args[1]
gene <- args[2]
out_dir <- file.path("../results", gene)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================
# Configurações
# ============================
gff_file <- "../inputs/schistosoma_mansoni.PRJEA36577.WBPS19.annotations.gff3"
protein_fasta <- "../inputs/schistosoma_mansoni.PRJEA36577.WBPS19.protein.fa"
cds_fasta <- "../inputs/schistosoma_mansoni.PRJEA36577.WBPS19.CDS_transcripts.fa"

# ============================
# 1. Ler GFF
# ============================
gff <- import.gff3(gff_file)
gff_df <- as.data.frame(gff)

# Corrigir nome: buscar por diferentes formatos
gene_patterns <- c(
  gene,
  paste0("gene:", gene),
  gsub("^gene:", "", gene),
  gsub("_", "-", gene),
  gsub(":", "_", gene)
)

gene_entry <- gff_df %>%
  filter(type == "gene") %>%
  filter(ID %in% gene_patterns)

if (nrow(gene_entry) == 0) {
  cat("ERRO: Gene", gene, "não encontrado no GFF\n")
  quit(status = 1)
}

# ============================
# 2. Buscar transcritos e isoformas
# ============================
transcripts <- gff_df %>%
  filter(type == "mRNA") %>%
  filter(str_detect(Parent, gene_entry$ID[1]))

isoform_info <- list()
if (nrow(transcripts) > 0) {
  for (tid in transcripts$ID) {
    tx <- transcripts %>% filter(ID == tid)
    tx_features <- gff_df %>% filter(Parent == tid)

    tx_exons <- tx_features %>% filter(type == "exon")
    tx_cds   <- tx_features %>% filter(type == "CDS")

    isoform_info[[tid]] <- tibble(
      transcript_id = tid,
      location = paste0(tx$seqnames[1], ":", tx$start[1], "-", tx$end[1]),
      strand = as.character(tx$strand[1]),
      n_exons = nrow(tx_exons),
      cds_length = sum(tx_cds$width, na.rm = TRUE)
    )
  }
  isoform_info <- bind_rows(isoform_info)
} else {
  isoform_info <- tibble(
    transcript_id = NA,
    location = NA,
    strand = NA,
    n_exons = NA,
    cds_length = NA
  )
}

# ============================
# 3. GO Terms
# ============================
if ("Ontology_term" %in% colnames(gene_entry)) {
  ontology <- gene_entry$Ontology_term[1]
} else {
  ontology <- NA
}

go_desc <- "Sem GO terms disponíveis"
if (!is.na(ontology) && length(ontology) > 0) {
  go_terms <- str_split(ontology, ",")[[1]] %>% str_trim()
  go_desc <- map_chr(go_terms, ~ {
    term <- GOTERM[[.x]]
    if (is.null(term)) return(paste(.x, "- descrição não encontrada"))
    paste0(.x, " - ", Term(term), " [", Ontology(term), "]")
  })
}


# ============================
# 4. Sequências
# ============================
seqs_prot <- read.fasta(protein_fasta, seqtype = "AA", as.string = TRUE)
prot_id <- grep(gsub("^gene:", "", gene_entry$ID[1]), names(seqs_prot), value = TRUE)[1]
prot_seq <- if (!is.na(prot_id)) seqs_prot[[prot_id]][1] else NA

seqs_cds <- read.fasta(cds_fasta, seqtype = "DNA", as.string = TRUE)
cds_id <- grep(gsub("^gene:", "", gene_entry$ID[1]), names(seqs_cds), value = TRUE)[1]
cds_seq <- if (!is.na(cds_id)) seqs_cds[[cds_id]][1] else NA

# ============================
# 5. Propriedades da proteína
# ============================
prot_summary <- NULL
if (!is.na(prot_seq)) {
  prot_summary <- tibble(
    length = nchar(prot_seq),
    mw = mw(prot_seq),
    pI = pI(prot_seq),
    hydrophobicity = hydrophobicity(prot_seq),
    charge = charge(prot_seq, pH = 7)
  )
}

# ============================
# 6. Exportar relatórios
# ============================
write_csv(isoform_info, file.path(out_dir, "isoforms.csv"))
if (!is.null(prot_summary)) write_csv(prot_summary, file.path(out_dir, "protein_properties.csv"))
if (!is.na(prot_seq)) write_lines(prot_seq, file.path(out_dir, "protein.fasta"))
if (!is.na(cds_seq)) write_lines(cds_seq, file.path(out_dir, "cds.fasta"))

sink(file.path(out_dir, "annotation_report.txt"))
cat("Relatório de Anotação -", gene, "\n")
cat("==============================\n\n")

cat("Descrição: ", gene_entry$description[1], "\n\n")

cat("Isoformas encontradas (transcritos):\n")
print(isoform_info)
cat("\n")

cat("GO Terms:\n")
print(go_desc)
cat("\n")

if (!is.null(prot_summary)) {
  cat("Propriedades da proteína principal:\n")
  print(prot_summary)
}
sink()

cat("prem_annot.R concluído para", gene, "\n")
