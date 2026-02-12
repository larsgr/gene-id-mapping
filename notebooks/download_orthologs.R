#!/usr/bin/env Rscript
# Download human ortholog data from Ensembl BioMart and NCBI.
#
# Outputs:
#   notebooks/ens_human_orthologs.tsv   - Ensembl BioMart orthologs
#   notebooks/ncbi_human_orthologs.tsv  - NCBI gene_orthologs (salmon-human)
#
# Usage:
#   Rscript notebooks/download_orthologs.R

library(biomaRt)
library(readr)

# Determine script directory from command-line args or fallback to getwd()
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) {
  script_dir <- dirname(normalizePath(script_path))
} else {
  script_dir <- getwd()
}

# ---------------------------------------------------------------------------
# 1. Ensembl BioMart: salmon-human orthologs
# ---------------------------------------------------------------------------
cat("Downloading Ensembl orthologs via BioMart...\n")

ensembl <- useEnsembl(biomart = "genes", dataset = "ssalar_gene_ensembl")

ens_orthologs <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "hsapiens_homolog_ensembl_gene",
    "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_orthology_type",
    "hsapiens_homolog_orthology_confidence"
  ),
  mart = ensembl
)

# Filter to rows that actually have a human ortholog
ens_orthologs <- ens_orthologs[ens_orthologs$hsapiens_homolog_ensembl_gene != "", ]

ens_out <- file.path(script_dir, "ens_human_orthologs.tsv")
write_tsv(ens_orthologs, ens_out)
cat(sprintf("  Wrote %d ortholog pairs to %s\n", nrow(ens_orthologs), ens_out))

# ---------------------------------------------------------------------------
# 2. NCBI gene_orthologs: salmon-human orthologs
# ---------------------------------------------------------------------------
cat("Downloading NCBI gene_orthologs...\n")

# Download gene_orthologs.gz
ortho_url <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz"
ortho_tmp <- tempfile(fileext = ".gz")
download.file(ortho_url, ortho_tmp, mode = "wb", quiet = TRUE)

# Read and filter for salmon (8030) <-> human (9606)
# Header line starts with #tax_id, so read it and fix the column name
ortho_all <- read_tsv(ortho_tmp, show_col_types = FALSE)
if ("#tax_id" %in% names(ortho_all)) {
  names(ortho_all)[names(ortho_all) == "#tax_id"] <- "tax_id"
}
# File stores pairs directionally: human (9606) in tax_id, salmon (8030) in Other_tax_id
salmon_human <- ortho_all[
  ortho_all$tax_id == 9606 & ortho_all$Other_tax_id == 8030,
]
cat(sprintf("  Found %d salmon-human ortholog pairs\n", nrow(salmon_human)))

# Download gene_info for human gene symbols
cat("Downloading NCBI human gene_info for symbols...\n")
gene_info_url <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
gene_info_tmp <- tempfile(fileext = ".gz")
download.file(gene_info_url, gene_info_tmp, mode = "wb", quiet = TRUE)

human_genes <- read_tsv(gene_info_tmp, show_col_types = FALSE)
if ("#tax_id" %in% names(human_genes)) {
  names(human_genes)[names(human_genes) == "#tax_id"] <- "tax_id"
}

# Build lookup: human GeneID -> Symbol
human_lookup <- human_genes[, c("GeneID", "Symbol")]
names(human_lookup) <- c("GeneID", "human_symbol")

# Join human symbols via GeneID (which is the human gene in this direction)
ncbi_orthologs <- merge(salmon_human, human_lookup, by = "GeneID", all.x = TRUE)

# Select and rename: Other_GeneID is the salmon gene, GeneID is the human gene
ncbi_orthologs <- ncbi_orthologs[, c(
  "Other_GeneID", "Other_tax_id", "relationship",
  "GeneID", "tax_id", "human_symbol"
)]
names(ncbi_orthologs) <- c(
  "salmon_gene_id", "salmon_tax_id", "relationship",
  "human_gene_id", "human_tax_id", "human_symbol"
)

ncbi_out <- file.path(script_dir, "ncbi_human_orthologs.tsv")
write_tsv(ncbi_orthologs, ncbi_out)
cat(sprintf("  Wrote %d ortholog pairs to %s\n", nrow(ncbi_orthologs), ncbi_out))

# Clean up temp files
unlink(c(ortho_tmp, gene_info_tmp))

cat("Done.\n")
