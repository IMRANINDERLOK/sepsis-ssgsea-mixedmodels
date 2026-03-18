# Check working directory
getwd()
# Install BiocManager (only once)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages (only once)
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE)

# Install CRAN packages (only once)
install.packages(c("tidyverse", "pheatmap"))

# Load libraries (every session)
library(GEOquery)
library(limma)
library(tidyverse)
library(pheatmap)
# Download GSE54514 dataset
gset <- getGEO("GSE54514", GSEMatrix = TRUE)

# Check how many platforms
length(gset)

# Extract first platform
gse <- gset[[1]]

# Extract expression matrix
expr <- exprs(gse)

# Extract metadata
pheno <- pData(gse)

# Check dimensions
dim(expr)
dim(pheno)

# Confirm matching between expression and metadata
all(colnames(expr) == rownames(pheno))
# ================================
# INSTALL REQUIRED PACKAGES (RUN ONCE)
# ================================

# Install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages
BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
BiocManager::install("limma", ask = FALSE, update = FALSE)

# Install CRAN packages
install.packages("tidyverse")
install.packages("pheatmap")
# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install GEOquery and limma
BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
BiocManager::install("limma", ask = FALSE, update = FALSE)
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE)
install.packages(c("tidyverse", "pheatmap"))
library(GEOquery)
library(limma)
library(tidyverse)
library(pheatmap)
gset <- getGEO("GSE54514", GSEMatrix = TRUE)
length(gset)

gse <- gset[[1]]
expr <- exprs(gse)
pheno <- pData(gse)

dim(expr)
dim(pheno)

all(colnames(expr) == rownames(pheno))
# ================================
# LOAD NON-NORMALIZED GEO FILE
# ================================

# Read expression file
expr_raw <- read.delim(
  "data/GSE54514_non-normalized.txt.gz",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Check dimensions
dim(expr_raw)

# See first few rows
head(expr_raw[, 1:5])
list.files("data")
# ================================
# READ NON-NORMALIZED EXPRESSION FILE
# ================================

expr_raw <- read.delim(
  "data/GSE54514_non-normalized.txt.gz",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Check dimensions
dim(expr_raw)
expr_raw <- read.delim(
  "data/GSE54514_non-normalized.txt.gz",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(expr_raw)
## =========================
## CHECK SCRIPT: expr_matrix
## =========================

cat("=== BASIC DIMENSIONS ===\n")
print(dim(expr_matrix))

cat("\n=== COLUMN NAME PREVIEW ===\n")
cat("First 10 column names:\n")
print(head(colnames(expr_matrix), 10))
cat("Last 10 column names:\n")
print(tail(colnames(expr_matrix), 10))

cat("\n=== DUPLICATE COLUMN NAMES ===\n")
dup_n <- sum(duplicated(colnames(expr_matrix)))
cat("Number of duplicated column names:", dup_n, "\n")
if (dup_n > 0) {
  cat("Example duplicated names:\n")
  print(unique(colnames(expr_matrix)[duplicated(colnames(expr_matrix))])[1:10])
}

cat("\n=== ARE COLUMN NAMES GSM IDs? ===\n")
gsm_like <- grepl("^GSM[0-9]+$", colnames(expr_matrix))
cat("Columns that look like GSM IDs:", sum(gsm_like), "out of", ncol(expr_matrix), "\n")
cat("Example GSM-like names:\n")
print(head(colnames(expr_matrix)[gsm_like], 10))

cat("\n=== CHECK FOR NON-NUMERIC COLUMNS (COMMON IMPORT ISSUE) ===\n")
non_numeric_cols <- which(!sapply(as.data.frame(expr_matrix), is.numeric))
cat("Non-numeric columns found:", length(non_numeric_cols), "\n")
if (length(non_numeric_cols) > 0) {
  cat("Indices of first non-numeric columns:\n")
  print(head(non_numeric_cols, 10))
  cat("Names of first non-numeric columns:\n")
  print(head(colnames(expr_matrix)[non_numeric_cols], 10))
}

cat("\n=== CHECK NA RATE (IMPORT QUALITY) ===\n")
na_rate <- mean(is.na(expr_matrix))
cat("Overall NA rate in expr_matrix:", round(na_rate, 6), "\n")

cat("\n=== CHECK RANGE (LOG2 OR NOT) ===\n")
rng <- range(expr_matrix, na.rm = TRUE)
cat("Range:", rng[1], "to", rng[2], "\n")
if (rng[2] > 100) {
  cat("NOTE: values look non-log scale (max > 100). We'll likely need log2(x+1) later.\n")
} else {
  cat("NOTE: values may already be log-scale.\n")
}

cat("\n=== QUICK CHECK: ARE THERE EXACTLY DOUBLE COLUMNS? (326 vs 163) ===\n")
cat("If GEO shows ~163 samples but you have", ncol(expr_matrix),
    "columns, it may be duplicated columns or you imported two matrices.\n")
cat("Try checking if 326/163 is close to 2:\n")
cat("326/163 =", ncol(expr_matrix)/163, "\n")
# Keep only expression intensity columns (drop Detection P-values)
is_detp <- grepl("\\.Detection Pval$", colnames(expr_matrix))

expr_int <- expr_matrix[, !is_detp]   # 48803 x 163
detP     <- expr_matrix[,  is_detp]   # 48803 x 163  (optional, keep for QC)

dim(expr_int)
dim(detP)
