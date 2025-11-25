#============================================================#
# Script: 03_process_GSE14520.R
# Purpose: Differential Expression Analysis for GSE14520 (Microarray)
# Author: Borbi Stephen
#============================================================#

library(GEOquery)
library(limma)
library(hgu133plus2.db)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(tibble)

#-------------------- Step 1: Load GEO Data --------------------#
gse3 <- getGEO("GSE14520", GSEMatrix = TRUE)
expression_set <- gse3[[1]]

# Extract expression matrix
expr_matrix <- exprs(expression_set)
probes <- rownames(expr_matrix)

#-------------------- Step 2: Map Probes to Gene Symbols --------#
gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = probes,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

# Add gene symbols to expression matrix
expr_matrix <- cbind(GeneSymbol = gene_symbols, expr_matrix)
expr_matrix <- expr_matrix[!is.na(expr_matrix[, "GeneSymbol"]), ]

 # Handle duplicated symbols by keeping row with highest variance
expr_matrix <- as.data.frame(expr_matrix)
expr_matrix <- expr_matrix %>%
  group_by(GeneSymbol) %>%
  mutate(Variance = apply(across(where(is.numeric)), 1, var, na.rm = TRUE)) %>%
  slice_max(order_by = Variance, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-Variance) %>%
  column_to_rownames("GeneSymbol")

# Ensure all data are numeric
rownames_expr <- rownames(expr_matrix)
expr_matrix <- as.data.frame(sapply(expr_matrix, as.numeric))
rownames(expr_matrix) <- rownames_expr
