########################################
# setup
########################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

libs <- c("GEOquery", "DESeq2", "limma", "org.Hs.eg.db", "biomaRt",
          "EnhancedVolcano", "ggplot2", "dplyr", "tibble", "hgu133plus2.db")

BiocManager::install(libs, ask = FALSE, update = TRUE)
lapply(libs, library, character.only = TRUE)

