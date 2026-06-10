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
gse3 <- getGEO("GSE14520", GSEMatrix = TRUE) # This does not download dataset locally 
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

#-------------------- Step 3: Prepare Phenotype Data ------------#
phenoData <- pData(expression_set)

phenoData <- phenoData %>%
  mutate(disease = case_when(
    !is.na(`Tissue:ch1`) & is.na(`tissue:ch1`) ~ `Tissue:ch1`,
    is.na(`Tissue:ch1`) & !is.na(`tissue:ch1`) ~ `tissue:ch1`,
    !is.na(`Tissue:ch1`) & !is.na(`tissue:ch1`) & `Tissue:ch1` == `tissue:ch1` ~ `Tissue:ch1`,
    !is.na(`Tissue:ch1`) & !is.na(`tissue:ch1`) ~ paste(`Tissue:ch1`, `tissue:ch1`, sep = "_"),
    TRUE ~ NA_character_
  ))

group <- factor(phenoData$disease)
levels(group) <- make.names(levels(group))

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

#-------------------- Step 4: PCA for Batch Check ---------------#
pca_res <- prcomp(t(expr_matrix), scale. = TRUE)
pca_df <- data.frame(pca_res$x, phenoData)

ggplot(pca_df, aes(PC1, PC2, color = phenoData$disease)) +
  geom_point(size = 3) +
  ggtitle("PCA Before Modeling: Check for Batch Effects") +
  theme_minimal()

#-------------------- Step 5: Fit Linear Model ------------------#
contrast.matrix <- makeContrasts(
  TumorVsNormal = Liver.Tumor.Tissue - Liver.Non.Tumor.Tissue,
  levels = design
)

fit <- lmFit(expr_matrix, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

#-------------------- Step 6: Extract DEGs ----------------------#
degs <- topTable(fit2, adjust = "fdr", number = Inf)
significant_degs <- degs[degs$adj.P.Val < 0.05 & abs(degs$logFC) > 1, ]

#-------------------- Step 7: Volcano Plot ----------------------#
EnhancedVolcano(
  degs,
  lab = rownames(degs),
  x = "logFC",
  y = "adj.P.Val",
  title = "Volcano Plot - Differentially Expressed Genes (GSE14520)",
  subtitle = "Significant genes (adj.p < 0.05 and |log2FC| > 1)",
  caption = "Data source: DEG analysis using limma",
  xlab = bquote(~Log[2]~ "fold change"),
  ylab = bquote(~-Log[10]~ "adjusted p-value"),
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 2.0,
  labSize = 3.0
)

#===============volcano plot
# Prepare data frame from DESeq2 results
vol_data <- as.data.frame(degs)
vol_data$gene <- rownames(vol_data)

# Define significance thresholds
pval_cutoff <- 0.05
lfc_cutoff  <- 1  # log2 fold change threshold

# Classify genes
vol_data$group <- "Not Significant"
vol_data$group[vol_data$logFC >  lfc_cutoff & vol_data$adj.P.Val < pval_cutoff] <- "Upregulated"
vol_data$group[vol_data$logFC < -lfc_cutoff & vol_data$adj.P.Val < pval_cutoff] <- "Downregulated"
vol_data$group <- factor(vol_data$group, levels = c("Upregulated", "Downregulated", "Not Significant"))

# Remove NA rows
vol_data <- vol_data[!is.na(vol_data$adj.P.Val) & !is.na(vol_data$logFC), ]

# Plot
ggplot(vol_data, aes(x = logFC, y = -log10(adj.P.Val), color = group)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Upregulated"     = "darkred",
    "Downregulated"   = "steelblue",
    "Not Significant" = "gray60"
  )) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(pval_cutoff),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  labs(
    title    = "Differentially Expressed Genes of GSE14520 Dataset",
    x        = expression(log[2]~"Fold Change"),
    y        = expression(-log[10]~"(Adjusted p-value)"),
    color    = "Expression"
  ) +
  theme_classic() +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position = "top"
  )
#-------------------- Step 8: Save Results ----------------------#
write.csv(significant_degs, "results/significantly_expressed_geneset/sig_genes3.csv")
