#--------------------  Load GEO Data -------------------#
gse2 <- getGEO("GSE142987", GSEMatrix = TRUE)[[1]]
pheno_data2 <- pData(gse2)

#-------------------- Get Supplementary Count Matrix ---#
# GEOquery does not always return raw counts directly
gse_2 <- getGEOSuppFiles("GSE142987")

# List the files and locate the count matrix
list.files("GSE142987")

# Load count matrix (adjust the path to your working directory)
count_data2 <- read.table("data/raw/GSE142987_sample_count_matrix.txt")

#-------------------- Data Formatting -------------------#
# Make the first column gene IDs and row names
rownames(count_data2) <- count_data2[, 1]
count_data2 <- count_data2[, -1]

# Assign sample names from phenotype data
row_name_phen2 <- rownames(pheno_data2)
colnames(count_data2) <- row_name_phen2

# Define condition factor (tumor vs. non_tumor)
pheno_data2$condition <- as.factor(pheno_data2$`disease state:ch1`)
pheno_data2$condition <- factor(pheno_data2$condition, 
                                levels = c("liver cancer patient", "healthy donor"), 
                                labels = c("tumor", "non_tumor"))

# Convert count matrix to numeric
rownames_dat2 <- rownames(count_data2)
count_data2 <- as.data.frame(sapply(count_data2, as.integer))
rownames(count_data2) <- rownames_dat2

# Remove any NA rows (often first row)
count_data2 <- count_data2[-1, ]

#-------------------- Create DESeq2 Object --------------#
dds2 <- DESeqDataSetFromMatrix(countData = count_data2,
                               colData = pheno_data2,
                               design = ~ condition)

#-------------------- Batch Effect Check ----------------#
vsd2 <- vst(dds2, blind = TRUE)
plotPCA(vsd2, intgroup = "condition") +
  ggtitle("PCA Plot of GSE142987 Samples by Condition")

#-------------------- Differential Expression ------------#
dds2$condition <- relevel(dds2$condition, ref = "non_tumor")

# Pre-filter: keep genes with at least 10 total reads
dds2 <- dds2[rowSums(counts(dds2)) >= 10, ]

# Run DESeq pipeline
dds2 <- DESeq(dds2)
res2 <- results(dds2)

#-------------------- Extract Significant Genes ----------#
sig_genes2_all <- res2[!is.na(res2$padj) & res2$padj < 0.05, ]
sig_genes2 <- rownames(res2)[which(res2$padj < 0.05 & !is.na(res2$padj))]

#-------------------- Visualization ----------------------#
plotMA(res2,
       main = "MA Plot - Differentially Expressed Genes (GSE142987)",
       colNonSig = "gray60",
       colSig = "darkred",
       colLine = "black")

#===============volcano plot
# Prepare data frame from DESeq2 results
vol_data <- as.data.frame(res2)
vol_data$gene <- rownames(vol_data)

# Define significance thresholds
pval_cutoff <- 0.05
lfc_cutoff  <- 1  # log2 fold change threshold

# Classify genes
vol_data$group <- "Not Significant"
vol_data$group[vol_data$log2FoldChange >  lfc_cutoff & vol_data$padj < pval_cutoff] <- "Upregulated"
vol_data$group[vol_data$log2FoldChange < -lfc_cutoff & vol_data$padj < pval_cutoff] <- "Downregulated"
vol_data$group <- factor(vol_data$group, levels = c("Upregulated", "Downregulated", "Not Significant"))

# Remove NA rows
vol_data <- vol_data[!is.na(vol_data$padj) & !is.na(vol_data$log2FoldChange), ]

# Plot
ggplot(vol_data, aes(x = log2FoldChange, y = -log10(padj), color = group)) +
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
    title    = "Differentially Expressed Genes of GSE142987 Dataset",
    x        = expression(log[2]~"Fold Change"),
    y        = expression(-log[10]~"(Adjusted p-value)"),
    color    = "Expression"
  ) +
  theme_classic() +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position = "top"
  )



#-------------------- Save Results -----------------------#
write.csv(sig_genes2, "results/significantly_expressed_geneset/sig_genes2.csv")
write.csv(sig_genes2_all, "results/significantly_expressed_geneset/sig_genes2_all.csv")

#-------------------- Annotate Gene Symbols -------------#
sig_genes2_all_name <- read.csv("results/significantly_expressed_geneset/sig_genes2_all.csv")

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
sig_genes2_all_name$Gene_Symbol <- gsub("\\.\\d+", "", sig_genes2_all_name$X)
ensembl_ids <- sig_genes2_all_name$Gene_Symbol

conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values = ensembl_ids,
                    mart = ensembl)

colnames(sig_genes2_all_name)[colnames(sig_genes2_all_name) == "Gene_Symbol"] <- "ensembl_gene_id"
sig_genes2_all_name_annotated <- left_join(sig_genes2_all_name, conversion, by = "ensembl_gene_id")

write.csv(sig_genes2_all_name_annotated, "results/significantly_expressed_geneset/sig_genes2_all_name_annotated.csv")
