# 01_preprocess_GSE146719.R
# Preprocessing and DESeq2 differential expression analysis for GSE146719

# Load libraries
library(GEOquery)
library(DESeq2)
library(dplyr)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(limma)
library(biomaRt)
library(tibble)
library(ggplot2)

#----------------- Dataset 1: GSE146719 ------------------#
gse1 <- getGEO("GSE146719", GSEMatrix = TRUE)[[1]]

count_data1 <- read.table("/home/nana-ansah-adomako/HCC/GSE146719_All.counts.RNA.csv")
pheno_data1 <- pData(gse1)

# Check if column names of count matrix match phenotype rownames
all(colnames(count_data1) %in% rownames(pheno_data1))

# Select first six rows of phenotype data
pheno_data1 <- pheno_data1[1:6,]

# Set rownames to sample names from phenotype data
rownames(pheno_data1) <- pheno_data1$source_name_ch1

# Ensure sample ordering matches
all(colnames(count_data1) == rownames(pheno_data1))

# Make first column gene names and first row sample names
rownames(count_data1) <- count_data1$V1
colnames(count_data1) <- count_data1[1, ]

# Remove first row and first column
count_data1 <- count_data1[-1, ]
count_data1 <- count_data1[, -1]

# Check order again
all(colnames(count_data1) == rownames(pheno_data1))

# Reorder columns to match phenotype order
colorder <- rownames(pheno_data1)
count_data1 <- count_data1[, colorder]

# Store rownames and colnames
rownames <- rownames(count_data1)
colnames <- colnames(count_data1)

# Convert counts to integers
count_data1 <- as.data.frame(sapply(count_data1, as.integer))
row.names(count_data1) <- rownames

# Define condition variable manually
dds1_condition <- as.factor(pheno_data1$`tissue type:ch1`)
pheno_data1$condition <- as.factor(c("tumor", "non_tumor"))

# Create DESeq2 dataset
dds1 <- DESeqDataSetFromMatrix(countData = count_data1,
                               colData = pheno_data1,
                               design = ~ condition)

#-------------- Batch Effect Check (PCA) --------------#
vsd <- vst(dds1, blind = TRUE)
plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA Plot of GSE146719 Samples by Condition")

#--------------------- Differential Expression ---------------------#

# Set reference level (non_tumor)
dds1$condition <- relevel(dds1$condition, ref = "non_tumor")

# Prefilter genes with fewer than 10 reads total
dds1 <- dds1[rowSums(counts(dds1)) >= 10,]

# Run DESeq2
dds1 <- DESeq(dds1)
res1 <- results(dds1)

# Extract significant genes
sig_genes1_all <- res1[!is.na(res1$padj) & res1$padj < 0.05, ]
sig_genes1 <- rownames(res1)[which(res1$padj < 0.05 & !is.na(res1$padj))]

# MA plot
plotMA(res1,
       main = "MA Plot - Differentially Expressed Genes of GSE146719",
       colNonSig = "gray60", colSig = "darkred",
       colLine = "grey40")

# Save results
write.csv(sig_genes1_all, "sig_genes1_all.csv")
write.csv(sig_genes1, "sig_genes1.csv")
