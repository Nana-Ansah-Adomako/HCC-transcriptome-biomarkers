#To find common genes across all highly DEG files
file1 <- read.csv("results/significantly_expressed_geneset/sig_genes1.csv")
file1 <- file1$x
file2 <- conversion$hgnc_symbol
file3 <- rownames(significant_degs)


gene_list <- intersect(intersect(file1,file2), file3)

write.csv(gene_list, "gene_list.csv")
