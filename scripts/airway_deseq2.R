library(DESeq2)
library(tidyverse)
library(airway)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(apeglm)


#Read sample Info
counts_df <- read.csv("airway_counts.csv", header = TRUE, sep = ",",  row.names = 1)
head(counts_df)
metadata <- read.csv("airway_metadata.csv", header = TRUE, sep = ",",  row.names = 1)
head(metadata)

#Make sure columns match and are in the same order
colnames(counts_df)
rownames(metadata)
colnames(counts_df) <- trimws(colnames(counts_df))
rownames(metadata) <- trimws(rownames(metadata))
all(colnames(counts_df) %in% rownames(metadata))
all(colnames(counts_df) == rownames(metadata))

#Filter out Samples with missing data
samples_to_keep <- !is.na(metadata$dex)
metadata <- metadata[samples_to_keep, ]
counts_df <- counts_df[, samples_to_keep]

#Filter out genes with low expression
genes_to_keep <- rowSums(counts_df) >= 10
counts_df <- counts_df[genes_to_keep, ]

#Construct a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_df,
                              colData = metadata,
                              design = ~ cell + dex)

#Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep]

#set the factor level, set reference level
dds$dex <- relevel(dds$dex, red = "untrt")

#run DESeq
dds <- DESeq(dds)
res <- results(dds)

#explore results
summary(res)

#change summary to summarize adjusted p-value < 0.05
res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)

#contrasts
res_contrast <- results(dds, contrast = c("dex", "trt", "untrt"))
summary(res_contrast)
res_contrast

# MA-Plot
plotMA(res)

#Top 20 genes
res_ordered <- res_contrast[order(res_contrast$padj), ]
res_ordered <- res_ordered[!is.na(res_ordered$padj), ]
top20_genes <- head(res_ordered, 20)
top20_genes

#Volcano Plot
EnhancedVolcano(res_contrast,
                lab = rownames(res_contrast),           
                x = 'log2FoldChange',                    
                y = 'padj',                              
                title = 'Volcano Plot: trt vs untrt',
                pCutoff = 0.05,                          
                )

#Heatmap
# Remove NAs
res_df <- as.data.frame(res_contrast)
res_df <- res_df[!is.na(res_df$padj), ]

# Select top 30 most significant genes
top_genes <- rownames(res_df[order(res_df$padj), ])[1:30]

# Use variance-stabilizing transformation (better for visualization)
vsd <- vst(dds, blind = FALSE)

# Extract matrix of normalized expression for top genes
top_counts <- assay(vsd)[top_genes, ]
annotation_col <- as.data.frame(colData(dds)[, c("dex", "cell")])
pheatmap(top_counts,
         scale = "row",
         annotation_col = annotation_col) 

## Plot dispersion estimates
plotDispEsts(dds)

# Perform variance-stabilizing transformation (VST) for visualization
vsd <- vst(dds, blind = FALSE)

# Plot PCA using DESeq2's built-in function 
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "dex") + 
  ggtitle("PCA: trt vs untrt") +
  theme_minimal()
