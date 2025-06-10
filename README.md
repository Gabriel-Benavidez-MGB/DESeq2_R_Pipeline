# DESeq2_R_Pipeline

This repository contains a simple and reproducible pipeline for performing Differential Gene Expression (DGE) analysis in R using the **DESeq2** package and the publicly available `airway` dataset.

---

## 📁 Folder Structure

DESeq2_R_Pipeline/
├── scripts/
│ └── airway_deseq2.R # Main analysis script
├── data/
│ ├── airway_counts.csv # Gene expression count matrix
│ └── airway_metadata.csv # Sample metadata
├── results/
│ ├── Dispersion_Plot.pdf
│ ├── HeatMap.pdf
│ ├── MA_Plot.pdf
│ ├── PCA_Plot.pdf
│ └── Volcano Plot.pdf
---
## 📋 Description

The pipeline performs the following steps:

1. Loads count data and metadata.
2. Constructs a `DESeqDataSet` using DESeq2.
3. Runs the differential expression analysis.
4. Generates visualizations:
   - **Volcano plot**
   - **MA plot**
   - **PCA plot**
   - **Heatmap**
   - **Dispersion plot**
---
## 🔧 Requirements
- R (≥ 4.0)
- `DESeq2`
- `tidyverse`
- `pheatmap`
- `ggplot2`
