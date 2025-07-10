# Load necessary libraries ----
library(Seurat)        # Single-cell RNA-seq analysis toolkit
library(tidyverse)     # Data manipulation and plotting
library(ggpubr)        # Enhanced ggplot2 visualizations
library(limma)         # Linear modeling
library(SingleR)       # Automated cell type annotation
library(enrichR)       # Gene set enrichment analysis
library(cowplot)       # Combining multiple ggplots
library(patchwork)     # Plot layout management
library(VennDiagram)   # Venn diagram generation
library(pheatmap)      # Heatmaps
library(RColorBrewer)  # Custom color palettes
library(msigdbr)       # MSigDB gene set database
library(celldex)       # Reference datasets for annotation
library(writexl)       # Writing to Excel files
library(ggrepel)       # Text label repulsion in plots
library(AUCell)        # Gene set scoring per cell
library(GSEABase)      # Handling gene set objects
library(Nebulosa)      # Expression density on embeddings

# Set working directories ----
base <- "C:/Users/julia/Box Sync/3 PhD Penn/2 Thesis/1 Projects/JHN_081_Analysis of 37418-251 patient/JHN081_scRNAseq"
objSaveDir <- file.path(base, "seurat_objects")
fig_base   <- file.path(base, "new_figures")
gene_dir   <- file.path(base, "gene_list")

# Step 0: Load Seurat object with AUCell scores ----
setwd(objSaveDir)
clean.seurat.1 <- readRDS("CleanSeurat_AUCell.RDS")

# Step 1: Identify marker genes for each cluster ----
# Parameters:
# - only.pos = TRUE : only returns upregulated genes
# - min.pct = 0.25  : genes expressed in at least 25% of cells in the cluster
# - logfc.threshold = 0.25 : minimum log fold-change

cluster0markers <- FindMarkers(clean.seurat.1, ident.1 = 0, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster0markers, file = "DEGs_Cluster0.csv")

cluster1markers <- FindMarkers(clean.seurat.1, ident.1 = 1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster1markers, file = "DEGs_Cluster1.csv")

cluster2markers <- FindMarkers(clean.seurat.1, ident.1 = 2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster2markers, file = "DEGs_Cluster2.csv")

cluster3markers <- FindMarkers(clean.seurat.1, ident.1 = 3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster3markers, file = "DEGs_Cluster3.csv")

cluster4markers <- FindMarkers(clean.seurat.1, ident.1 = 4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster4markers, file = "DEGs_Cluster4.csv")

# Step 2: DotPlot of functional modules across clusters ----
# Shows expression of key marker sets grouped by function (naive, cytotoxic, exhausted, etc.)

DotPlot(
  clean.seurat.1,
  features = list(
    "CD4/CD8" = c("rna_CD4", "rna_CD8A", "rna_CD8B"),
    "Naive" = c("CCR7", "SELL", "LEF1", "TCF7", "IL7R"),
    "Memory" = c("rna_CD27", "rna_CD44", "rna_CD28"),
    "Cytotoxic" = c("GZMB", "PRF1", "GNLY", "IFNG", "KLRG1"),
    "Exhausted" = c("PDCD1", "CTLA4", "TOX", "LAG3"),
    "Proliferating" = c("MKI67"),
    "Activation" = c("IL2RA", "CD69", "HLA-DRA", "TNFRSF4", "TNFRSF9", "ICOS"),
    "Migration" = c("S1PR1", "CXCR3", "CXCR4", "CCR5", "ITGAE", "ITGB7"),
    "Innate/Inflammatory" = c("TLR4", "IL1R2", "CSF2", "CD33")
  ),
  dot.min = 0.25,
  cols = c("grey", "blue")  # expression intensity scale
) +
  RotatedAxis()

# You may want to save this plot:
ggsave("DotPlot_TcellMarkers.pdf", width = 12, height = 6)

