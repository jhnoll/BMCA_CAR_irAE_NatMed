# Load required libraries
library(Seurat)        # Core toolkit for single-cell RNA-seq analysis
library(tidyverse)     # Data manipulation and plotting
library(ggpubr)        # Publication-ready ggplot2 enhancements
library(limma)         # Linear modeling (originally for microarrays)
library(SingleR)       # Cell type annotation using reference datasets
library(enrichR)       # Enrichment analysis via Enrichr API
library(cowplot)       # Flexible plot composition
library(patchwork)     # Combine ggplot-based plots
library(VennDiagram)   # Venn diagram visualization
library(pheatmap)      # Pretty heatmaps
library(RColorBrewer)  # Color palettes
library(msigdbr)       # Gene set collections for MSigDB
library(celldex)       # Reference data for SingleR
library(dplyr)
library(ggplot2)
library(writexl)       # Excel output

# -------------------------------
# Set working directories & load data
# -------------------------------
base <- "C:/Users/julia/Box Sync/3 PhD Penn/2 Thesis/1 Projects/JHN_081_Analysis of 37418-251 patient/JHN081_scRNAseq"
objSaveDir <- file.path(base, "seurat_objects")
fig_base <- file.path(base, "figures")

setwd(objSaveDir)
Tcells.filtered.ciltacel.11 <- readRDS("Tcells.filtered.ciltacel.11_integration_AUCell.RDS")

# -------------------------------
# Step 1: Subset to clusters of interest
# -------------------------------
# Keep clusters 0–8 and 11 for further refinement
clusters_keep <- c(0,1,2,3,4,5,6,7,8,11)
clean.seurat <- subset(Tcells.filtered.ciltacel.11, subset = seurat_clusters %in% clusters_keep)

# Visualize composition by compartment
DimPlot(clean.seurat, split.by = "pbmc.or.csf")

# -------------------------------
# Step 2: Re-run PCA, UMAP, and clustering
# -------------------------------
clean.seurat.1 <- ScaleData(clean.seurat, verbose = FALSE)
clean.seurat.1 <- RunPCA(clean.seurat.1, npcs = 15, verbose = FALSE)
clean.seurat.1 <- RunUMAP(clean.seurat.1, reduction = "pca", dims = 1:15)
clean.seurat.1 <- FindNeighbors(clean.seurat.1, reduction = "pca", dims = 1:15, k.param = 15)
clean.seurat.1 <- FindClusters(clean.seurat.1, resolution = 0.2)

# -------------------------------
# Step 3: Refine to only clusters 0–3
# -------------------------------
clusters_keep <- c(0,1,2,3)
clean.seurat.1 <- subset(clean.seurat.1, subset = seurat_clusters %in% clusters_keep)

# Visualize new clustering results
DimPlot(clean.seurat.1, split.by = "pbmc.or.csf")
ggsave("new_figures/Composition by CSF and PBMC.pdf", width = 16, height = 7)

DimPlot(clean.seurat.1, split.by = "seurat_clusters")

# -------------------------------
# Step 4: Plot gene expression across clusters
# -------------------------------
FeaturePlot(clean.seurat.1, features = c("rna_CD3E", "rna_CD4", "rna_CD8A", "ciltacel"), ncol = 2)

# Lineage and activation markers
VlnPlot(clean.seurat.1, features = c("rna_CD3E", "rna_CD4", "rna_CD8A", "rna_CD27"),
        group.by = "seurat_clusters", pt.size = 0)

# Myeloid contamination check
VlnPlot(clean.seurat.1, features = c("rna_CD14", "rna_CD163", "rna_MPO"),
        group.by = "seurat_clusters", pt.size = 0)

# Cytotoxicity and proliferation
VlnPlot(clean.seurat.1, features = c("rna_GZMB", "rna_GZMA", "rna_GZMK", "rna_GZMH", "rna_KI67"),
        group.by = "seurat_clusters", pt.size = 0)

# QC metrics and ciltacel detection
VlnPlot(clean.seurat.1, features = c("ciltacel", "nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "seurat_clusters", pt.size = 0)

# -------------------------------
# Step 5: Identify subpopulations within cluster 2 based on CD8A expression
# -------------------------------
clean.seurat.1_cluster2 <- WhichCells(clean.seurat.1, idents = 2)
clean.seurat.1_cluster2 <- clean.seurat.1[, clean.seurat.1_cluster2]

# Visualize cluster 2 subset
DimPlot(clean.seurat.1_cluster2)

# Split cluster 2 cells based on CD8A scaled expression
cluster2.8 <- colnames(clean.seurat.1_cluster2@assays$integrated$scale.data)[
  clean.seurat.1_cluster2@assays$integrated$scale.data["CD8A", ] >= 1]
cluster2.4 <- colnames(clean.seurat.1_cluster2@assays$integrated$scale.data)[
  clean.seurat.1_cluster2@assays$integrated$scale.data["CD8A", ] < 1]

# Subset to new groups for visualization
cluster2.8 <- clean.seurat.1[, cluster2.8]
cluster2.4 <- clean.seurat.1[, cluster2.4]

# Visualize CD4/CD8 in CD8A-high cluster
FeaturePlot(cluster2.8, features = c("rna_CD4", "CD8A"), ncol = 2)

# -------------------------------
# Step 6: Manually reassign cluster identities
# -------------------------------
# Reclassify CD8-high cells from cluster 2 as cluster 4
cluster4_cells <- colnames(cluster2.8)

idents <- Idents(clean.seurat.1)
levels(idents) <- c(levels(idents), "4")  # Add cluster 4 label
idents[cluster4_cells] <- "4"

clean.seurat.1 <- SetIdent(clean.seurat.1, value = idents)
table(Idents(clean.seurat.1))
ncol(clean.seurat.1)

# -------------------------------
# Step 7: Final visualizations
# -------------------------------

# UMAP with final cluster assignments
DimPlot(clean.seurat.1, label = TRUE)
ggsave("new_figures/DimPlot_ciltacel_seurat_clusters.pdf", width = 9, height = 7)

# UMAP by sample
DimPlot(clean.seurat.1, split.by = "orig.ident")
ggsave("new_figures/DimPlot_ciltacel_by_sample.pdf", width = 20, height = 7)

# UMAP by compartment
DimPlot(clean.seurat.1, split.by = "pbmc.or.csf", label = TRUE)
ggsave("new_figures/DimPlot_ciltacel_PBMC-or-CSF.pdf", width = 13, height = 7)

# Feature plots for CD4/CD8
FeaturePlot(clean.seurat.1, features = c("rna_CD4", "rna_CD8A"), ncol = 2)
ggsave("new_figures/FeaturePlot_ciltacel_CD4_CD8A.pdf", width = 8, height = 7)

# UMAPs by CD4/CD8 grouping (e.g., manual annotation column)
DimPlot(clean.seurat.1, group.by = "CD4vsCD8")
ggsave("new_figures/DimPlot_ciltacel_CD4vsCD8.pdf", width = 9, height = 7)

# UMAP colored by PBMC vs CSF with custom colors
DimPlot(clean.seurat.1, group.by = "pbmc.or.csf", cols = c("purple", "orange"))
ggsave("new_figures/DimPlot_ciltacel_groupPBMC-or-CSF.pdf", width = 13, height = 7)
