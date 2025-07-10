# Load necessary libraries ----
library(Seurat)        # Single-cell RNA-seq analysis
library(tidyverse)     # Data wrangling and visualization (includes ggplot2, dplyr, etc.)
library(ggpubr)        # Publication-quality ggplot2 enhancements
library(limma)         # Linear models for expression analysis
library(SingleR)       # Automated cell type annotation
library(enrichR)       # Enrichment analysis
library(cowplot)       # Plot composition
library(patchwork)     # Combine plots
library(VennDiagram)   # Venn diagram generation
library(pheatmap)      # Heatmap visualization
library(RColorBrewer)  # Color palettes
library(msigdbr)       # MSigDB gene sets
library(celldex)       # SingleR references
library(dplyr)
library(ggplot2)
library(writexl)       # Write Excel files
library(ggrepel)       # Repel overlapping text
library(AUCell)        # Gene set activity via AUC
library(GSEABase)      # Gene set objects
library(Nebulosa)      # Density-based plots on UMAPs

# Set working directories ----
base <- "C:/Users/julia/Box Sync/3 PhD Penn/2 Thesis/1 Projects/JHN_081_Analysis of 37418-251 patient/JHN081_scRNAseq"
objSaveDir <- file.path(base, "seurat_objects")
fig_base   <- file.path(base, "new_figures")
gene_dir   <- file.path(base, "gene_list")

# Step 0: Load Seurat object with AUCell scores and annotated clusters ----
setwd(objSaveDir)
clean.seurat.1 <- readRDS("CleanSeurat.RDS")

# Step 1: UMAP with annotated clusters ----
DimPlot(clean.seurat.1, group.by = "cluster_annotation")
ggsave("new_figures/annotatedDimPlot.pdf", width = 8, height = 8)

# Step 2: Define and apply new cluster labels ----
# If identities were not previously annotated, manually assign using current Idents
cluster.labels <- c(
  "0" = "Activated CD4 T cells",
  "1" = "Cytotoxic CD8 T cells",
  "2" = "Proliferating CD4 T cells",
  "3" = "Interferon-stimulated T-rm like CD4 T cells",
  "4" = "Proliferating CD8 T cells"
)

# Add new annotation column based on Seurat cluster ID
clean.seurat.1$cluster_annotation <- case_when(
  Idents(clean.seurat.1) == 0 ~ cluster.labels["0"],
  Idents(clean.seurat.1) == 1 ~ cluster.labels["1"],
  Idents(clean.seurat.1) == 2 ~ cluster.labels["2"],
  Idents(clean.seurat.1) == 3 ~ cluster.labels["3"],
  Idents(clean.seurat.1) == 4 ~ cluster.labels["4"]
)

# Step 3A: Bar plot — composition by compartment (PBMC vs CSF) ----
abundance_table <- clean.seurat.1@meta.data %>%
  as.data.frame() %>%
  group_by(pbmc.or.csf, cluster_annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pbmc.or.csf) %>%
  mutate(freq = n / sum(n))

ggplot(abundance_table, aes(x = pbmc.or.csf, y = freq, fill = cluster_annotation)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.75) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(),
    legend.key = element_rect(fill = NA),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ylab("Proportion of Cells") + xlab("Time") +
  ggtitle("Composition by CSF and PBMC") +
  scale_y_continuous(labels = scales::percent_format())

ggsave("new_figures/T cell composition by CSF and PBMC.pdf", width = 15, height = 12)
write.csv(abundance_table, "new_figures/T cell composition by CSF and PBMC.csv")

# Step 3B: Bar plot — flipped orientation (cluster × PBMC/CSF) ----
abundance_table <- clean.seurat.1@meta.data %>%
  as.data.frame() %>%
  group_by(cluster_annotation, pbmc.or.csf) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster_annotation) %>%
  mutate(freq = n / sum(n))

ggplot(abundance_table, aes(x = cluster_annotation, y = freq, fill = pbmc.or.csf)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.75) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(),
    legend.key = element_rect(fill = NA),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylab("Proportion of Cells") + xlab("T Cell Subtype") +
  ggtitle("Composition by CSF and PBMC") +
  scale_y_continuous(labels = scales::percent_format())

ggsave("new_figures/Composition by CSF and PBMC.pdf", width = 15, height = 12)
write.csv(abundance_table, "new_figures/Composition by CSF and PBMC.csv")

# Step 3C: Bar plot — composition by sample and cluster ----
abundance_table <- clean.seurat.1@meta.data %>%
  as.data.frame() %>%
  group_by(orig.ident, cluster_annotation, pbmc.or.csf) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster_annotation) %>%
  mutate(freq = n / sum(n))

ggplot(abundance_table, aes(x = orig.ident, y = freq, fill = cluster_annotation)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.75) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(),
    legend.key = element_rect(fill = NA),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylab("Proportion of Cells") + xlab("Sample") +
  ggtitle("Composition by Sample and PBMC/CSF") +
  scale_y_continuous(labels = scales::percent_format())

ggsave("new_figures/Composition by CSF and PBMC.pdf", width = 15, height = 12)
write.csv(abundance_table, "new_figures/Composition of samples by T cell clusters and CSF and PBMC.csv")
