# Load necessary libraries ----
library(Seurat)        # Single-cell RNA-seq toolkit
library(tidyverse)     # Data wrangling and visualization (ggplot2, dplyr, etc.)
library(ggpubr)        # Publication-ready ggplots
library(limma)         # Differential expression modeling
library(SingleR)       # Reference-based cell type annotation
library(enrichR)       # Enrichment analysis via Enrichr
library(cowplot)       # Layout control for ggplots
library(patchwork)     # Combine ggplot2 objects using operators
library(VennDiagram)   # Venn diagram generation
library(pheatmap)      # Heatmap visualization
library(RColorBrewer)  # Color palettes
library(msigdbr)       # Gene sets from MSigDB
library(celldex)       # SingleR reference datasets
library(dplyr)
library(ggplot2)
library(writexl)       # Write Excel files
library(ggrepel)       # Repel overlapping ggplot labels
library(AUCell)        # AUCell scoring of gene set activity
library(GSEABase)      # Gene set representation
library(Nebulosa)      # Density plotting on UMAPs
library(ggpubr)
library(pheatmap)

# Set working directories ----
base <- "C:/Users/julia/Box Sync/3 PhD Penn/2 Thesis/1 Projects/JHN_081_Analysis of 37418-251 patient/JHN081_scRNAseq"
objSaveDir <- file.path(base, "seurat_objects")
fig_base   <- file.path(base, "new_figures")
gene_dir   <- file.path(base, "gene_list")

# Step 0: Load integrated Seurat object (AUCell-scored) ----
setwd(objSaveDir)
clean.seurat.1 <- readRDS("CleanSeurat.RDS")

# Step 1: Define key marker features for density plotting and heatmaps ----
all_features <- c(
  # Cytotoxicity
  "rna_GZMA", "rna_GZMB", "rna_GZMH", "rna_GZMK", "rna_PRF1", "rna_NKG7", "rna_LAMP1", "rna_IFNG",
  
  # Activation / Costimulation
  "rna_CD69", "rna_CD38", "rna_CD27", "rna_TNFRSF9",
  
  # Proliferation
  "rna_MKI67", "rna_TYMS", "rna_PCNA",
  
  # Naive / Memory / Homing
  "rna_CCR7", "rna_SELL", "rna_IL7R",
  
  # Inhibitory Receptors
  "rna_PDCD1", "rna_LAG3", "rna_TIGIT", "rna_CTLA4",
  
  # T Cell Lineage / Transcription Factors
  "rna_CD3E", "rna_CD4", "rna_CD8A", "rna_CD8B",
  "rna_RUNX3", "rna_TBX21", "rna_EOMES",
  
  # Interferon Response
  "rna_IFIT1", "rna_IFIT3", "rna_ISG15", "rna_OAS1", "rna_OAS3", "rna_MX1",
  "rna_CXCL9", "rna_CXCL10", "rna_CXCL11", "rna_IRF1", "rna_IRF7", "rna_STAT1",
  
  # Tissue-Resident Memory (Trm) Markers
  "rna_ITGAE", "rna_ITGA1", "rna_CXCR6", "rna_CD101",
  "rna_ZNF683", "rna_BHLHE40", "rna_PRDM1",
  "rna_S1PR1", "rna_KLF2",
  
  # Chemokines / Receptors
  "rna_CCL3", "rna_CCL4", "rna_CCL5", "rna_CXCL13",
  "rna_CCR5", "rna_CXCR3", "rna_CXCR4", "rna_CXCR5", "rna_IL15",
  
  # Cytokine Signaling
  "rna_IL15", "rna_IL15RA"
)

# Step 2: Generate density plots in chunks of 12 genes per page ----
n_per_page <- 12
n_pages <- ceiling(length(all_features) / n_per_page)

for (i in seq_len(n_pages)) {
  features_subset <- all_features[((i - 1) * n_per_page + 1):min(i * n_per_page, length(all_features))]
  
  pdf(paste0("density_plot_page_", i, ".pdf"), width = 12, height = 8)
  print(plot_density(clean.seurat.1, features = features_subset, joint = FALSE))
  dev.off()
}

# Step 3: Heatmap of average expression across clusters ----
DefaultAssay(clean.seurat.1) <- "RNA"
Idents(clean.seurat.1) <- clean.seurat.1@meta.data$cluster_annotation

# Remove "rna_" prefix for readability in heatmap
all_features_stripped <- gsub("^rna_", "", all_features)

# Calculate average expression per cluster (log-normalized)
avg_expr <- AverageExpression(
  clean.seurat.1,
  features = all_features_stripped,
  group.by = "cluster_annotation",
  slot = "data"
)$RNA

# Z-score scale across rows (genes)
scaled_expr <- t(scale(t(avg_expr)))

# Plot heatmap
pheatmap(
  scaled_expr,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize = 9,
  cellwidth = 30,
  cellheight = 10,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
)

# Step 4: Targeted density plots comparing PBMC vs CSF ----

# Subset the Seurat object by sample origin
pbmc_cells <- subset(clean.seurat.1, subset = pbmc.or.csf == "PBMC")
csf_cells  <- subset(clean.seurat.1, subset = pbmc.or.csf == "CSF")

# Define additional cytotoxic / exhaustion / Trm genes to visualize
genes <- c(
  "rna_GZMA", "rna_GZMB", "rna_GZMK", "rna_GZMM", "rna_GZMH",
  "rna_PRF1", "rna_GNLY",
  "rna_IFNG", "rna_TNF", "rna_IL2",
  "rna_TBX21", "rna_EOMES",
  "rna_FASLG", "rna_TNFSF10", "rna_KLRK1", "rna_KLRD1", "rna_KLRG1",
  "rna_LAMP1", "rna_RAB27A",
  "rna_TOX", "rna_FOXP3"
)

# Loop through and create side-by-side density plots for PBMC vs CSF
for (gene in genes) {
  p_pbmc <- plot_density(pbmc_cells, features = gene) + ggtitle(paste("PBMC:", gene))
  p_csf  <- plot_density(csf_cells,  features = gene) + ggtitle(paste("CSF:", gene))
  
  combined_plot <- p_pbmc + p_csf + plot_layout(ncol = 2)
  
  ggsave(
    filename = paste0("PBMCvsCSF_", gene, ".pdf"),
    plot = combined_plot,
    width = 10, height = 4, dpi = 300
  )
}

# Step 5: Save overall gene density plots in one PDF ----
plot_density(clean.seurat.1, features = genes)
ggsave(filename = "DensityPlots.pdf", width = 20, height = 20)

# Step 6: Save processed object ----
saveRDS(clean.seurat.1, file = "CleanSeurat.RDS")

# (Optional) Load it again later
clean.seurat.1 <- readRDS("CleanSeurat.RDS")
