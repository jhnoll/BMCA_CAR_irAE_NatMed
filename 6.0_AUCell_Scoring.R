# Load necessary libraries ----
library(Seurat)        # Core toolkit for single-cell RNA-seq analysis
library(tidyverse)     # Includes dplyr, ggplot2, etc.
library(ggpubr)        # Enhances ggplot2 for publication
library(limma)         # Linear modeling for RNA-seq/microarrays
library(SingleR)       # Reference-based cell type annotation
library(enrichR)       # Enrichment analysis via Enrichr
library(cowplot)       # Plot composition and layout
library(patchwork)     # Combine ggplots with operators
library(VennDiagram)   # Venn diagrams
library(pheatmap)      # Pretty heatmaps
library(RColorBrewer)  # Color palettes
library(msigdbr)       # Gene sets from MSigDB
library(celldex)       # Reference datasets for SingleR
library(writexl)       # Write to Excel
library(ggrepel)       # Better text labels in ggplot2
library(AUCell)        # AUCell scoring of pathway activity
library(GSEABase)      # For working with gene sets
library(Nebulosa)      # Density plots on UMAPs

# Set working directories ----
base <- "C:/Users/julia/Box Sync/3 PhD Penn/2 Thesis/1 Projects/JHN_081_Analysis of 37418-251 patient/JHN081_scRNAseq"
objSaveDir <- file.path(base, "seurat_objects")
fig_base   <- file.path(base, "new_figures")
gene_dir   <- file.path(base, "gene_list")

# Step 0: Load processed Seurat object ----
setwd(objSaveDir)
clean.seurat.1 <- readRDS("CleanSeurat.RDS")

# Step 1: Load gene sets from GMT files ----
gmt_files <- list.files(path = gene_dir, pattern = "\\.gmt$", full.names = TRUE)

named_gene_sets <- lapply(gmt_files, function(f) {
  gs <- getGmt(f)
  name_prefix <- tools::file_path_sans_ext(basename(f))  # get file base name
  names(gs) <- paste0(name_prefix)                       # assign set name
  return(gs)
})

# Combine all gene sets into a GeneSetCollection
all_gene_sets <- do.call(c, named_gene_sets)
all_gene_sets <- GeneSetCollection(unlist(all_gene_sets))

# Optional: Add custom gene sets here (if defined)
# all_gene_sets <- GeneSetCollection(c(as.list(all_gene_sets), as.list(custom_gsc)))

# Step 2: Compute AUCell pathway activity scores ----

# Get RNA assay matrix
exprMatrix <- as.matrix(GetAssayData(clean.seurat.1, assay = "RNA", slot = "data"))

# Build gene rankings per cell
cells_rankings <- AUCell_buildRankings(exprMatrix)

# Calculate AUC scores for each gene set across all cells
cells_AUC <- AUCell_calcAUC(all_gene_sets, cells_rankings)

# Step 3: Add AUCell scores back to Seurat metadata ----
auc_matrix <- as.data.frame(t(getAUC(cells_AUC)))
clean.seurat.1 <- AddMetaData(clean.seurat.1, auc_matrix)

# Step 4: Density plots comparing PBMC vs CSF (per gene set) ----

# Subset by compartment
pbmc_cells <- subset(clean.seurat.1, subset = pbmc.or.csf == "PBMC")
csf_cells  <- subset(clean.seurat.1, subset = pbmc.or.csf == "CSF")

# List of gene sets
lists <- colnames(auc_matrix)

# Plot and save density plots side-by-side for each gene set
for (list in lists) {
  p_pbmc <- plot_density(pbmc_cells, features = list) + ggtitle(paste("PBMC:", list))
  p_csf  <- plot_density(csf_cells,  features = list) + ggtitle(paste("CSF:", list))
  
  combined_plot <- p_pbmc + p_csf + plot_layout(ncol = 2)
  
  ggsave(
    filename = paste0("AUCell_PBMCvsCSF_", list, ".pdf"),
    plot = combined_plot,
    width = 10, height = 4, dpi = 300
  )
}

# Step 5: Save updated Seurat object ----
saveRDS(clean.seurat.1, file = "CleanSeurat_AUCell.RDS")
clean.seurat.1 <- readRDS("CleanSeurat_AUCell.RDS")

# Step 6: Summarize AUCell scores by cluster and compartment ----

auc_summary <- clean.seurat.1@meta.data %>%
  as.data.frame() %>%
  group_by(cluster_annotation, pbmc.or.csf) %>%
  summarise(across(
    .cols = all_of(colnames(auc_matrix)),
    .fns = mean,
    .names = "avg_{.col}"
  ), .groups = "drop")

write.csv(auc_summary, file = "new_figures/AUCell score by T cell clusters and CSF and PBMC.csv")

# Step 7: Boxplot comparisons for one cluster only ----

# Define the gene sets to plot
pathways <- colnames(clean.seurat.1@meta.data)[11:29]

# Filter to specific cluster (e.g., "CSF-Specific CD4 T cells")
plot_data <- clean.seurat.1@meta.data %>%
  filter(cluster_annotation == "CSF-Specific CD4 T cells") %>%
  dplyr::select(pbmc.or.csf, all_of(pathways)) %>%
  pivot_longer(cols = all_of(pathways), names_to = "Pathway", values_to = "Score")

# Plot boxplots per pathway, split by PBMC vs CSF
ggplot(plot_data, aes(x = pbmc.or.csf, y = Score, fill = pbmc.or.csf)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.4) +
  facet_wrap(~ Pathway, scales = "free_y", ncol = 4) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     comparisons = list(c("PBMC", "CSF")),
                     label.y.npc = "top", size = 3) +
  labs(
    title = "Pathway Scores in CSF-specific T cells (PBMC vs CSF)",
    x = "Compartment", y = "Score"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Step 8: Per-pathway boxplots across all clusters ----

# Create output folder
out_dir <- "AUCell_Plots"
dir.create(out_dir, showWarnings = FALSE)

for (pathway in pathways) {
  plot_data <- clean.seurat.1@meta.data %>%
    dplyr::select(pbmc.or.csf, cluster_annotation, !!sym(pathway)) %>%
    dplyr::rename(Score = !!sym(pathway))  # Dynamically rename
  
  p <- ggplot(plot_data, aes(x = pbmc.or.csf, y = Score, fill = pbmc.or.csf)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
    facet_wrap(~ cluster_annotation, scales = "free_y", ncol = 3) +
    stat_compare_means(
      method = "wilcox.test", label = "p.format", label.y.npc = "top", size = 2.8
    ) +
    labs(
      title = paste0("AUCell Score: ", pathway),
      x = "Compartment", y = "AUCell Score"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  ggsave(file.path(out_dir, paste0(pathway, "_AUCell_boxplot.pdf")),
         plot = p, width = 8.5, height = 6)
}

# Step 9: Boxplots for CSF-only cells by cluster ----

plot_data <- clean.seurat.1@meta.data %>%
  filter(pbmc.or.csf == "CSF") %>%
  dplyr::select(cluster_annotation, all_of(pathways)) %>%
  pivot_longer(cols = all_of(pathways), names_to = "Pathway", values_to = "Score")

output_dir <- "CSF_Pathway_Plots"
dir.create(output_dir, showWarnings = FALSE)

unique_pathways <- unique(plot_data$Pathway)

for (p in unique_pathways) {
  plot <- ggplot(filter(plot_data, Pathway == p), aes(x = cluster_annotation, y = Score)) +
    geom_boxplot(outlier.shape = NA, fill = "gray80") +
    geom_jitter(width = 0.2, size = 0.8, alpha = 0.6) +
    theme_minimal(base_size = 12) +
    labs(
      title = p,
      x = "Cluster Annotation", y = "Pathway Score"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(filename = file.path(output_dir, paste0(gsub("[^a-zA-Z0-9]", "_", p), ".pdf")),
         plot = plot, width = 6, height = 4)
}
