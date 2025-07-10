# Load required libraries ----
library(flowCore)        # Core package for flow cytometry data handling
library(flowWorkspace)   # Workspace management and analysis for flow cytometry
library(Seurat)          # Toolkit for single-cell genomics, applicable to flow cytometry
library(ggplot2)         # Data visualization package
library(tidyverse)       # Collection of data manipulation and visualization tools
library(dplyr)           # Data manipulation functions
library(FlowSOM)        # FlowSOM clustering and visualization
library(reshape2)        # Reshaping data
library(pheatmap)        # Heatmap visualization
library(wadeTools)       # Specific tools for flow cytometry data

# Step 0: Define directories ----
proj_base = ""
# Adjust the username for work computer

# Set up directories for data and figures
data_base = tight(proj_base, "/FlowJo_Excel/CD3")  # FlowJo CSV export directory
fig_base = tight(proj_base, "/figures/CD3")       # Directory for saving pictures
int_base = tight(proj_base,"/intermediate_subjects/CD3")  # Directory for intermediate files

# Step 1: Run FlowSOM ----
setwd(proj_base)  # Set working directory to the project base

set.seed(0)  # Set seed for reproducibility

# Perform FlowSOM clustering on the transformed data
fSOM.15 <- FlowSOM(
  t(as.matrix(seurat.concat@assays$flow$data)),  # Transpose and convert Seurat data to matrix
  compensate = FALSE,   # No compensation for fluorescence
  transform = FALSE,    # No transformation applied
  scale = FALSE,        # No scaling applied
  nClus = 15           # Number of clusters
)

# Add FlowSOM clustering results to Seurat object
seurat.concat[["FlowSOM15_clusters"]] <- as.character(GetClusters(fSOM.15))
seurat.concat[["FlowSOM_metaclusters_15"]] <- GetMetaclusters(fSOM.15)
# Plot UMAP colored by FlowSOM metaclusters
DimPlot(seurat.concat, group.by = "FlowSOM_metaclusters_15", label = TRUE, raster = TRUE)
ggsave("UMAP_FlowSOM.pdf", width = 8, height = 6)  # Save UMAP plot
# Save FlowSOM results for future use
saveRDS(fSOM.15, file = "fSOM.15.RDS")
fSOM.15 <- readRDS(file = "intermediate_subjects/fSOM.15.RDS")

# Step 2: Analyze Cluster Proportions ----
# Calculate cluster proportions for each sample
fSOM.15_cluster_proportions <- seurat.concat@meta.data %>% 
  group_by(SampleName, FlowSOM_metaclusters_15) %>% 
  summarise(n = n()) %>%  # Count cells in each cluster
  mutate(freq = n / sum(n))  # Calculate frequency of each cluster

# Rename clusters for clarity
levels(fSOM.15_cluster_proportions$FlowSOM_metaclusters_15) <- paste0("fSOM.15_cluster", levels(fSOM.15_cluster_proportions$FlowSOM_metaclusters_15))

# Reshape data for heatmap
fSOM.15_cluster_proportions <- dcast(fSOM.15_cluster_proportions, SampleName ~ FlowSOM_metaclusters_15, value.var = "freq")
fSOM.15_cluster_proportions[is.na(fSOM.15_cluster_proportions)] <- 0
write.csv(fSOM.15_cluster_proportions, "fSOM_15_cluster_proportions.csv", row.names = TRUE)

fSOM.15_cluster_proportions <-read.csv("fSOM_15_cluster_proportions.csv")

rownames(fSOM.15_cluster_proportions) <- fSOM.15_cluster_proportions$SampleName
fSOM.15_cluster_proportions <- fSOM.15_cluster_proportions[,-1]
pheatmap(fSOM.15_cluster_proportions, scale = "column", cluster_rows = FALSE)

# Step 3: Analyze Cluster MFI (Mean Fluorescence Intensity) ----
# Calculate mean fluorescence intensity for each cluster
fSOM.15_cluster_MFI <- GetMetaclusterMFIs(fSOM.15, colsUsed = FALSE, prettyColnames = FALSE)

# Plot heatmap of MFI values
pheatmap(fSOM.15_cluster_MFI)
saveRDS(fSOM.15_cluster_MFI, file = "intermediate_subjects/fSOM.15_cluster_MFI.RDS")

# Reload MFI data and scale for better visualization
fSOM.15_cluster_MFI_scaled_t <- t(scale(t(fSOM.15_cluster_MFI)))
pheatmap(fSOM.15_cluster_MFI_scaled_t, cluster_rows = FALSE)
write.csv(fSOM.15_cluster_MFI_scaled_t, "fSOM_15_cluster_MFI.csv", row.names = TRUE)


