# Load required libraries ----
library(flowCore)        # Core package for flow cytometry data handling
library(flowWorkspace)   # Workspace management and analysis for flow cytometry
library(Seurat)          # Toolkit for single-cell genomics, applicable to flow cytometry
library(Matrix)          # Sparse and dense matrix operations
library(ggplot2)         # Data visualization package
library(tidyverse)       # Collection of data manipulation and visualization tools
library(dplyr)           # Data manipulation functions
library(gtools)          # Assorted R programming tools
library(wadeTools)       # Specific tools for flow cytometry data
library(gridExtra)       # Functions to arrange multiple plots

# Step 0: Import Data ----
# Define directories for the project
proj_base = ""

# Set up directories for data and figures
data_base = tight(proj_base, "/FlowJo_Excel/CD3")  # FlowJo CSV export directory
fig_base = tight(proj_base, "/figures/CD3")       # Directory for saving pictures
int_base = tight(proj_base,"/intermediate_subjects/CD3")  # Directory for intermediate files

setwd(int_base)  # Set working directory to int directory

# Load transformed data matrix and metadata from previous steps
transformed_data.t <- readRDS(file = "flowMatrix_transformed.rds")
metadata <- readRDS(file = "metadata.rds")

# Step 1: Create Seurat Object ----
# Create a Seurat object with filtered data and associated metadata
seurat.concat <- CreateSeuratObject(counts = transformed_data.t, data = transformed_data.t, assay = "flow", meta.data = metadata)
seurat.concat[["all.cells"]] <- "Cell"  # Add an identifier for all cells
# Set variable features for further analysis
VariableFeatures(seurat.concat) <- rownames(seurat.concat@assays$flow)

# Step 2: Scale Data ----
# Scale the data (center and scale across features)
seurat.concat@assays$flow$scale.data <- as.matrix(t(scale(t(as.matrix(seurat.concat@assays$flow$counts)), center = TRUE, scale = TRUE)))

# Step 3: Perform UMAP for Visualization and Clustering ----
# Run PCA for dimensionality reduction
seurat.concat <- RunPCA(seurat.concat, verbose = TRUE)
# Identify neighbors and clusters based on PCA dimensions
seurat.concat <- FindNeighbors(seurat.concat, dims = 1:5, k.param = 20, verbose = TRUE)
seurat.concat <- FindClusters(seurat.concat, verbose = TRUE)
# Perform UMAP for visualization, setting the number of neighbors and dimensions
seurat.concat <- RunUMAP(seurat.concat, dims = 1:5, verbose = TRUE, n.neighbors = 20)

# Save UMAP plot
setwd(fig_base)  # Change working directory to figures directory
DimPlot(seurat.concat, reduction = "umap", pt.size = 1, label = TRUE, raster = TRUE)  # Plot UMAP
ggsave("UMAP.pdf", width = 8, height = 6)  # Save UMAP plot

rownames(seurat.concat) <- sub(".*\\.\\.\\.\\.", "", rownames(seurat.concat))

FeaturePlot(seurat.concat, features = rownames(seurat.concat), reduction = "umap", pt.size = 1, min.cutoff = "q05", max.cutoff = "q95", ncol = 4, raster=TRUE)
ggsave("UMAP_FeaturePlot.pdf", width = 7 * 4 * 0.75, height = 8 * 3 * 0.75)  # Save feature plot

# Step 4: Create Density Plot for UMAP ----
# Extract metadata for plotting
dat <- data.frame(seurat.concat@meta.data)

# Add UMAP coordinates to the metadata
dat$UMAP1 <- seurat.concat@reductions$umap@cell.embeddings[,1]
dat$UMAP2 <- seurat.concat@reductions$umap@cell.embeddings[,2]


# Plot 2D density of UMAP embeddings, split by sample
ggplot(dat, aes(x = UMAP1, y = UMAP2)) +
  geom_density_2d_filled(contour_var = "ndensity") +  # Filled contour plot
  geom_density_2d(linewidth = 0.25, colour = "black") +  # Line contour
  facet_wrap(~ SampleName, nrow = 3, switch = "y") +  # Split by unique sample ID
  scale_x_continuous(expand = c(0, 0)) +  # Remove padding from x-axis
  scale_y_continuous(expand = c(0, 0))  # Remove padding from y-axis

ggsave("density_UMAP-ID.pdf", width = 30, height = 15, limitsize = FALSE)


ggplot(dat, aes(x = UMAP1, y = UMAP2)) +
  geom_density_2d_filled(contour_var = "ndensity") +  # Filled contour plot
  geom_density_2d(linewidth = 0.25, colour = "black") +  # Line contour
  facet_wrap(~ Immunecomp, nrow = 3, switch = "y") +  # Split by unique sample ID
  scale_x_continuous(expand = c(0, 0)) +  # Remove padding from x-axis
  scale_y_continuous(expand = c(0, 0))  # Remove padding from y-axis

ggsave("density_UMAP-ImmuneComp.pdf", width = 30, height = 30, limitsize = FALSE)


ggplot(dat, aes(x = UMAP1, y = UMAP2)) +
  geom_density_2d_filled(contour_var = "ndensity") +  # Filled contour plot
  geom_density_2d(linewidth = 0.25, colour = "black") +  # Line contour
  facet_wrap(~ ALC, nrow = 3, switch = "y") +  # Split by unique sample ID
  scale_x_continuous(expand = c(0, 0)) +  # Remove padding from x-axis
  scale_y_continuous(expand = c(0, 0))  # Remove padding from y-axis

ggsave("density_UMAP-ALC.pdf", width = 30, height = 30, limitsize = FALSE)

dat_sampled <- dat %>%
  group_by(SampleName) %>%
  slice_sample(n = 5000) %>%
  ungroup()
table(dat_sampled$SampleName)


ggplot(dat_sampled, aes(x = UMAP1, y = UMAP2)) +
  # Blue points in the background
  geom_point(data = subset(dat_sampled, guzmano_vs_hiALC_vs_loALC == "Lower"), 
             aes(color = guzmano_vs_hiALC_vs_loALC), size = 1, alpha = 0.3, fill = "seagreen") + 
  # Orange points in the middleground
  geom_point(data = subset(dat_sampled, guzmano_vs_hiALC_vs_loALC == "Higher"), 
             aes(color = guzmano_vs_hiALC_vs_loALC), size = 1,shape = 21, fill = "#F6B26B", alpha = 0.3) +
  # Red points on top with larger size and striking color
  geom_point(data = subset(dat_sampled, guzmano_vs_hiALC_vs_loALC == "251"), 
             aes(color = guzmano_vs_hiALC_vs_loALC), size = 1, shape = 21, fill = "#7F3F98", alpha = 1) +
  scale_color_manual(values = c("Lower" = "seagreen","Higher"= "#F6B26B", "251" = "#7F3F98")) +
  theme_minimal() +
  labs(title = "UMAP 251 vs hiALC vs loALC", color = "ALC")

ggsave("UMAP-251vsHighALCvsLowALC.pdf", 
       width = 30, 
       height = 20, 
       units = "cm",  # Specify the units (e.g., cm, inches, etc.)
       device = "pdf", 
       dpi = 300,  # High resolution
       limitsize = FALSE)


# Step 5: Save and Load Intermediate Files ----
# Save Seurat object to an intermediate file for future use
setwd(int_base)  # Change working directory to intermediate files directory
saveRDS(seurat.concat, file = "seurat_object.rds")  # Save Seurat object

# Load Seurat object from an intermediate file if needed
seurat.concat <- readRDS(file = "seurat_object.rds")




