# Load required packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Set paths for input/output
base <- "/scRNAseq_analysis/"  # base directory for the project (adjust as needed)
objSaveDir <- file.path(base, "seurat_objects")  # where Seurat RDS files will be saved
fig_base <- file.path(base, "figures")  # where QC and UMAP plots will go
data.dir <- "data/cellranger"  # relative path to Cell Ranger outputs

# Create folders if they don't already exist
dir.create(objSaveDir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_base, "individual_samples"), recursive = TRUE, showWarnings = FALSE)

# Get sample folder names (expects Cell Ranger output folders to start with "HTS")
samplenames <- list.files(data.dir, pattern = "^HTS")
seurat.objs <- list()

# Loop over each sample and create a Seurat object
for (current.samplename in samplenames) {
  message("Processing: ", current.samplename)
  
  # Construct full path to Cell Ranger output matrix
  current.filename <- paste0(data.dir, "/", current.samplename, "/outs/filtered_feature_bc_matrix.h5")
  import.obj <- Read10X_h5(current.filename)
  
  # Create Seurat object with basic filtering
  seurat.obj <- CreateSeuratObject(
    counts = import.obj,
    min.cells = 3,
    min.features = 200,
    project = current.samplename
  )
  
  # Calculate mitochondrial % per cell
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")
  
  # Save violin plots before filtering
  p <- VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
  ggsave(paste0(fig_base, "/individual_samples/", current.samplename, "_qc_vlnplot.pdf"), p, width=6, height=6)
  
  # Save histograms before filtering
  pdf(paste0(fig_base, "/individual_samples/", current.samplename, "_unfiltered_qc_histograms.pdf"), width=6, height=8)
  par(mfrow=c(3,1))
  hist(seurat.obj$nFeature_RNA, breaks=200, main = "nFeature_RNA")
  hist(seurat.obj$nCount_RNA, breaks=200, main = "nCount_RNA")
  hist(seurat.obj$percent.mt, breaks=200, main = "percent.mt")
  dev.off()
  par(mfrow=c(1,1))
  
  # Filter out low-quality cells
  seurat.obj <- subset(seurat.obj, subset = percent.mt < 10 & nFeature_RNA > 200 & nFeature_RNA < 10000)
  
  # Save violin plots after filtering
  p <- VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
  ggsave(paste0(fig_base, "/individual_samples/", current.samplename, "_filtered_qc_vlnplot.pdf"), p, width=6, height=6)
  
  # Save histograms after filtering
  pdf(paste0(fig_base, "/individual_samples/", current.samplename, "_filtered_qc_histograms.pdf"), width=6, height=8)
  par(mfrow=c(3,1))
  hist(seurat.obj$nFeature_RNA, breaks=200, main = "nFeature_RNA")
  hist(seurat.obj$nCount_RNA, breaks=200, main = "nCount_RNA")
  hist(seurat.obj$percent.mt, breaks=200, main = "percent.mt")
  dev.off()
  
  # Standard Seurat preprocessing
  seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
  seurat.obj <- ScaleData(seurat.obj)
  seurat.obj <- RunPCA(seurat.obj)
  seurat.obj <- FindNeighbors(seurat.obj, dims=1:15)
  seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)
  seurat.obj <- RunUMAP(seurat.obj, dims=1:15)
  
  # Save UMAP plots and marker expression
  ggsave(paste0(fig_base, "/individual_samples/", current.samplename, "_initial_UMAP_clusters.pdf"),
         DimPlot(seurat.obj), width=5.5, height=5)

  ggsave(paste0(fig_base, "/individual_samples/", current.samplename, "_initial_UMAP_Tcell-cellcycle-markers.pdf"),
         FeaturePlot(seurat.obj, features=c("CD3E", "CD14", "CD4", "CD8A", "MKI67", "TOP2A", "TCF7", "TBX21", "ciltacel"), ncol=3),
         width=12, height=7)

  ggsave(paste0(fig_base, "/individual_samples/", current.samplename, "_initial_UMAP_Tcell-qc-metrics.pdf"),
         FeaturePlot(seurat.obj, features=c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol=3),
         width=10, height=3)
  
  # Store object in list
  seurat.objs[[current.samplename]] <- seurat.obj
}

# Save individual sample objects (not integrated)
saveRDS(seurat.objs, file = file.path(objSaveDir, "seurat.objs_nonintegrated.RDS"))

# Integration workflow across all samples
integration.features <- SelectIntegrationFeatures(object.list = seurat.objs)

start.time <- Sys.time()
anchors <- FindIntegrationAnchors(object.list = seurat.objs, anchor.features = integration.features)
end.time <- Sys.time()
print(end.time - start.time)

saveRDS(anchors, file = file.path(objSaveDir, "anchors_integrate_all_samples.RDS"))

# Free up memory before integrating
rm(seurat.objs)

# Run data integration
Tcells.allsamples <- IntegrateData(anchorset = anchors, dims = 1:50)
DefaultAssay(Tcells.allsamples) <- "integrated"

# Run downstream clustering and dimensionality reduction
Tcells.allsamples <- ScaleData(Tcells.allsamples, verbose = FALSE)
Tcells.allsamples <- RunPCA(Tcells.allsamples, npcs = 15, verbose = FALSE)
Tcells.allsamples <- RunUMAP(Tcells.allsamples, reduction = "pca", dims = 1:15)
Tcells.allsamples <- FindNeighbors(Tcells.allsamples, reduction = "pca", dims = 1:15, k.param = 15)
Tcells.allsamples <- FindClusters(Tcells.allsamples, resolution = 0.5)

# Plot UMAP + marker gene expression (these are not saved to file by default)
DimPlot(Tcells.allsamples, label = TRUE)
FeaturePlot(Tcells.allsamples, features = c("CD3E", "CD14", "CD4", "CD8A", "MKI67", "TOP2A", "TCF7", "TBX21", "ciltacel"), ncol = 3)
FeaturePlot(Tcells.allsamples, features = c("rna_CD8A", "rna_CD4", "rna_CAR19", "rna_FOXP3"))
FeaturePlot(Tcells.allsamples, features = c("rna_GZMK", "rna_GZMB", "rna_PRF1", "rna_GZMA"))
FeaturePlot(Tcells.allsamples, features = c("rna_IL7R", "rna_CCR7", "rna_TCF7", "rna_CCL5"))
FeaturePlot(Tcells.allsamples, features = c("rna_PCNA", "rna_MKI67", "rna_TOP2A", "rna_CCNB1"))
FeaturePlot(Tcells.allsamples, features = c("rna_IRF7", "rna_RSAD2", "rna_ISG15", "rna_OAS1"))

# Show violin plots for key features
VlnPlot(Tcells.allsamples, features = c("rna_CD3E", "rna_CD4", "rna_CD8", "rna_CD14", "rna_CD19", "ciltacel", "nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)

# Save final integrated object
saveRDS(Tcells.allsamples, file = file.path(objSaveDir, "Tcells.allsamples_initial-integration.RDS"))
