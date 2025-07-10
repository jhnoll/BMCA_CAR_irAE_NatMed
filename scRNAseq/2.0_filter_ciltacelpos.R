library(Seurat)

# Load integrated object containing all cells, and list of individual sample Seurat objects
Tcells.allsamples <- readRDS(file = "seurat_objs/Tcells.allsamples_initial-integration.RDS")
seurat.objs <- readRDS(file = "seurat_objs/seurat.objs_nonintegrated.RDS")

# ---------------------------
# Step 1: Visualize all cells
# ---------------------------

# UMAP plot colored by cluster
pdf("figures/integration_allcells/UMAP_allcells_bycluster.pdf", width = 8, height = 8)
DimPlot(Tcells.allsamples, label = TRUE)
dev.off()

# UMAP plots of key marker genes to evaluate cluster identity
pdf("figures/integration_allcells/UMAP_allcells_markers.pdf", width = 20, height = 20)
FeaturePlot(Tcells.allsamples, features = c("CD3E", "CD14", "CD4", "CD8A", "MKI67", "TOP2A", "TCF7", "TBX21", "ciltacel"), ncol = 3)
dev.off()

# Violin plots of QC metrics and lineage markers across clusters
pdf("figures/integration_allcells/VlnPlot_allcells_bycluster.pdf", width = 16, height = 16)
VlnPlot(Tcells.allsamples, features = c("rna_CD3E", "rna_CD4", "rna_CD8A", "rna_CD14", "rna_CD19", "ciltacel",
                                        "nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
dev.off()

# -------------------------------
# Step 2: Define filtering criteria
# -------------------------------

# Dimensions of full dataset
dim(Tcells.allsamples)

# Cells per sample
table(Tcells.allsamples$orig.ident)

# Cells per cluster × sample
table(Tcells.allsamples$orig.ident, Tcells.allsamples$seurat_clusters)

# Manually define cluster IDs to retain (CD3+ T cells)
# Removed: 3, 5, 8, 12, 13, 15 due to low CD3E, CD14+, or low-quality
# Retained: 0, 1, 2, 4, 6, 7, 9, 10, 12, 14
clusters.keep <- c("0", "1", "2", "4", "6", "7", "9", "10", "12", "14")

# -----------------------------------
# Step 3: Identify ciltacel+ cells
# -----------------------------------

# Extract cells with non-zero ciltacel transgene from each RNA matrix
cells.with.ciltacel <- c(
  colnames(Tcells.allsamples@assays$RNA$counts.HTS_SO188_01_251_CSF)[Tcells.allsamples@assays$RNA$counts.HTS_SO188_01_251_CSF["ciltacel", ] > 0],
  colnames(Tcells.allsamples@assays$RNA$counts.HTS_SO188_02_439_CSF)[Tcells.allsamples@assays$RNA$counts.HTS_SO188_02_439_CSF["ciltacel", ] > 0],
  colnames(Tcells.allsamples@assays$RNA$counts.HTS_SO188_03_470_CSF)[Tcells.allsamples@assays$RNA$counts.HTS_SO188_03_470_CSF["ciltacel", ] > 0],
  colnames(Tcells.allsamples@assays$RNA$counts.HTS_SO188_04_492_CSF)[Tcells.allsamples@assays$RNA$counts.HTS_SO188_04_492_CSF["ciltacel", ] > 0],
  colnames(Tcells.allsamples@assays$RNA$counts.HTS_SO188_05_251_PBMC)[Tcells.allsamples@assays$RNA$counts.HTS_SO188_05_251_PBMC["ciltacel", ] > 0],
  colnames(Tcells.allsamples@assays$RNA$counts.HTS_SO188_06_439_PBMC)[Tcells.allsamples@assays$RNA$counts.HTS_SO188_06_439_PBMC["ciltacel", ] > 0],
  colnames(Tcells.allsamples@assays$RNA$counts.HTS_SO188_07_470_PBMC)[Tcells.allsamples@assays$RNA$counts.HTS_SO188_07_470_PBMC["ciltacel", ] > 0],
  colnames(Tcells.allsamples@assays$RNA$counts.HTS_SO188_08_492_PBMC)[Tcells.allsamples@assays$RNA$counts.HTS_SO188_08_492_PBMC["ciltacel", ] > 0]
)

# Keep only cells that are both in T cell clusters and ciltacel+
cells.in.clusters.keep <- names(Tcells.allsamples$seurat_clusters)[Tcells.allsamples$seurat_clusters %in% clusters.keep]
cells.keep <- intersect(cells.with.ciltacel, cells.in.clusters.keep)

# Clean up memory
rm(Tcells.allsamples)

# -----------------------------------
# Step 4: Subset original sample objects to cells.keep
# -----------------------------------

for (i in 1:8) {
  samplename <- names(seurat.objs)[i]
  current.cells.keep <- grep(paste0("_", i), cells.keep, value = TRUE)
  current.cells.keep <- sub(paste0("_", i), "", current.cells.keep)
  seurat.objs[[samplename]] <- subset(seurat.objs[[samplename]], cells = sub("_.*", "", cells.keep))
}

# Check retained cell counts
sapply(seurat.objs, ncol)
sum(sapply(seurat.objs, ncol))  # Should be ~43,694

# -----------------------------------
# Step 5: Re-integrate only ciltacel+ T cells
# -----------------------------------

# Select integration features and find anchors
integration.features <- SelectIntegrationFeatures(object.list = seurat.objs)
start.time <- Sys.time()
anchors <- FindIntegrationAnchors(object.list = seurat.objs, anchor.features = integration.features)
end.time <- Sys.time()
end.time - start.time

saveRDS(anchors, file = "seurat_objs/anchors_integrate_filtered_ciltacel.RDS")

# Clean up memory
rm(seurat.objs)

# Handle integration issue for small samples: k.weight must be ≤ min cells/sample
# Sample 3 has only 60 cells, so we set k.weight = 60
Tcells.filtered.ciltacel <- IntegrateData(anchorset = anchors, dims = 1:50, k.weight = 60)

DefaultAssay(Tcells.filtered.ciltacel) <- "integrated"

# Standard post-integration workflow
Tcells.filtered.ciltacel <- ScaleData(Tcells.filtered.ciltacel, verbose = FALSE)
Tcells.filtered.ciltacel <- RunPCA(Tcells.filtered.ciltacel, npcs = 15, verbose = FALSE)
Tcells.filtered.ciltacel <- RunUMAP(Tcells.filtered.ciltacel, reduction = "pca", dims = 1:15)
Tcells.filtered.ciltacel <- FindNeighbors(Tcells.filtered.ciltacel, reduction = "pca", dims = 1:15, k.param = 15)
Tcells.filtered.ciltacel <- FindClusters(Tcells.filtered.ciltacel, resolution = 0.5)

# -----------------------------------
# Step 6: Visualization
# -----------------------------------

DimPlot(Tcells.filtered.ciltacel, label = TRUE)

FeaturePlot(Tcells.filtered.ciltacel, features = c("CD3E", "CD14", "CD4", "CD8A", "MKI67", "TOP2A", "TCF7", "TBX21", "ciltacel"), ncol = 3)
FeaturePlot(Tcells.filtered.ciltacel, features = c("rna_CD8A", "rna_CD4", "rna_cltacel", "rna_FOXP3"))
FeaturePlot(Tcells.filtered.ciltacel, features = c("rna_GZMK", "rna_GZMB", "rna_PRF1", "rna_GZMA"))
FeaturePlot(Tcells.filtered.ciltacel, features = c("rna_IL7R", "rna_CCR7", "rna_TCF7", "rna_CCL5"))
FeaturePlot(Tcells.filtered.ciltacel, features = c("rna_PCNA", "rna_MKI67", "rna_TOP2A", "rna_CCNB1"))
FeaturePlot(Tcells.filtered.ciltacel, features = c("rna_IRF7", "rna_RSAD2", "rna_ISG15", "rna_OAS1"))

# Violin plots by cluster and sample
VlnPlot(Tcells.filtered.ciltacel, features = c("rna_CD3E", "rna_CD4", "rna_CD8", "rna_CD14", "rna_CD19", "ciltacel", "nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)

VlnPlot(Tcells.filtered.ciltacel, features = c("rna_CD4", "rna_CD8", "rna_GZMA", "rna_GZMK", "GZMB"), group.by = "orig.ident", pt.size = 0)

# Save final filtered integrated object
saveRDS(Tcells.filtered.ciltacel, file = "seurat_objs/Tcells.filtered.ciltacel_integration.RDS")
