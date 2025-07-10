# BMCA_CAR_irAE_NatMed
BMCA-CAR_irAE_NatMed

Code for unbiased analysis of flow cytometry data 

1. Read and Transform Data

Import FCS files from FlowJo.
Apply an asinh transformation (adjust cofactor as needed).
Convert to Seurat for Clustering & UMAP

2. Construct a Seurat object from transformed data.
   
Normalize, scale, and perform PCA.
Run UMAP for visualization.
FlowSOM Clustering & Cluster Analysis

3. Apply FlowSOM for clustering.
   
Compare cluster abundances across conditions.


Code for scRNAseq analysis 

**1.0_integrate_all_samples** 
Integrates all filtered single-cell samples into a unified Seurat object using standard workflows.

**2.0_filter_ciltacelpos**
Filters the dataset to retain only CAR-positive (Cilta-cel) T cells based on transgene expression.

**3.0_adjust_clusters** 
Performs clustering and UMAP dimensionality reduction, followed by manual annotation of clusters.

**4.0_Generate_NebulosaPlots** 
Generates Nebulosa density plots to visualize expression of selected genes across clusters.

**5.0_Calculating_abundances** 
Computes cell-type frequencies by cluster and compartment (PBMC vs CSF), and visualizes compositions.

**6.0_AUCell_Scoring** 
Scores gene sets per cell using AUCell and compares pathway activity across clusters and compartments.

**7.0_Cluster_Identification** 
Identifies marker genes per cluster and visualizes expression of functional gene modules via DotPlot.

