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
