# Load required libraries ----
library(flowCore)        # Core package for handling flow cytometry data
library(flowVS)          # Tools for analyzing and visualizing flow cytometry data
library(flowViz)         # Flow cytometry visualization package
library(flowWorkspace)   # Flow cytometry data analysis and workspace management
library(Seurat)          # Single-cell data analysis toolkit
library(Matrix)          # Sparse and dense matrix operations
library(ggplot2)         # Data visualization package
library(tidyverse)       # Collection of data manipulation and visualization tools
library(dplyr)           # Data manipulation functions
library(gtools)          # Assorted R programming tools
library(wadeTools)       # Specific tools for flow cytometry data
library(gridExtra)       # Functions to arrange multiple plots
library(viridis)         # Color scales for scientific data visualization

# Step 0: Define work directories ----
# Set up project base and directory paths for input and output files
proj_base = ""

# Define directories for data, pictures, gated FCS files, and intermediate subjects
data_base = tight(proj_base, "/FlowJo_Excel/CD3")  # FlowJo CSV export directory
fig_base = tight(proj_base, "/figures/CD3")       # Directory for saving pictures
int_base = tight(proj_base,"/intermediate_subjects/CD3")  # Directory for intermediate files

setwd(data_base)  # Set working directory to FlowJo CSV export directory

# Step 1: Read in FlowJo exported CSV with CD4 and CD8 T cells ----
# Load CSV files containing single-cell MFI data for CD4 and CD8 T cells
file_paths <- list.files(path = data_base, pattern = "*.csv", full.names = TRUE)  # List all CSV files in the directory

data_list <- lapply(file_paths, read.csv)  # Read each CSV file into a list of data frames

# Remove unwanted columns related to forward scatter (FSC) and side scatter (SSC)
columns_to_remove <- c("FSC.A", "FSC.H", "SSC.A", "SSC.B.A", "SSC.B.H", "SSC.H", "Time", "Comp.AF.A","Comp.LIVE.DEAD.Blue.A....LD.Blue", "Comp.BV750.A....CD3" ,"Comp.BV480.A....CD19","Comp.BUV563.A....CD14","Comp.BUV496.A....CD16")
data_list <- lapply(data_list, function(df) df[, !colnames(df) %in% columns_to_remove])  # Remove specific columns from each data frame

# Extract sample names from file paths using regex and assign names to the data list
sample_names <- sub(".*/(\\d{3})\\.csv$", "\\1", file_paths)  # Extract sample names based on pattern
names(data_list) <- sample_names  # Name the data list based on sample names

# Step 2: Filter and concatenate data ----
# Filter data to include a maximum of 100,000 events per sample
data_list_filtered <- lapply(data_list, function(ff) {
  set.seed(0)  # Set seed for reproducibility
  if(nrow(ff) > 100000) {  # If more than 100,000 events
    ff <- ff[sample(1:nrow(ff), 100000),]  # Randomly select 100,000 events
  }
  return(ff)  # Return filtered data
})

# Concatenate all filtered data into a single matrix
common_colnames <- colnames(data_list_filtered[[1]])
data_list_filtered <- lapply(data_list_filtered, function(df) {
  colnames(df) <- common_colnames
  df
})

flowmatrix.concat <- Reduce(rbind, data_list_filtered)  # Combine data frames by row
flowmatrix.concat.t <- t(flowmatrix.concat)  # Transpose the matrix

# Create unique cell identifiers for the concatenated data
Cells <- unlist(lapply(names(data_list_filtered), function(x) {
  paste0("Cell_", x, "_", 1:nrow(data_list_filtered[[x]]))  # Generate unique cell names
}))

# Assign cell names to columns of the transposed matrix
colnames(flowmatrix.concat.t) <- Cells  # Set cell names as column names
setwd(proj_base)
ic_table <- read.csv("JHN081_metadata.csv")

# Create metadata with sample names, immune complications, ALC
SampleName <- rep(names(data_list_filtered), times = sapply(data_list_filtered, nrow))  # Repeat sample names for each event
Immunecomp <- rep(ic_table$Immune.complications, times = sapply(data_list_filtered, nrow))  # Repeat sample names for each event
ALC <- rep(ic_table$ALC.Expansion.CO.0.78, times = sapply(data_list_filtered, nrow))

metadata <- cbind(Cells, SampleName, Immunecomp, ALC)  # Combine cell names and sample names into metadata
metadata <- data.frame(metadata)  # Convert metadata to a data frame
rownames(metadata) <- metadata$Cells  # Set cell names as row names in the metadata

# Step 3: Transformation ----
# Apply asinh transformation to the flow cytometry data using per-marker cofactors
mean_values <- apply(flowmatrix.concat, 2, mean, na.rm = TRUE)  # Calculate mean for each marker
cofactors <- mean_values / 10  # Calculate cofactors for asinh transformation

# Define asinh transformation function
asinh_transform <- function(x, cofactor) {
  return(asinh(x / cofactor))  # Apply asinh transformation
}

# Apply asinh transformation to each marker in the data
transformed_data <- mapply(asinh_transform, 
                           as.data.frame(flowmatrix.concat),  # Data frame for transformation
                           cofactors,  # Cofactors for each marker
                           SIMPLIFY = FALSE)  # Simplify result to a list

# Concatenate transformed data and transpose it for consistency
transformed_data <- do.call(cbind, transformed_data)  # Combine the list into a data frame
transformed_data.t <- t(transformed_data)  # Transpose transformed data
colnames(transformed_data.t) <- Cells  # Assign cell names to columns


# Step 4: Visualize original and transformed data ----
# Compare original and transformed data using histograms
par(mfrow = c(2, 1))  # Arrange plots in a 2-row layout
hist(flowmatrix.concat[, 12], main = "Original Data", xlab = "Value", col = "lightblue")  # Histogram for original data
hist(transformed_data[, 12], main = "Transformed Data", xlab = "Value", col = "lightgreen")  # Histogram for transformed data
title(colnames(transformed_data)[12])  # Add title based on the marker name

# Step 5: Save histograms for each marker ----
# Check if the "histograms" directory exists and create it if not
if (!dir.exists("histograms")) {
  dir.create("histograms")  # Create directory
}

# Get marker names
marker_names <- colnames(transformed_data)  # Extract marker names from transformed data

# Loop through each marker and save histograms for original and transformed data
for (i in 1:ncol(flowmatrix.concat)) {
  marker_name <- marker_names[i]  # Get marker name
  
  # Original Data Histogram
  png(filename = paste0("histograms/Original_", marker_name, ".png"))  # Define file path
  hist(flowmatrix.concat[, i], main = paste("Original Data -", marker_name),  # Histogram for original data
       xlab = "Value", col = "lightblue", breaks = 50)  # Set histogram parameters
  dev.off()  # Close the file
  
  # Transformed Data Histogram
  png(filename = paste0("histograms/Transformed_", marker_name, ".png"))  # Define file path
  hist(transformed_data[, i], main = paste("Transformed Data -", marker_name),  # Histogram for transformed data
       xlab = "Value", col = "lightgreen", breaks = 50)  # Set histogram parameters
  dev.off()  # Close the file
}


# Step 6: Save results ----
# Save the transposed flow matrix and transformed matrix to file
setwd(int_base)
saveRDS(flowmatrix.concat.t, file = "flowMatrix.rds")  # Save original data matrix
saveRDS(transformed_data.t, file = "flowMatrix_transformed.rds")  # Save transformed data matrix
saveRDS(metadata, file = "metadata.rds")  # Save metadata
