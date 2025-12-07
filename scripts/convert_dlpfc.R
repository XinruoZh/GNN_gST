# File: scripts/convert_dlpfc.R
library(zellkonverter)
library(SingleCellExperiment)
library(HDF5Array) # Required for this data format

# Set working directory to where the data is
setwd("/ocean/projects/cis250160p/xzhaoa/GNN_gST/data/DLPFC")

# 1. Load the HDF5-backed experiment (It's a directory, not a file)
print("Loading HDF5 dataset from directory 'sce_DLPFC_annotated'...")
sce <- loadHDF5SummarizedExperiment("sce_DLPFC_annotated")

print("Object loaded:")
print(sce)

# 2. Save as .h5ad
print("Converting to H5AD (this may take a minute)...")
writeH5AD(sce, "DLPFC_snRNAseq_ref.h5ad")

print("Success! Saved as DLPFC_snRNAseq_ref.h5ad")