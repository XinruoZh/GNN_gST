# File: get_dlpfc.R

# conda install -c bioconda bioconductor-spatiallibd bioconductor-zellkonverter
# Load the libraries (installed via Conda)
library(spatialLIBD)
library(zellkonverter)

# Set working directory
setwd("/ocean/projects/cis250160p/xzhaoa/GNN_gST/data/DLPFC")

print("Downloading all 12 DLPFC samples... (This may take a few minutes)")

# 1. Fetch data (uses the cache in your home dir by default)
spe <- fetch_data(type = "spe")

print("Download complete. Converting to AnnData...")

# 2. Save as .h5ad for Python
# The output file will be approx 600MB - 1GB
writeH5AD(spe, "DLPFC_all_samples.h5ad")

print("Success! File saved as 'DLPFC_all_samples.h5ad'")