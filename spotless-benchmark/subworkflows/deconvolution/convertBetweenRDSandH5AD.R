#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(SeuratDisk)
library(Seurat)
library(tools)

par <- list(output_file_ext = "")
args <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)
par[names(args)] <- args

file_name <- tools::file_path_sans_ext(basename(par$input_path))
ext <- tools::file_ext(par$input_path)

if (tolower(ext) == "rds"){
  cat(paste0(">>> Processing: ", par$input_path, "\n"))
  input_obj <- readRDS(par$input_path)
  
  if (!inherits(input_obj, "Seurat")){
    
    # 1. Extract Counts
    counts <- NULL
    if ("counts" %in% names(input_obj)) counts <- input_obj$counts
    else if ("raw_counts" %in% names(input_obj)) counts <- input_obj$raw_counts
    
    if (is.null(counts)) stop("CRITICAL ERROR: No 'counts' found in list.")
    
    seurat_obj <- CreateSeuratObject(counts = counts)
    
    # 2. Extract Coordinates
    coords <- NULL
    possible_names <- c("coordinates", "spatial", "locs", "locations", "centroids", "offsets")
    
    for (name in possible_names) {
        if (name %in% names(input_obj)) {
            cat(paste0(">>> Found coordinates in slot: '", name, "'\n"))
            coords <- as.data.frame(input_obj[[name]])
            break
        }
    }
    
    # --- FALLBACK: Generate Grid if Missing ---
    if (is.null(coords)) {
        cat("WARNING: No spatial coordinates found in file. Generating synthetic grid...\n")
        
        n_spots <- ncol(seurat_obj)
        # Calculate grid dimensions (approx square)
        rows <- ceiling(sqrt(n_spots))
        cols <- ceiling(n_spots / rows)
        
        # Generate x, y coordinates
        x <- rep(1:cols, times = rows)[1:n_spots]
        y <- rep(1:rows, each = cols)[1:n_spots]
        
        coords <- data.frame(x = x, y = y)
        rownames(coords) <- colnames(seurat_obj)
        cat(paste0(">>> Created dummy grid: ", rows, "x", cols, "\n"))
    }
    
    # Standardize Column Names (ensure x, y)
    colnames(coords) <- tolower(colnames(coords))
    if (!all(c("x", "y") %in% colnames(coords))) {
         if (ncol(coords) >= 2 && is.numeric(coords[,1]) && is.numeric(coords[,2])) {
             colnames(coords)[1:2] <- c("x", "y")
         }
    }

    # Add Metadata
    if (nrow(coords) == ncol(seurat_obj)) {
        seurat_obj <- AddMetaData(seurat_obj, metadata = coords)
        cat(">>> Success: Coordinates added to metadata.\n")
    } else {
        stop("CRITICAL ERROR: Coordinate dimensions do not match count matrix.")
    }

    # Add generic metadata
    if ("meta.data" %in% names(input_obj)) {
         seurat_obj <- AddMetaData(seurat_obj, metadata = as.data.frame(input_obj$meta.data))
    }

  } else {
    seurat_obj <- input_obj
  }
  
  rm(input_obj)
  DefaultAssay(seurat_obj) <- names(seurat_obj@assays)[grep("RNA|Spatial", names(seurat_obj@assays))[1]]
  
  if (file.exists(paste0(file_name, ".h5ad"))) file.remove(paste0(file_name, ".h5ad"))
  
  cat(">>> Writing .h5Seurat...\n")
  SaveH5Seurat(seurat_obj, filename = paste0(file_name, ".h5seurat"), overwrite = TRUE)
  
  cat(">>> Converting to .h5ad...\n")
  Convert(paste0(file_name, ".h5seurat"), dest = "h5ad", overwrite = TRUE)
  file.remove(paste0(file_name, ".h5seurat"))

} else if (tolower(ext) == "h5ad"){
  h5seurat_path <- file.path(dirname(Sys.readlink(par$input_path)), paste0(file_name, ".h5seurat")) 
  if (!file.exists(h5seurat_path)) Convert(par$input_path, dest = "h5seurat")
  saveRDS(LoadH5Seurat(h5seurat_path), paste0(file_name, par$output_file_ext, ".rds"))
  file.remove(h5seurat_path)
}
