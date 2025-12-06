import argparse
import time
import pandas as pd
import numpy as np
import scanpy as sc
import torch
import os

# Import your modified GraphST
from GraphST import GraphST

def parse_args():
    parser = argparse.ArgumentParser(description="Run GraphST for Spotless Benchmark")
    parser.add_argument("sc_input", type=str, help="Path to single-cell H5AD file")
    parser.add_argument("sp_input", type=str, help="Path to spatial H5AD file")
    parser.add_argument("annot", type=str, help="Column name for cell type annotations in sc_input")
    parser.add_argument("--datatype", type=str, default="10X", 
                        help="GraphST datatype/model enum (e.g., '10X', 'GAT', 'SGC', 'GATv2')")
    parser.add_argument("--device", type=str, default="auto", help="Compute device (cpu or cuda)")
    parser.add_argument("--epochs", type=int, default=600, help="Number of training epochs")
    parser.add_argument("--output", type=str, default="proportions.tsv", help="Output file path")
    return parser.parse_args()

def main():
    args = parse_args()

    # 1. Device Setup
    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)
    print(f"Using device: {device}")

    # 2. Load Data
    print(f"Loading single-cell data from {args.sc_input}...")
    adata_sc = sc.read_h5ad(args.sc_input)
    
    print(f"Loading spatial data from {args.sp_input}...")
    adata_sp = sc.read_h5ad(args.sp_input)

    # Ensure annotation column exists
    if args.annot not in adata_sc.obs.columns:
        raise ValueError(f"Annotation column '{args.annot}' not found in single-cell data.")

    # 3. Preprocessing (Standard Scanpy/GraphST checks)
    # GraphST's built-in preprocess expects 'highly_variable' genes or raw counts.
    # Spotless data is usually already preprocessed, but we ensure basic consistency.
    adata_sc.var_names_make_unique()
    adata_sp.var_names_make_unique()

    # Intersection of genes (Crucial for GraphST)
    intersect = np.intersect1d(adata_sc.var_names, adata_sp.var_names)
    adata_sc = adata_sc[:, intersect].copy()
    adata_sp = adata_sp[:, intersect].copy()
    print(f"Data intersection: {len(intersect)} common genes.")

    # 4. Initialize and Train GraphST
    # We pass the 'datatype' argument to trigger your custom GAT/SGC logic
    model = GraphST(
        adata=adata_sp,
        adata_sc=adata_sc,
        device=device,
        epochs=args.epochs,
        datatype=args.datatype, # This triggers your modified code logic
        deconvolution=True      # Essential for mapping
    )

    print(f"Training GraphST model (Type: {args.datatype})...")
    start_time = time.time()
    
    # Train and learn mapping
    adata_sp, adata_sc = model.train_map()
    
    end_time = time.time()
    print(f"Training finished in {end_time - start_time:.2f} seconds.")

    # 5. Extract and Aggregate Results
    # GraphST outputs a (Spot x Cell) matrix in adata.obsm['map_matrix']
    # We need to sum these by Cell Type to get (Spot x CellType) proportions.
    
    map_matrix = adata_sp.obsm['map_matrix'] # Shape: (n_spots, n_cells)
    
    # Get cell types corresponding to the columns (cells) of map_matrix
    # Note: GraphST aligns adata_sc internally, but usually preserves order. 
    # To be safe, we rely on the input adata_sc observations.
    sc_labels = adata_sc.obs[args.annot].values
    
    # Create DataFrame for aggregation
    df_map = pd.DataFrame(map_matrix, index=adata_sp.obs_names)
    
    # Aggregate columns by cell type
    # We create a dummy dataframe to map cell_index -> cell_type
    print("Aggregating cell-to-spot mapping into cell type proportions...")
    unique_labels = np.unique(sc_labels)
    proportions = pd.DataFrame(index=adata_sp.obs_names, columns=unique_labels)
    
    # Efficient summation per category
    for label in unique_labels:
        # Find indices of cells belonging to this type
        indices = np.where(sc_labels == label)[0]
        if len(indices) > 0:
            # Sum the probabilities for these cells
            proportions[label] = map_matrix[:, indices].sum(axis=1)
        else:
            proportions[label] = 0.0

    # 6. Save Output
    # Normalize rows to sum to 1 (optional but recommended for proportions)
    proportions = proportions.div(proportions.sum(axis=1), axis=0).fillna(0)
    
    print(f"Saving results to {args.output}...")
    proportions.to_csv(args.output, sep="\t", index=True, header=True)

if __name__ == "__main__":
    main()