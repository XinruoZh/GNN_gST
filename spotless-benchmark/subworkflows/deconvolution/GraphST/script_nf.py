import argparse
import time
import pandas as pd
import numpy as np
import scanpy as sc
import torch
import os
from scipy.sparse import issparse

# --- Import Fixes ---
from GraphST.GraphST import GraphST
from GraphST.preprocess import construct_interaction

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

    # --- DATA FIX: Ensure 'spatial' key exists in obsm ---
    if 'spatial' not in adata_sp.obsm:
        print("WARNING: 'spatial' key missing in adata.obsm. Attempting to fix...")
        found = False
        
        # Strategy A: Check alternative obsm keys
        for key in ['X_spatial', 'spatial_coordinates', 'coordinates']:
            if key in adata_sp.obsm:
                adata_sp.obsm['spatial'] = adata_sp.obsm[key]
                print(f" -> Found coordinates in obsm['{key}']. Copied to obsm['spatial'].")
                found = True
                break
        
        # Strategy B: Check obs columns (x/y)
        if not found:
            pairs = [('x', 'y'), ('X', 'Y'), ('x_centroid', 'y_centroid'), ('pxl_col_in_fullres', 'pxl_row_in_fullres')]
            for x_col, y_col in pairs:
                obs_lower = {k.lower(): k for k in adata_sp.obs.columns}
                if x_col.lower() in obs_lower and y_col.lower() in obs_lower:
                    real_x = obs_lower[x_col.lower()]
                    real_y = obs_lower[y_col.lower()]
                    adata_sp.obsm['spatial'] = adata_sp.obs[[real_x, real_y]].values
                    print(f" -> Found coordinates in obs columns '{real_x}' and '{real_y}'. Created obsm['spatial'].")
                    found = True
                    break
                    
        if not found:
             print("Available obsm keys:", list(adata_sp.obsm.keys()))
             print("Available obs columns:", list(adata_sp.obs.columns))
             raise KeyError("Could not locate spatial coordinates. GraphST requires adata.obsm['spatial'].")

    # Ensure annotation column exists
    if args.annot not in adata_sc.obs.columns:
        raise ValueError(f"Annotation column '{args.annot}' not found in single-cell data.")

    # 3. Preprocessing (Intersection)
    adata_sc.var_names_make_unique()
    adata_sp.var_names_make_unique()

    intersect = np.intersect1d(adata_sc.var_names, adata_sp.var_names)
    adata_sc = adata_sc[:, intersect].copy()
    adata_sp = adata_sp[:, intersect].copy()
    print(f"Data intersection: {len(intersect)} common genes.")

    # --- GRAPHST MANUAL PREP ---
    # 1. Build Adjacency Matrix (adj)
    print("Constructing spatial interaction graph...")
    construct_interaction(adata_sp)
    
    # 2. Set 'feat' (Force Full Dimensions, No PCA)
    # We use the full expression matrix matching the intersection
    print("Setting up feature matrices (Full Dimension)...")
    if issparse(adata_sp.X):
        adata_sp.obsm['feat'] = adata_sp.X.toarray()
    else:
        adata_sp.obsm['feat'] = adata_sp.X.copy()

    # 3. Set 'feat_a' (Neighbor Aggregation)
    # This replaces the internal get_feature() call which we skipped.
    if 'adj' in adata_sp.obsm:
        adj = adata_sp.obsm['adj']
        # Calculate feat_a = adj * feat
        if issparse(adj):
            feat_a = adj.dot(adata_sp.obsm['feat'])
        else:
            feat_a = np.dot(adj, adata_sp.obsm['feat'])
        
        adata_sp.obsm['feat_a'] = feat_a
    else:
        raise KeyError("Adjacency matrix 'adj' not found after construct_interaction().")

    # 4. Initialize and Train GraphST
    model = GraphST(
        adata=adata_sp,
        adata_sc=adata_sc,
        device=device,
        epochs=args.epochs,
        datatype=args.datatype,
        deconvolution=True
    )

    print(f"Training GraphST model (Type: {args.datatype})...")
    start_time = time.time()
    
    # Train and learn mapping
    adata_sp, adata_sc = model.train_map()
    
    end_time = time.time()
    print(f"Training finished in {end_time - start_time:.2f} seconds.")

    # 5. Extract and Aggregate Results
    map_matrix = adata_sp.obsm['map_matrix'] 
    sc_labels = adata_sc.obs[args.annot].values
    
    print("Aggregating cell-to-spot mapping into cell type proportions...")
    unique_labels = np.unique(sc_labels)
    proportions = pd.DataFrame(index=adata_sp.obs_names, columns=unique_labels)
    
    for label in unique_labels:
        indices = np.where(sc_labels == label)[0]
        if len(indices) > 0:
            proportions[label] = map_matrix[:, indices].sum(axis=1)
        else:
            proportions[label] = 0.0

    # 6. Save Output
    proportions = proportions.div(proportions.sum(axis=1), axis=0).fillna(0)
    
    print(f"Saving results to {args.output}...")
    proportions.to_csv(args.output, sep="\t", index=True, header=True)

if __name__ == "__main__":
    main()
