#!/bin/bash
#SBATCH --job-name=run_gst_fixed
#SBATCH --partition=GPU-shared
#SBATCH --gres=gpu:v100-32:1
#SBATCH --time=0:30:00                 # INCREASED TIME for a stable test run
#SBATCH --account=cis250160p
#SBATCH --mem=32G                      # Reduced memory to 32G (safer for shared partition)
#SBATCH --chdir=/ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark
#SBATCH --output=/ocean/projects/cis250160p/xzhaoa/GNN_gST/log/%x_%j.log
#SBATCH --error=/ocean/projects/cis250160p/xzhaoa/GNN_gST/log/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=xinruoz@andrew.cmu.edu

# --- PATHS (Simplified) -------------------------------------------
ROOTDIR="/ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark"
PYTHON_SCRIPT="$ROOTDIR/subworkflows/deconvolution/GraphST/script_nf.py"
DATA_DIR="$ROOTDIR/standards"

# --- ENVIRONMENT SETUP (The Fix) ----------------------------------
# 1. Clean environment
module purge
conda deactivate 2>/dev/null || true 

# 2. Load and activate the host environment (where POT is installed)
module load anaconda3/2024.10-1
source activate your_env_name         # <-- CRITICAL: REPLACE with your actual environment name

mkdir -p log

# ------------------------------------------------------------------
# RUN GNN CODE DIRECTLY (Bypassing Nextflow/Singularity)
# ------------------------------------------------------------------
# Note: Nextflow usually converts .rds to .h5ad first. We assume the .h5ad files are present
# from your previous successful conversion processes.
echo ">>> Running GraphST (GCN) directly on V100"

# Define the dataset pairs (The list you provided)
pairs=(
  "gold_standard_1.rds gold_standard_1/*.rds"
)

# Define paths (Using variables from your fixed script)
ROOTDIR="/ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark"
PYTHON_SCRIPT="$ROOTDIR/subworkflows/deconvolution/GraphST/script_nf.py"
STANDARDS="$ROOTDIR/standards"

# Loop through each dataset definition pair
for pair in "${pairs[@]}"; do
    # Read the single-cell reference (REF) and the spatial wildcard pattern (SP_PATTERN)
    read -r REF_RDS SP_PATTERN <<< "$pair"
    
    # 1. Determine the reference H5AD file name
    # Example: gold_standard_1.rds -> gold_standard_1.h5ad
    REF_H5AD="${REF_RDS/.rds/.h5ad}"

    echo "--- Starting analysis for Reference: $REF_H5AD ---"

    # 2. Loop through all individual spatial files matching the wildcard
    # The shell expands "$STANDARDS/$SP_PATTERN" to all matching paths
    for SP_FILE_RDS in $STANDARDS/$SP_PATTERN; do
        
        # 3. Derive the required H5AD file name and the output identifier
        SP_FILE_H5AD="${SP_FILE_RDS/.rds/.h5ad}"
        
        # Extract the base filename (e.g., fov4.h5ad) for the output file
        OUTPUT_ID=$(basename "$SP_FILE_H5AD" .h5ad)
        
        echo "   -> Processing spatial sample: $OUTPUT_ID"

        # 4. EXECUTE PYTHON COMMAND (Targeting the V100 via CUDA)
        python "$PYTHON_SCRIPT" \
            "$STANDARDS/$REF_H5AD" \
            "$SP_FILE_H5AD" \
            celltype \
            --datatype 10X \
            --device cuda \
            --epochs 2 \
            --output "proportions_graphst_${OUTPUT_ID}.tsv"
            
        # Optional: Add error checking here if necessary
    done
done

echo "Loop finished. Check output files in the current directory."