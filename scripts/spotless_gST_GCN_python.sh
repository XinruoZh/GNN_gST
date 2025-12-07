#!/bin/bash
#SBATCH --job-name=run_gst_fixed
#SBATCH --partition=GPU-shared
#SBATCH --gres=gpu:v100-32:1
#SBATCH --time=0:10:00               # INCREASED TIME for a stable test run
#SBATCH --account=cis250160p
#SBATCH --mem=32G                    # Reduced memory to 32G (safer for shared partition)
#SBATCH --chdir=/ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark
#SBATCH --output=/ocean/projects/cis250160p/xzhaoa/GNN_gST/log/%x_%j.log
#SBATCH --error=/ocean/projects/cis250160p/xzhaoa/GNN_gST/log/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=xinruoz@andrew.cmu.edu

# --- PATHS -------------------------------------------
ROOTDIR="/ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark"
PYTHON_SCRIPT="$ROOTDIR/subworkflows/deconvolution/GraphST/script_nf.py"
DATA_DIR="$ROOTDIR/standards"

# --- ENVIRONMENT SETUP ----------------------------------
module purge
module load anaconda3/2024.10-1

# Activate the environment located at the path you provided
source activate /ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark/subworkflows/deconvolution/GraphST/graphst_env

mkdir -p log

# ------------------------------------------------------------------
# RUN GNN CODE DIRECTLY (Bypassing Nextflow/Singularity)
# ------------------------------------------------------------------
echo ">>> Running GraphST (GCN) directly on V100"

# Define the dataset pairs
pairs=(
  "gold_standard_1.rds gold_standard_1/*.rds"
)

# Loop through each dataset definition pair
for pair in "${pairs[@]}"; do
    read -r REF_RDS SP_PATTERN <<< "$pair"
    
    # 1. Determine the reference H5AD file name
    REF_H5AD="${REF_RDS/.rds/.h5ad}"

    echo "--- Starting analysis for Reference: $REF_H5AD ---"

    # 2. Loop through all individual spatial files matching the wildcard
    for SP_FILE_RDS in $DATA_DIR/$SP_PATTERN; do
        
        # 3. Derive the required H5AD file name and the output identifier
        SP_FILE_H5AD="${SP_FILE_RDS/.rds/.h5ad}"
        OUTPUT_ID=$(basename "$SP_FILE_H5AD" .h5ad)
        
        echo "   -> Processing spatial sample: $OUTPUT_ID"

        # 4. EXECUTE PYTHON COMMAND
        python "$PYTHON_SCRIPT" \
            "$DATA_DIR/$REF_H5AD" \
            "$SP_FILE_H5AD" \
            celltype \
            --datatype 10X \
            --device cuda \
            --epochs 2 \
            --output "proportions_graphst_${OUTPUT_ID}.tsv"
            
    done
done

echo "Loop finished. Check output files in the current directory."