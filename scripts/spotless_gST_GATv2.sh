#!/bin/bash
#SBATCH --job-name=run_gst_gatv2
#SBATCH --partition=GPU-shared
#SBATCH --gres=gpu:v100-32:1
#SBATCH --time=0:10:00
#SBATCH --account=cis250160p
#SBATCH --mem=60G
#SBATCH --chdir=/ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark
#SBATCH --output=/ocean/projects/cis250160p/xzhaoa/GNN_gST/log/%x_%j_gatv2.log
#SBATCH --error=/ocean/projects/cis250160p/xzhaoa/GNN_gST/log/%x_%j_gatv2.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=xinruoz@andrew.cmu.edu

# ------------------------------------------------------------------
# SETUP
# ------------------------------------------------------------------
cd /ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark
module purge
module load anaconda3/2024.10-1
module load nextflow/21.10.6

# Ensure Singularity is available (Bridges-2 default)
if ! command -v singularity &> /dev/null; then
    echo "WARNING: Singularity not found in path. Trying module load..."
    module load singularity
fi

mkdir -p log

# ------------------------------------------------------------------
# CACHE CONFIGURATION
# ------------------------------------------------------------------
# Using the specific cache folder for this project
export NXF_SINGULARITY_CACHEDIR="/ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark/.nf_singularity_cache"
mkdir -p $NXF_SINGULARITY_CACHEDIR

# Use Local NVMe for temporary build files (Faster performance)
export APPTAINER_TMPDIR="$LOCAL/apptainer_tmp"
export SINGULARITY_TMPDIR="$LOCAL/apptainer_tmp"
export APPTAINER_CACHEDIR="$LOCAL/apptainer_cache"
export SINGULARITY_CACHEDIR="$LOCAL/apptainer_cache"

mkdir -p "$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_CACHEDIR"

# ------------------------------------------------------------------
# RUN PIPELINE
# ------------------------------------------------------------------
ROOTDIR="/ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark"
REFERENCE="$ROOTDIR/standards/reference"
STANDARDS="$ROOTDIR/standards"

# Define the dataset pair
pair="gold_standard_1.rds gold_standard_1/*.rds"
read -r r s <<< "$pair"

echo ">>> Running GraphST (GATv2) on $r with spatial $s"

# RUN COMMAND
# -profile singularity: Matches your working script (avoids VSC cluster config issues in 'hpc' profile)
# --graphst_model_type 10X: Runs the standard GCN version
nextflow run main.nf -resume \
    -profile singularity \
    --rootdir "$ROOTDIR" \
    --mode run_dataset \
    --sc_input "$REFERENCE/$r" \
    --sp_input "$STANDARDS/$s" \
    --annot celltype \
    --methods graphst_custom \
    --graphst_model_type GATv2 \
    --epochs 2 \
    --gpu \
    -with-report log/report_graphst_gatv2.html \
    -with-timeline log/timeline_graphst_gatv2.html \
    -with-trace log/trace_graphst_gatv2.txt