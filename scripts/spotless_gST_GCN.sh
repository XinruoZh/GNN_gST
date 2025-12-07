#!/bin/bash
#SBATCH --job-name=run_gst
#SBATCH --partition=GPU-shared
#SBATCH --gres=gpu:v100-32:1
#SBATCH --time=0:10:00
#SBATCH --account=cis250160p
#SBATCH --mem=60G
#SBATCH --chdir=/ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark
#SBATCH --output=/ocean/projects/cis250160p/xzhaoa/GNN_gST/logs/%x_%j.log
#SBATCH --error=/ocean/projects/cis250160p/xzhaoa/GNN_gST/logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=xinruoz@andrew.cmu.edu

# ------------------------------------------------------------------
# SETUP
# ------------------------------------------------------------------
module purge
conda deactivate 2>/dev/null || true
module load nextflow

# If Singularity/Apptainer isn't on PATH, try loading it
if ! command -v singularity &> /dev/null && ! command -v apptainer &> /dev/null; then
    echo "WARNING: Singularity/Apptainer not found in PATH. Trying module load singularity..."
    module load singularity 2>/dev/null || echo "Singularity module not found; relying on container support in profile."
fi

mkdir -p logs

# ------------------------------------------------------------------
# CACHE CONFIGURATION
# ------------------------------------------------------------------
ROOTDIR="/ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark"
cd /ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark
# Nextflow Singularity cache (project-specific)
export NXF_SINGULARITY_CACHEDIR="$ROOTDIR/.nf_singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

# Use node-local NVMe for Apptainer/Singularity temp + cache (if $LOCAL is defined)
export APPTAINER_TMPDIR="${LOCAL:-/tmp}/apptainer_tmp"
export SINGULARITY_TMPDIR="${LOCAL:-/tmp}/apptainer_tmp"
export APPTAINER_CACHEDIR="${LOCAL:-/tmp}/apptainer_cache"
export SINGULARITY_CACHEDIR="${LOCAL:-/tmp}/apptainer_cache"

mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

# ------------------------------------------------------------------
# RUN PIPELINE
# ------------------------------------------------------------------
REFERENCE="$ROOTDIR/standards/reference"
STANDARDS="$ROOTDIR/standards"

# Define the dataset pair (reference + spatial glob)
pair="gold_standard_1.rds gold_standard_1/*.rds"
read -r r s <<< "$pair"

echo ">>> Running GraphST (GCN) on reference: $r with spatial files: $s"

nextflow run main.nf -resume \
    -profile singularity \
    --rootdir "$ROOTDIR" \
    --mode run_dataset \
    --sc_input "$REFERENCE/$r" \
    --sp_input "$STANDARDS/$s" \
    --annot celltype \
    --methods graphst_custom \
    --graphst_model_type 10X \
    --epochs 2 \
    --gpu \
    -with-report   "logs/report_graphst.html" \
    -with-timeline "logs/timeline_graphst.html" \
    -with-trace    "logs/trace_graphst.txt"
