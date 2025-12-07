#!/bin/bash
#SBATCH --job-name=run_gst_fixed
#SBATCH --partition=GPU-shared
#SBATCH --gres=gpu:v100-32:1
#SBATCH --time=2:00:00
#SBATCH --account=cis250160p
#SBATCH --mem=60G
#SBATCH --chdir=/ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark
#SBATCH --output=/ocean/projects/cis250160p/xzhaoa/GNN_gST/logs/%x_%j.log
#SBATCH --error=/ocean/projects/cis250160p/xzhaoa/GNN_gST/logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=xinruoz@andrew.cmu.edu

cd /ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark

# --- SETUP ---
module purge
module load nextflow

# 1. LOAD ANACONDA (CRITICAL FIX)
# Nextflow needs this to access your 'graphst_env' defined in the config
module load anaconda3/2024.10-1

# Load Singularity (continue if it fails/is pre-loaded)
module load singularity 2>/dev/null || true

mkdir -p logs

# --- PATHS ---
ROOTDIR="/ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark"
export NXF_SINGULARITY_CACHEDIR="$ROOTDIR/.nf_singularity_cache"
export APPTAINER_TMPDIR="${LOCAL:-/tmp}/apptainer_tmp"
export APPTAINER_CACHEDIR="${LOCAL:-/tmp}/apptainer_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

# --- RUN PIPELINE ---
# --- RUN PIPELINE ---
# Fix file locations
STANDARDS_ROOT="$ROOTDIR/standards"

# Check config
if [ ! -f "fix_paths.config" ]; then
    echo "Error: fix_paths.config not found!"
    exit 1
fi

# Define all dataset pairs
# Define all dataset pairs
# Indices are shifted because Cerebellum has two parts (Silver 2 & 3)
# Define all dataset pairs based on your file list
pairs=(
  # --- GOLD STANDARDS ---
  "gold_standard_1.rds gold_standard_1/*.rds"
  "gold_standard_2.rds gold_standard_2/*.rds"
  "gold_standard_3_19celltypes.rds gold_standard_3/*.rds"

  # --- SILVER STANDARDS ---
  # 1. Brain Cortex
  "silver_standard_1_brain_cortex.rds silver_standard_1-*/*.rds"
  
  # 2. Cerebellum (Cell)
  "silver_standard_2_cerebellum_cell.rds silver_standard_2-*/*.rds"
  
  # 3. Cerebellum (Nucleus) - Now included since it exists in your folder
  "silver_standard_3_cerebellum_nucleus.rds silver_standard_3-*/*.rds"

  # 4. Hippocampus
  "silver_standard_4_hippocampus.rds silver_standard_4-*/*.rds"
  
  # 5. Kidney
  "silver_standard_5_kidney.rds silver_standard_5-*/*.rds"
  
  # 6. SCC P5
  "silver_standard_6_scc_p5.rds silver_standard_6-*/*.rds"
)

# LOOP through each pair
for pair in "${pairs[@]}"; do
    read -r r s <<< "$pair"
    
    # Dynamic path construction
    if [ -f "$STANDARDS_ROOT/$r" ]; then
        SC_FILE="$STANDARDS_ROOT/$r"
    else
        SC_FILE="$STANDARDS_ROOT/reference/$r"
    fi
    
    # Generate a unique ID for logs (remove .rds and spaces)
    ID=$(basename "$r" .rds)
    
    echo "========================================================"
    echo ">>> PROCESSING: $ID"
    echo ">>> SC Reference: $SC_FILE"
    echo ">>> SP Pattern:   $STANDARDS_ROOT/$s"
    echo "========================================================"

    # Run Nextflow for this dataset
    # We use a separate trace/report file for each dataset to avoid overwriting
    nextflow run main.nf -resume \
        -profile singularity \
        -c fix_paths.config \
        --rootdir "$ROOTDIR" \
        --mode run_dataset \
        --sc_input "$SC_FILE" \
        --sp_input "$STANDARDS_ROOT/$s" \
        --annot celltype \
        --methods graphst_custom \
        --graphst_model_type 10X \
        --epochs 600 \
        --gpu \
        -with-report   "logs/report_graphst_${ID}.html" \
        -with-timeline "logs/timeline_graphst_${ID}.html" \
        -with-trace    "logs/trace_graphst_${ID}.txt"
        
    # Clean work directory between large runs to save space (Optional but recommended)
    # nextflow clean -f -k last
done