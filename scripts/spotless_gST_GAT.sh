#!/bin/bash
#SBATCH --job-name=run_gst_gat
#SBATCH --partition=GPU-shared
#SBATCH --gres=gpu:h100-80:1
#SBATCH --time=03:59:00
#SBATCH --account=cis250160p
#SBATCH --mem=60G
#SBATCH --chdir=/ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark
#SBATCH --output=/ocean/projects/cis250160p/xzhaoa/GNN_gST/logs/%x_%j.log
#SBATCH --error=/ocean/projects/cis250160p/xzhaoa/GNN_gST/logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=xinruoz@andrew.cmu.edu

cd /ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark

# --- 1. DEFINE METHOD PARAMETERS HERE ---
MODEL_TYPE="GAT"

# Define Output Directory
OUTDIR="${ROOTDIR:-/ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark}/results_${MODEL_TYPE}"

# --- SETUP --------------------------------------------------------
module purge
module load nextflow
module load anaconda3/2024.10-1
module load singularity 2>/dev/null || true

mkdir -p logs
mkdir -p "$OUTDIR"

# --- CACHE CONFIGURATION ------------------------------------------
ROOTDIR="/ocean/projects/cis250160p/xzhaoa/GNN_gST/spotless-benchmark"
export NXF_SINGULARITY_CACHEDIR="$ROOTDIR/.nf_singularity_cache"
export APPTAINER_TMPDIR="${LOCAL:-/tmp}/apptainer_tmp"
export APPTAINER_CACHEDIR="${LOCAL:-/tmp}/apptainer_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

# --- RUN PIPELINE -------------------------------------------------
STANDARDS_ROOT="$ROOTDIR/standards"

# Define all dataset pairs
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

echo "========================================================"
echo ">>> RUNNING BATCH: $MODEL_TYPE"
echo ">>> OUTPUT DIR:    $OUTDIR"
echo "========================================================"

# LOOP through each pair
for pair in "${pairs[@]}"; do
    read -r r s <<< "$pair"
    
    if [ -f "$STANDARDS_ROOT/$r" ]; then
        SC_FILE="$STANDARDS_ROOT/$r"
    else
        SC_FILE="$STANDARDS_ROOT/reference/$r"
    fi
    
    ID=$(basename "$r" .rds)
    
    echo ">>> Processing Dataset: $ID"

    nextflow run main.nf -resume \
        -profile singularity \
        -c fix_paths.config \
        -w "work_${MODEL_TYPE}" \
        --rootdir "$ROOTDIR" \
        --mode run_dataset \
        --sc_input "$SC_FILE" \
        --sp_input "$STANDARDS_ROOT/$s" \
        --annot celltype \
        --outdir.props "$OUTDIR" \
        --outdir.metrics "$OUTDIR" \
        --methods graphst_custom \
        --graphst_model_type "$MODEL_TYPE" \
        --epochs 2 \
        --gpu \
        -with-report   "logs/report_${MODEL_TYPE}_${ID}.html" \
        -with-timeline "logs/timeline_${MODEL_TYPE}_${ID}.html" \
        -with-trace    "logs/trace_${MODEL_TYPE}_${ID}.txt"
        
done