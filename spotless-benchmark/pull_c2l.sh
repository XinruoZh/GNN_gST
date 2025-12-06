#!/bin/bash
#SBATCH --job-name=pull_spotless
#SBATCH --partition=RM-shared
#SBATCH --time=02:30:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2000M
#SBATCH --account=bio230007p
#SBATCH --output=pull_spotless_%j.out
#SBATCH --error=pull_spotless_%j.err

set -euo pipefail

module purge
module load nextflow   # brings apptainer

PROJECT_DIR="/ocean/projects/cis250160p/xzhaoa/st_cell2location/spotless-benchmark"
cd "$PROJECT_DIR"

echo "Job started on: $(hostname)"
echo "Time: $(date)"
echo "==========================================="

CACHE_DIR="$PROJECT_DIR/.nf_singularity_cache"
mkdir -p "$CACHE_DIR"

# Use project space for Apptainer cache/temp
export APPTAINER_CACHEDIR="$PROJECT_DIR/.apptainer_cache"
export APPTAINER_TMPDIR="$PROJECT_DIR/.apptainer_tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

pull_img () {
    name=$1
    docker_url=$2
    sif_path="$CACHE_DIR/$name.sif"

    echo "-------------------------------------------"
    echo "Checking: $name"
    echo "Docker:   $docker_url"
    echo "SIF path: $sif_path"
    echo "-------------------------------------------"

    if [[ -f "$sif_path" ]]; then
        echo "[SKIP] Already exists: $sif_path"
    else
        echo "[PULL] Downloading $docker_url → $sif_path"
        apptainer pull "$sif_path" "docker://$docker_url"
    fi
}

echo "Pulling Spotless containers..."
echo ""

# only the two you need right now
pull_img "csangara-seuratdisk-latest"       "csangara/seuratdisk:latest"
pull_img "csangara-sp_cell2location-latest" "csangara/sp_cell2location:latest"

echo ""
echo "==========================================="
echo "Cache contents:"
ls -lh "$CACHE_DIR"
echo "Finished at: $(date)"


