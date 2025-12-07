#!/bin/bash
#SBATCH --job-name=convert_dlpfc
#SBATCH --partition=GPU-shared        # Must use GPU partition for cis250160p
#SBATCH --gpus=v100-32:1              # Request 1 GPU to get access to the node
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10            # Request reasonable CPU power
#SBATCH --mem=60G                     # Request 60GB RAM (plenty for conversion)
#SBATCH --time=00:30:00
#SBATCH --account=cis250160p
#SBATCH --output=convert_dlpfc.log
#SBATCH --error=convert_dlpfc.err

# Load modules
module purge
module load anaconda3/2024.10-1

# Activate environment
source activate download_env

# Run conversion
cd /ocean/projects/cis250160p/xzhaoa/GNN_gST
echo "Starting conversion on GPU node..."
Rscript scripts/convert_dlpfc.R
echo "Conversion finished."