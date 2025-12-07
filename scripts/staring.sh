module load anaconda3/2024.10-1
module load cuda/12.4.0
conda activate c2l_gnn

srun --partition=GPU-shared --gres=gpu:v100-32:1 --time=0:40:00 --account=cis250160p --pty bash

model = GraphST(adata, datatype='GATv2', device=device)
['GATv2', 'GAT', 'SGC']

pairs=(
  "gold_standard_1.rds gold_standard_1/*.rds"
  "gold_standard_2.rds gold_standard_2/*.rds"
  "gold_standard_3_19celltypes.rds gold_standard_3/*.rds"
  "silver_standard_1_brain_cortex.rds silver_standard_1-*/*.rds"
  "silver_standard_2_cerebellum.rds silver_standard_2-*/*.rds"
  "silver_standard_3_hippocampus.rds silver_standard_3-*/*.rds"
  "silver_standard_4_kidney.rds silver_standard_4-*/*.rds"
  "silver_standard_5_scc_p5.rds silver_standard_5-*/*.rds"
  "silver_standard_6_melanoma.rds silver_standard_6-*/*.rds"
)
