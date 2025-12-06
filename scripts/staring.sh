module load anaconda3/2024.10-1
module load cuda/12.4.0
conda activate c2l_gnn

srun --partition=GPU-shared --gres=gpu:v100-32:1 --time=0:40:00 --account=cis250160p --pty bash

model = GraphST(adata, datatype='GATv2', device=device)
['GATv2', 'GAT', 'SGC']
