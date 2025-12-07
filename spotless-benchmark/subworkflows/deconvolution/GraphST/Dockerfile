FROM python:3.8-slim

# 1. Install system dependencies
RUN apt-get update && apt-get install -y build-essential git procps

# 2. Install PyTorch (Specific version recommended for stability)
# This installs PyTorch 2.0.1 with CUDA 11.8 support (common on HPCs like Bridges-2)
# If Bridges-2 uses a different CUDA version, change 'cu118' to 'cu121' etc.
RUN pip install torch==2.0.1 torchvision==0.15.2 torchaudio==2.0.2 --index-url https://download.pytorch.org/whl/cu118

# 3. Install Torch Geometric Dependencies (The Critical Part!)
# We use the -f flag to find the pre-compiled wheels for our specific Torch version
RUN pip install torch-scatter torch-sparse torch-cluster torch-spline-conv \
    -f https://data.pyg.org/whl/torch-2.0.1+cu118.html

# 4. Install Torch Geometric itself (Contains GAT, SGC, GATv2)
RUN pip install torch-geometric

# 5. Install GraphST dependencies
RUN pip install scanpy==1.9.1 anndata==0.8.0 pandas numpy scipy scikit-learn matplotlib==3.6.0

# 6. Add your modified GraphST code
# This copies everything from your current folder (./) to /app/GraphST inside the container
COPY . /app/GraphST

# 7. Install the custom GraphST
WORKDIR /app/GraphST
RUN pip install .

# 8. Set setup for execution
WORKDIR /app
ENV PYTHONPATH="${PYTHONPATH}:/app"