# MM/GBSA Analysis Docker Image
# Based on CUDA-enabled OpenMM environment

FROM nvidia/cuda:11.8-devel-ubuntu20.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV CONDA_AUTO_UPDATE_CONDA=false

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    cmake \
    pkg-config \
    libffi-dev \
    libssl-dev \
    libbz2-dev \
    libreadline-dev \
    libsqlite3-dev \
    libncurses5-dev \
    libncursesw5-dev \
    xz-utils \
    tk-dev \
    libxml2-dev \
    libxmlsec1-dev \
    libffi-dev \
    liblzma-dev \
    libgdbm-compat-dev \
    libnss3-dev \
    libatlas-base-dev \
    libblas-dev \
    liblapack-dev \
    libhdf5-dev \
    libnetcdf-dev \
    libgdal-dev \
    libproj-dev \
    libgeos-dev \
    libgsl-dev \
    libfftw3-dev \
    libgmp-dev \
    libmpfr-dev \
    libmpc-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh

# Add conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

# Create conda environment
RUN conda create -n mmgbsa python=3.9 -y

# Activate environment and install packages
SHELL ["conda", "run", "-n", "mmgbsa", "/bin/bash", "-c"]

# Install conda packages
RUN conda install -c conda-forge -c openeye -c omnia \
    openmm=8.0 \
    mdtraj=1.9.8 \
    openff-toolkit=0.14.0 \
    openmmforcefields=1.0.0 \
    numpy=1.21.0 \
    pandas=1.3.0 \
    scipy=1.7.0 \
    matplotlib=3.5.0 \
    seaborn=0.11.0 \
    pyyaml=6.0 \
    jupyter \
    ipython \
    -y

# Install additional pip packages
RUN pip install \
    pathlib2 \
    openff-forcefields \
    nglview \
    mbuild \
    foyer

# Set working directory
WORKDIR /workspace

# Copy MM/GBSA package files
COPY mmgbsa_v3.py /workspace/
COPY mmgbsa_v3_entropy.py /workspace/
COPY NormalModeAnalysis.py /workspace/
COPY per_residue_decompose.py /workspace/
COPY mmgbsa_runner.py /workspace/
COPY requirements.txt /workspace/
COPY README.md /workspace/

# Create test directory and copy test files
RUN mkdir -p /workspace/test
COPY test/ /workspace/test/

# Create output directory
RUN mkdir -p /workspace/mmgbsa_results

# Set environment variables for OpenMM
ENV OPENMM_DEFAULT_PLATFORM=CUDA
ENV CUDA_VISIBLE_DEVICES=0

# Create entrypoint script
RUN echo '#!/bin/bash\n\
source /opt/conda/etc/profile.d/conda.sh\n\
conda activate mmgbsa\n\
exec "$@"' > /entrypoint.sh \
    && chmod +x /entrypoint.sh

# Set entrypoint
ENTRYPOINT ["/entrypoint.sh"]

# Default command
CMD ["python", "mmgbsa_runner.py", "--help"]

# Expose port for Jupyter (optional)
EXPOSE 8888

# Add labels
LABEL maintainer="MM/GBSA Development Team"
LABEL description="MM/GBSA Analysis Package with CUDA Support"
LABEL version="0.0.4"
LABEL org.opencontainers.image.source="https://github.com/your-repo/mmgbsa" 