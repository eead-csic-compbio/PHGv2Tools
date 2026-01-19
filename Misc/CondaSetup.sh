#!/bin/bash

# PHGv2Tools Environment Setup Script
# This script creates a conda environment with all required dependencies

ENV_NAME="phgtools"

echo "=== PHGv2Tools Environment Setup ==="

# Check if we're already inside the target environment
if [[ "$CONDA_DEFAULT_ENV" == "$ENV_NAME" ]]; then
    echo "You are already inside the '$ENV_NAME' environment."
    echo "Installing/updating dependencies in the current environment..."
    
    # Install dependencies via conda (includes proper C libraries)
    echo "Installing dependencies via conda..."
    conda install matplotlib pandas scipy numpy seaborn openpyxl tqdm rich pillow ncurses -y -c conda-forge
    
    # Check for samtools
    if ! command -v samtools &> /dev/null; then
        echo "Installing samtools..."
        conda install samtools -y -c bioconda -c conda-forge
    fi
    
    # Get the script directory to find PHGv2Tools
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PHGTOOLS_DIR="$(dirname "$SCRIPT_DIR")"
    
    # Install PHGv2Tools in editable mode if pyproject.toml exists
    if [[ -f "$PHGTOOLS_DIR/pyproject.toml" ]]; then
        echo "Installing PHGv2Tools package..."
        pip install -e "$PHGTOOLS_DIR" --no-deps
    else
        echo "Warning: pyproject.toml not found at $PHGTOOLS_DIR"
        echo "You may need to install PHGv2Tools manually with: pip install -e /path/to/PHGv2Tools"
    fi
    
    echo ""
    echo "=== Setup Complete ==="
    echo "To verify installation, run: phgtools --check-setup"
    exit 0
fi

# Check if the environment already exists and remove it if it does
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "The environment $ENV_NAME already exists. Removing it..."
    conda env remove -n $ENV_NAME -y
fi

# Create a new conda environment with Python and samtools
echo "Creating conda environment: $ENV_NAME"
conda create -n $ENV_NAME python=3.10 samtools ncurses -y -c bioconda -c conda-forge

# Activate the environment
echo "Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $ENV_NAME

# Install dependencies via conda (includes proper C libraries)
echo "Installing dependencies via conda..."
conda install matplotlib pandas scipy numpy seaborn openpyxl tqdm rich pillow -y -c conda-forge

# Get the script directory to find PHGv2Tools
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PHGTOOLS_DIR="$(dirname "$SCRIPT_DIR")"

# Install PHGv2Tools in editable mode if pyproject.toml exists
if [[ -f "$PHGTOOLS_DIR/pyproject.toml" ]]; then
    echo "Installing PHGv2Tools package..."
    pip install -e "$PHGTOOLS_DIR" --no-deps
else
    echo "Warning: pyproject.toml not found at $PHGTOOLS_DIR"
    echo "You may need to install PHGv2Tools manually with: pip install -e /path/to/PHGv2Tools"
fi

echo ""
echo "=== Setup Complete ==="
echo "Environment name: $ENV_NAME"
echo "To activate the environment, run: conda activate $ENV_NAME"
echo "To verify installation, run: phgtools --check-setup"

