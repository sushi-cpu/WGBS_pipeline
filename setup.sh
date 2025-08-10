#!/bin/bash
set -euo pipefail

echo "────────────────────────────────────────────"
echo "🌱 WGBS Pipeline Setup Script (Linux & macOS)"
echo "────────────────────────────────────────────"

# Detect OS
OS_TYPE="$(uname)"
ARCH_TYPE="$(uname -m)"
echo "🔍 Detected OS: $OS_TYPE, Architecture: $ARCH_TYPE"

# Set installer URL based on OS
if [[ "$OS_TYPE" == "Linux" ]]; then
    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
elif [[ "$OS_TYPE" == "Darwin" ]]; then
    if [[ "$ARCH_TYPE" == "arm64" ]]; then
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
    else
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
    fi
else
    echo "❌ Unsupported OS: $OS_TYPE"
    exit 1
fi

# 1. Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "🛠️  Conda not found. Installing Miniconda..."
    INSTALLER_NAME="${MINICONDA_URL##*/}"
    wget "$MINICONDA_URL" -O "$INSTALLER_NAME"
    bash "$INSTALLER_NAME" -b -p "$HOME/miniconda"
    eval "$($HOME/miniconda/bin/conda shell.bash hook)"
    conda init
    source ~/.bashrc || source ~/.zshrc || true
    echo "✅ Miniconda installed."
else
    echo "✅ Conda already installed."
fi

# 2. Setup the conda environment
echo "📦 Creating Conda environment for WGBS pipeline..."
conda env create -f wgbs-pipeline.yml || echo "ℹ️ Environment may already exist."

echo "✅ Environment 'wgbs-pipeline' is ready."

echo "🌱 WGBS Pipeline Setup Complete!
You can now activate the environment using: 
conda activate wgbs-pipeline"
