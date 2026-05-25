#!/usr/bin/env bash
set -euo pipefail
echo "===== Step H2w2 v3 WSL / fastp preflight ====="
echo "WSL user: $(whoami)"
echo "WSL home: $HOME"
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniforge3/etc/profile.d/conda.sh"
elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/anaconda3/etc/profile.d/conda.sh"
else
  echo "ERROR: conda.sh not found under ~/miniforge3, ~/miniconda3, or ~/anaconda3"
  exit 91
fi
conda activate fastp_env
echo "Conda env: $CONDA_DEFAULT_ENV"
echo "fastp path:"
which fastp
echo "fastp version:"
fastp --version
