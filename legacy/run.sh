#!/bin/bash

#SBATCH --job-name=boltz2.0.3
#SBATCH --time=24:00:00
#SBATCH -p SPBe_gpu

source /prd/pkgs/miniconda/conda3-py38/bashrc
conda activate boltz2.0.3

python screening.py
echo 'fin'