#!/bin/bash
#SBATCH --job-name=train_bpnet
#SBATCH --output=train_bpnet_%j.out
#SBATCH --error=train_bpnet_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=256G
#SBATCH --partition=scu-gpu
#SBATCH --cpus-per-task=12  # adjust if needed

# Activate your conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bpnet

# Run the command
bpnet train dataspec_test.yml ./output  --num-workers 12 --gpu=1 --override='train.epochs=10' --overwrite
