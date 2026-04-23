#!/bin/bash

#SBATCH -p compute
#SBATCH -w node15
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=32G
#SBATCH -J analysis
#SBATCH -o analysis_%j.out
#SBATCH -e analysis_%j.err
#SBATCH --time=48:00:00

# Load conda environment support
#source ~/anaconda3/etc/profile.d/conda.sh

# Activate your conda environment (already created)
#conda activate swalih_gatk      # change name if different

# Run the R script
Rscript sft_new.R

