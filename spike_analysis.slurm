#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --job-name="Spike Analysis"
#SBATCH --nodelist=ai0[1-5]
#SBATCH -o slurm/spike_analysis.out
#SBATCH -e slurm/spike_analysis.error
#SBATCH --partition=gpu
#SBATCH --gres=gpu:4

source /etc/profile.d/modules.sh

python spike_analysis.py
