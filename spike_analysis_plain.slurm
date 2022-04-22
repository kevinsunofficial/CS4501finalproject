#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --job-name="Spike Analysis Plain"
#SBATCH -o slurm/spike_analysis_plain.out
#SBATCH -e slurm/spike_analysis_plain.error
#SBATCH --partition=gpu
#SBATCH --gres=gpu:4

source /etc/profile.d/modules.sh

python spike_analysis_plain.py
