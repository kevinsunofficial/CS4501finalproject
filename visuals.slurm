#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --job-name="Spike visuals"
#SBATCH -o slurm/spike_visuals.out
#SBATCH -e slurm/spike_visuals.error
#SBATCH --partition=gpu
#SBATCH --gres=gpu:4

source /etc/profile.d/modules.sh

python visuals.py