#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --job-name="Spike Analysis"
#SBATCH --nodelist=ai[01-02,04],lynx[08-12]
#SBATCH -o slurm/spike_analysis.out
#SBATCH -e slurm/spike_analysis.error
#SBATCH --partition=gpu
#SBATCH --gres=gpu:8

source /etc/profile.d/modules.sh

python spike_analysis.py
