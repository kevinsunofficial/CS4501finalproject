#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH --job-name="Spike Analysis"
#SBATCH --nodes=8
#SBATCH -o slurm/spike_analysis.out
#SBATCH -e slurm/spike_analysis.error
#SBATCH --partition=main

source /etc/profile.d/modules.sh

python spike_analysis.py
