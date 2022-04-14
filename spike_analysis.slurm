#!/bin/bash
#SBATCH -t 72:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH -o slurm/spike_analysis.out
#SBATCH -e slurm/spike_analysis.error
#SBATCH -p parellel

