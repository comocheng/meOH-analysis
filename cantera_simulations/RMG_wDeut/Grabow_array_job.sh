#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=20GB
#SBATCH --exclude=c5003
#SBATCH --job-name=grabDeut
#SBATCH --output=logs/deut_grabow.%a.log
#SBATCH --error=logs/deut_grabow.%a.slurm.log
#SBATCH --partition=short

#an array for the job.
#SBATCH --array=1-107


####################################################
source activate cantera
python -u Grabow_reactor_script.py