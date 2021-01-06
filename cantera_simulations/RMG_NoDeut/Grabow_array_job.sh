#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=20GB
#SBATCH --exclude=c5003
#SBATCH --job-name=PFR
#SBATCH --output=nodeut_grabow.%a.log
#SBATCH --error=nodeut_grabow.%a.slurm.log
#SBATCH --partition=short

#an array for the job.
#SBATCH --array=1-36


####################################################
source activate rmg_env
python -u Grabow_reactor_nodeut.py