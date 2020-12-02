#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --mem=20GB
#SBATCH --exclude=c5003
#SBATCH --job-name=PFR
#SBATCH --output=grabow.%a.log
#SBATCH --error=grabow.%a.slurm.log
#SBATCH --partition=west

#an array for the job.
#SBATCH --array=1-36


####################################################
source activate rmg_env
export RMGpy=/home/blais.ch/RMG_Base_env/RMG-Py
python -u Grabow_reactor.py