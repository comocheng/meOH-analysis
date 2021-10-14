#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=2:00:00
#SBATCH --mem=20GB
#SBATCH --exclude=c5003
#SBATCH --job-name=graaf
#SBATCH --output=logs_graaf/graaf.%a.log
#SBATCH --error=logs_graaf/graaf.%a.slurm.log
#SBATCH --partition=short
#SBATCH --mail-user=blais.ch@northeastern.edu 
#SBATCH --mail-type=FAIL,END

#an array for the job. usually 107, testing use 2
#SBATCH --array=0-216


####################################################
source activate rmg_env
python -u Graaf_sbr_script.py