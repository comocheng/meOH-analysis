#!/bin/bash
#SBATCH --job-name=parametric_uncertainty_analysis1
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH --nodes=1
#SBATCH --partition=west
#SBATCH --exclude=c5003
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --mail-user="$USER"@northeastern.edu 
#SBATCH --mail-type=FAIL,END

source activate rmg_julia_env
python-jl /scratch/westgroup/methanol/meOH-analysis/cantera_simulations/uncertainty/run_graaf_slurm_scripts.py


