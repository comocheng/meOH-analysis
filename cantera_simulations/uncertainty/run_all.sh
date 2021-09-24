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


python /scratch/westgroup/methanol/meOH-synthesis/perturbed_runs/run_slurm_scripts.py

