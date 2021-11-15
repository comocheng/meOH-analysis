#!/bin/bash
#SBATCH --job-name=mk_parametric_uncertainty_plots
#SBATCH --error=error_plots.log
#SBATCH --mail-user=$USER@northeastern.edu 
#SBATCH --mail-type=FAIL,END

python make_plots.py

