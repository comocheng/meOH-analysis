#!/bin/bash
#SBATCH --job-name=mk_parametric_uncertainty_plots
#SBATCH --error=error_plots.log
#SBATCH --mail-user=blais.ch@northeastern.edu 
#SBATCH --mail-type=FAIL,END

python-jl make_plots_density.py

