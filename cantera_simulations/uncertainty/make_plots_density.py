import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize 
from matplotlib import cm
from scipy.stats import gaussian_kde
import os
import glob
import re

# rmg_runs_dir = "/home/moon/methanol/perturb_5000/"
rmg_runs_dir = "/scratch/westgroup/methanol/perturb_5000/"

csv_files = glob.glob(os.path.join(rmg_runs_dir, 'run_****', 'cantera', 'ct_analysis.csv'))
# print(csv_files)

T = 'all'
P = 3039750.0000000005
V = 3.32416e-05

all_fields = [
#               'T (K) final',
              'N2', 'Ne',
              'H2(2)', 'CO(3)', 'CO2(4)', 'H2O(5)', 'CH2O(6)', 'HCOOH(7)', 'CH3OH(8)',
              'HCOOCH3(9)', 'CH4(24)', 'X(1)', 'H*(10)', 'O*(11)', 'OH*(12)',
              'H2O*(13)', 'CO*(14)', 'CO2*(15)', 'HCO*(16)', 'HCOO*(17)', 'COOH*(18)',
              'HCOOH*(19)', 'CH2O*(20)', 'CH3O*(21)', 'CH3O2*(22)', 'CH3OH*(23)',
              'H2X(25)', 'COX2(26)', 'SX(56)(27)', 'CHOX2(28)', 'CH3OX(29)',
              'CH2OX2(30)', 'SX(31)', 'SX(32)', 'SX(33)', 'C2HO2X(34)', 'CH4X(35)',
              'CH3X(36)', 'SX(37)', 'SX(62)(38)', 'SX(58)(39)', 'CHO2X2(40)',
              'C2H4OX(41)', 'SX(42)', 'C2H5OX(43)', 'SX(44)', 'C2H5OX(45)',
              #'SX(116)'
]

Imax = 5000
for field in all_fields:
    
    field_to_plot_fname = re.sub(r'[^a-zA-Z0-9]', '', field)

    # get SS MeOH concentration vs temperature
    first_run = True
    for i, csv_file in enumerate(csv_files):
        df = pd.read_csv(csv_file)
        if first_run:
            first_run = False
            T = df["T (K)"].values
            x_methanol = df[field].values
        else: 
            T = np.append(T,df["T (K)"].values)
            x_methanol = np.append(x_methanol,df[field].values)
        
        
    if i>Imax:
        break
            
    # Calculate the point density
    xy = np.vstack([T,x_methanol])
    z = gaussian_kde(xy)(xy) 
    
    # Make scatterplot
    fig, ax = plt.subplots()
    cax = ax.scatter(T, x_methanol, c=z, s=300)
       
    # plot the density graph
    plt.xlabel('Temperature (K)')
    plt.ylabel("Mole Fraction " + field)
    
    # make a bar that contains a colormap to the point density. 
    norm = Normalize(vmin = np.min(z), vmax = np.max(z))
    cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
    cbar.ax.set_ylabel('Point Density')
    
    plt.savefig(f'{field_to_plot_fname}_{Imax}.png')
    plt.clf()


# columns for reference:
# ['Unnamed: 0', 'time (s)', 'T (K)', 'P (Pa)', 'V (m^3/s)',
#        'x_CO initial', 'x_CO2 initial', 'x_H2 initial', 'x_H2O initial',
#        'CO2/(CO2+CO)', '(CO+CO2/H2)', 'T (C) final', 'Rtol', 'Atol',
#        'reactor type', 'energy on?', 'catalyst weight (kg)', 'N2', 'Ne',
#        'H2(2)', 'CO(3)', 'CO2(4)', 'H2O(5)', 'CH2O(6)', 'HCOOH(7)', 'CH3OH(8)',
#        'HCOOCH3(9)', 'CH4(24)', 'X(1)', 'H*(10)', 'O*(11)', 'OH*(12)',
#        'H2O*(13)', 'CO*(14)', 'CO2*(15)', 'HCO*(16)', 'HCOO*(17)', 'COOH*(18)',
#        'HCOOH*(19)', 'CH2O*(20)', 'CH3O*(21)', 'CH3O2*(22)', 'CH3OH*(23)',
#        'H2X(25)', 'COX2(26)', 'SX(56)(27)', 'CHOX2(28)', 'CH3OX(29)',
#        'CH2OX2(30)', 'SX(31)', 'SX(32)', 'SX(33)', 'C2HO2X(34)', 'CH4X(35)',
#        'CH3X(36)', 'SX(37)', 'SX(62)(38)', 'SX(58)(39)', 'CHO2X2(40)',
#        'C2H4OX(41)', 'SX(42)', 'C2H5OX(43)', 'SX(44)', 'C2H5OX(45)',
#        'SX(116)']
