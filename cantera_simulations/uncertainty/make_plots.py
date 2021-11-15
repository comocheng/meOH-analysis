import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
              'N2', 'Ne',
              'H2(2)', 'CO(3)', 'CO2(4)', 'H2O(5)', 'CH2O(6)', 'HCOOH(7)', 'CH3OH(8)',
              'HCOOCH3(9)', 'CH4(24)', 'X(1)', 'H*(10)', 'O*(11)', 'OH*(12)',
              'H2O*(13)', 'CO*(14)', 'CO2*(15)', 'HCO*(16)', 'HCOO*(17)', 'COOH*(18)',
              'HCOOH*(19)', 'CH2O*(20)', 'CH3O*(21)', 'CH3O2*(22)', 'CH3OH*(23)',
              'H2X(25)', 'COX2(26)', 'SX(56)(27)', 'CHOX2(28)', 'CH3OX(29)',
              'CH2OX2(30)', 'SX(31)', 'SX(32)', 'SX(33)', 'C2HO2X(34)', 'CH4X(35)',
              'CH3X(36)', 'SX(37)', 'SX(62)(38)', 'SX(58)(39)', 'CHO2X2(40)',
              'C2H4OX(41)', 'SX(42)', 'C2H5OX(43)', 'SX(44)', 'C2H5OX(45)',
              'SX(116)'
]

Imax = 5000
for field in all_fields:

    #field_to_plot_fname = re.sub(r'[^a-zA-Z0-9]', '', field_to_plot)
    field_to_plot_fname = re.sub(r'[^a-zA-Z0-9]', '', field)

    # Plot SS MeOH concentration vs temperature
    for i, csv_file in enumerate(csv_files):
        df = pd.read_csv(csv_file)
        # select the part of the dataframe we are interested in
        # since we are not doing different pressures and Volume flowrates currently, 
        # commented these parts out.
        # data_slice = df[df["P (Pa)"] == P]
        # data_slice = data_slice[data_slice["V (m^3/s)"] == V]
        data_slice = df
        T = data_slice["T (K)"].values
        x_methanol = data_slice[field].values
        plt.scatter(T, x_methanol)
    
        if i>Imax:
            break

    plt.xlabel('T')
    plt.ylabel(field)
    plt.savefig(f'{field_to_plot_fname}_{Imax}.png')
    plt.clf()