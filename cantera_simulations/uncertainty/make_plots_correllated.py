from min_sbr import MinSBR
import os
import pandas as pd
import plotly.express as px
import platform
import time
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize 
from matplotlib import cm
import seaborn as sns


# get perturbations used for correllated run (make sure seed and dims are the same)
from torch.quasirandom import SobolEngine
# Define the uncertainty ranges based on the paper
DELTA_ALPHA_MAX = 0.15
# DELTA_E0_MAX = 30  # 3 eV is about 30 kJ/mol
DELTA_E0_MAX_J_MOL = 30000  # 3 eV is about 30000 J/mol

# only perturb vdw by +/- 0.2 eV
DELTA_E0_MAX_J_MOL_VDW = 20000

# Define the number of perturbations to run
N = 5000

# Create the pseudo randoms
sobol = SobolEngine(dimension=80, scramble=True, seed=100)
x_sobol = sobol.draw(N)


# Make list of lists for sobol perturbations
delta_pt_c = []
delta_pt_o = []
delta_pt_h = []
delta_pt_vdw = []

for row in x_sobol:
    delta_pt_c.append(float(DELTA_E0_MAX_J_MOL - 2.0 * row[0] * DELTA_E0_MAX_J_MOL)/9.6e4)
    delta_pt_o.append(float(DELTA_E0_MAX_J_MOL - 2.0 * row[1] * DELTA_E0_MAX_J_MOL)/9.6e4)
    delta_pt_h.append(float(DELTA_E0_MAX_J_MOL - 2.0 * row[2] * DELTA_E0_MAX_J_MOL)/9.6e4)
    delta_pt_vdw.append(float(DELTA_E0_MAX_J_MOL_VDW  - 2.0 * row[3] * DELTA_E0_MAX_J_MOL_VDW)/9.6e4)

index=0   


rmg_runs_dir = "/scratch/westgroup/methanol/perturb_5000_correllated/"

csv_files = glob.glob(os.path.join(rmg_runs_dir, 'run_0***', 'cantera', 'ct_analysis.csv'))

P = 3039750.0000000005
V = 3.32416e-05
T = []
x_methanol = []
perturb_c = []
perturb_o = []
perturb_h = []
perturb_vdw = []


# list of species from model
# 'N2', 'Ne',
# 'H2(2)', 'CO(3)', 'CO2(4)', 'H2O(5)', 'CH2O(6)', 'HCOOH(7)', 'CH3OH(8)',
# 'HCOOCH3(9)', 'CH4(24)', 'X(1)', 'H*(10)', 'O*(11)', 'OH*(12)',
# 'H2O*(13)', 'CO*(14)', 'CO2*(15)', 'HCO*(16)', 'HCOO*(17)', 'COOH*(18)',
# 'HCOOH*(19)', 'CH2O*(20)', 'CH3O*(21)', 'CH3O2*(22)', 'CH3OH*(23)',
# 'H2X(25)', 'COX2(26)', 'SX(56)(27)', 'CHOX2(28)', 'CH3OX(29)',
# 'CH2OX2(30)', 'SX(31)', 'SX(32)', 'SX(33)', 'C2HO2X(34)', 'CH4X(35)',
# 'CH3X(36)', 'SX(37)', 'SX(62)(38)', 'SX(58)(39)', 'CHO2X2(40)',
# 'C2H4OX(41)', 'SX(42)', 'C2H5OX(43)', 'SX(44)', 'C2H5OX(45)',

# set to 5000 to do all runs (does take some time)
Imax = 1000

# currently just plot one species
field = "CH3OH(8)"

# get points from all CSV files 
for i, csv_file in enumerate(csv_files):
    
    df = pd.read_csv(csv_file)
    # get run id from filename for referencing in sobol sequence
    match = re.compile('run_....')
    field_to_plot_fname = match.search(csv_file)
    run_id = int(field_to_plot_fname.group().replace("run_",""))
    
    if run_id < Imax:
        value = 500
        index = abs(df['T (K)'] - value).idxmin()

        
        T.append(df['T (K)'][index])
        x_methanol.append(df[field][index])
        
        
        perturb_c.append(delta_pt_c[run_id])
        perturb_o.append(delta_pt_o[run_id])
        perturb_h.append(delta_pt_h[run_id])
        perturb_vdw.append(delta_pt_vdw[run_id])


    if len(T) != len(perturb_c):
        print("failed, vectors are uneven")
        break

        
# get meoh TOF
Pres = 3039750.0 # Pa
Temp = value # K
Vflow = 3.32416e-5 # m^3/s
R = 8.3145
sites = 1.307404e-3
F = ((Pres*Vflow)/(R*Temp*sites))
x_methanol_TOF = [s*F for s in x_methanol]

# make dot plot 
fig, ax = plt.subplots()
cax=ax.scatter(perturb_vdw,perturb_c,c=x_methanol_TOF, s=50)
plt.xlabel('Pt_c')
plt.ylabel('Pt_vdw')
norm = Normalize(vmin = np.min(x_methanol), vmax = np.max(x_methanol))
cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
cbar.ax.set_ylabel('max meoh')
plt.savefig(f'plots/{field}_dots.png')

# make an interpolated contour plot
from scipy.interpolate import griddata

x = perturb_c
y = perturb_o
z = x_methanol_TOF

# ny, nx = 100, 200

xmin = min(x)
xmax = max(x)

ymin = min(y)
ymax = max(y)

# x = np.r_[x,xmin,xmax]
# y = np.r_[y,ymax,ymin]
# z = np.r_[z,z[0],z[-1]]
# xi = np.linspace(xmin, xmax, nx)
# yi = np.linspace(ymin, ymax, ny)

# xy = np.vstack((x, y)).T


def plot_contour(x,y,z,resolution = 100,contour_method='linear'):
    resolution = str(resolution)+'j'
    X,Y = np.mgrid[min(x):max(x):complex(resolution),   min(y):max(y):complex(resolution)]
    points = [[a,b] for a,b in zip(x,y)]
    Z = griddata(points, z, (X, Y), method=contour_method,fill_value=0.0)
    return X,Y,Z

X,Y,Z = plot_contour(x,y,z,resolution = 100,contour_method='linear')

fig, ax = plt.subplots()
ax.contourf(X,Y,Z)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

norm = Normalize(vmin = np.min(z), vmax = np.max(z))
cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
cbar.ax.set_ylabel('MeOH TOF (1/s)',rotation=270,labelpad=15)

ax.set_xlabel("Carbon Perturbation (eV)")
ax.set_ylabel("Oxygen Perturbation (eV)")

plt.savefig(f'plots/{field}_contour.png')







