###############################################
# Grabow reactor batch script
# Chris Blais
# Northeastern University
# runs through all reactor conditions
###############################################
import pandas as pd
import numpy as np
import time
import cantera as ct
from matplotlib import pyplot as plt
import csv
import math
import os
import sys
import re
import itertools
import logging
from collections import defaultdict

# fet 
array_i = int(os.getenv('SLURM_ARRAY_TASK_ID'))

# Grabow model and RMG input files
cti_file_rmg = 'chem_annotated.cti'

# Reactor settings for run
# Temps = [483.7-273.25, 499.3-273.25, 516.7-273.25]
# Pressures = [15,30,50]
# volume_flows = [0.00424,0.0106,0.02544]

# Temps = [483.7,499.3,516.7]
Temps = [800]
# Pressures = [15,30,50]
Pressures = [75]
# volume_flows = [0.00424,0.0106,0.02544]
volume_flows = [0.00424]

# CO+CO2/H2
H2_fraction = [0.8,0.5,0.95,0.75]

#CO2/CO
CO_CO2_ratio = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9] 
settings  = list(itertools.product(Temps,
                                   Pressures,
                                   volume_flows,
                                   H2_fraction,
                                   CO_CO2_ratio
                                  ))



# Mole fractions for run (min max and mid from Graaf data)
# X_cos = [0.053, 0.1365, 0.220] 
# X_co2s = [0.261,0.1205, 0.261]
# X_h2s = [0.625, 0.7625, 0.900]


#######################################################################
# Input Parameters for combustor
#######################################################################

#constants
pi = math.pi

# set initial temps, pressures, concentrations
temp = settings[array_i][0] # kelvin
pressure = settings[array_i][1]*ct.one_atm # Pascals


X_h2 = settings[array_i][3]

if X_h2 == 0.75:
    X_h2o = 0.05
else:
    X_h2o = 0
    
X_co = (1-(X_h2+X_h2o))*(settings[array_i][4])
X_co2 = (1-(X_h2+X_h2o))*(1-settings[array_i][4])


#normalize mole fractions just in case
# X_co = X_co/(X_co+X_co2+X_h2)
# X_co2= X_co2/(X_co+X_co2+X_h2)
# X_h2 = X_h2/(X_co+X_co2+X_h2)

mw_co = 28.01e-3       # [kg/mol]
mw_co2 = 44.01e-3      # [kg/mol]
mw_h2 = 2.016e-3       # [kg/mol]
mw_h2o = 18.01528e-3   # [kg/mol]

co2_ratio = X_co2/(X_co + X_co2)
h2_ratio = (X_co2+X_co)/X_h2

# CO/CO2/H2/H2: typical is 
concentrations_rmg = {'CO(3)': X_co,'CO2(4)': X_co2, 'H2(2)':X_h2}

# initialize cantera gas and surface
gas= ct.Solution(cti_file_rmg,'gas')

# surf_grab = ct.Interface(cti_file,'surface1_grab', [gas_grab])
surf = ct.Interface(cti_file_rmg,'surface1', [gas])

# gas_grab.TPX = 
gas.TPX = temp, pressure, concentrations_rmg
surf.TP = temp, pressure

# create gas inlet
inlet = ct.Reservoir(gas)

#create gas outlet
exhaust = ct.Reservoir(gas)

# Reactor volume
rradius = 35e-3
rlength = 70e-3
rvol = (rradius**2)*pi*rlength

# Catalyst Surface Area
site_density = surf.site_density*1000                #[mol/m^2]cantera uses kmol/m^2, convert to mol/m^2
cat_weight = 4.24e-3                                 #[kg]
cat_site_per_wt = (300*1e-6)*1000                    #[mol/kg] 1e-6mol/micromole, 1000g/kg
cat_area = site_density/(cat_weight*cat_site_per_wt) #[m^2]

# reactor initialization
r = ct.IdealGasReactor(gas, energy='off')
rsurf = ct.ReactorSurface(surf, r, A=cat_area)
r.volume = rvol
surf.coverages = 'X(1):1.0'

# flow controllers (Graaf measured flow at 293.15 and 1 atm)
one_atm = ct.one_atm
FC_temp = 293.15
volume_flow = settings[array_i][2]                            # [m^3/s]
molar_flow = volume_flow*one_atm/(8.3145*293.15)               # [mol/s]
mass_flow = molar_flow*(X_co*mw_co+X_co2*mw_co2+X_h2*mw_h2+X_h2o*mw_h2o)   # [kg/s]
mfc = ct.MassFlowController(inlet, r, mdot=mass_flow)

# initialize reactor network
sim = ct.ReactorNet([r])
# set relative and absolute tolerances on the simulation
sim.rtol = 1.0e-11
sim.atol = 1.0e-22


# for i in range(0, surf.n_reactions):
#      if 'H2X(50)' in str(surf.reaction_equation(i)):
#         print(surf.reaction_equation(i), i)
#         surf.set_multiplier(0.0,i)

#################################################
# Run single reactor 
#################################################

output_filename = f'Grabow_Results_RMG_H2_{X_h2}_COCO2_{X_co}_{X_co2}_area{cat_area}.csv'
outfile = open(output_filename,'w')
writer = csv.writer(outfile)
writer.writerow(['T (C)', 'P (atm)', 'V (M^3/s)', 'X_co initial','X_co2 initial','X_h2 initial','X_h2o initial',
                 'CO2/(CO2+CO)','(CO+CO2/H2)', 'T (C) final'] + gas.species_names + surf.species_names,)

t = 0.0
dt = 0.01

# run the simulation

while t < 1000000.0:
    t += dt
    sim.advance(t)
   


writer.writerow([temp, pressure, volume_flow, X_co, X_co2, X_h2, X_h2o, co2_ratio, h2_ratio, gas.T] +
                 list(gas.X) + list(surf.X),)