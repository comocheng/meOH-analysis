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
c = int(os.getenv('SLURM_ARRAY_TASK_ID'))

# Grabow model and RMG input files
cti_file_rmg = '../base/cantera/chem_annotated.cti'


Temps = [483.7-273.25, 499.3-273.25, 516.7-273.25]
Pressures = [15,30,50]
volume_flows = [0.00424,0.0106,0.02544]

settings  = list(itertools.product(Temps,
                                   Pressures,
                                   volume_flows,
                                  ))

# Mole fractions for run1
X_cos = 0.065
X_co2s = 0.261
X_h2s = 0.674

#######################################################################
# Input Parameters for combustor
#######################################################################

# RMG input files
cti_file_rmg = '../base/cantera/chem_annotated.cti'

#constants
pi = math.pi

# set initial temps, pressures, concentrations
temp = settings[array_i][0] # kelvin
pressure = settings[array_i][1]*ct.one_atm # Pascals
X_co = X_cos
X_co2= X_co2s
X_h2o = 0.8
X_h2 = X_h2s
mw_co = 28.01e-3  # [kg/mol]
mw_co2 = 44.01e-3 # [kg/mol]
mw_h2o = 18.01e-3 # [kg/mol]

# CO/CO2/H2/H2O: typical is 
concentrations_rmg = {'CO(3)': X_co,'CO2(4)': X_co2, 'H2O(5)':X_h2o}

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
cat_area = site_density/(cat_weight*cat_site_per_wt) #[m^3]

# reactor initialization
r = ct.IdealGasReactor(gas, energy='on')
rsurf = ct.ReactorSurface(surf, r, A=cat_area)
r.volume = rvol

# flow controllers (Graaf measured flow at 293.15 and 1 atm)
volume_flow = settings[array_i][2]                            # [m^3/s]
molar_flow = volume_flow*pressure/(8.3145*temp)               # [mol/s]
mass_flow = molar_flow*(X_co*mw_co+X_co2*mw_co2+X_h2o*mw_h2o) # [kg/s]
mfc = ct.MassFlowController(inlet, r, mdot=mass_flow)

# initialize reactor network
sim = ct.ReactorNet([r])
# set relative and absolute tolerances on the simulation
sim.rtol = 1.0e-9
sim.atol = 1.0e-21


#################################################
# Run single reactor 
#################################################

output_filename = f'Grabow_Results_RMG'
outfile = open(output_filename,'w')
writer = csv.writer(outfile)
writer.writerow(['T (C)', 'P (atm)', 'V (M^3/s)', 'X_co initial','X_co initial','X_co initial', 'T (C) final'] +
                gas.species_names + surf.species_names,)

t = 0.0
dt = 0.1

# run the simulation
sim.advance_to_steady_state(atol=1e-8)
writer.writerow([temp, pressure, volume_flow, 'X_co initial','X_co initial','X_co initial', gas.T] +
                list(gas.X) + list(surf.X),)