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


def save_flux_diagrams(*phases, suffix=''):
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
    """
    Saves the flux diagrams. The filenames have a suffix if provided,
    so you can keep them separate and not over-write.
    """
    for element in 'CHONX':
        for phase_object in phases:
            phase = phase_object.name

            diagram = ct.ReactionPathDiagram(phase_object, element)
            diagram.title = f'Reaction path diagram following {element} in {phase}'
            diagram.label_threshold = 0.01

            dot_file = f"{suffix}/reaction_path_{element}_{phase}.dot"
            img_file = f"{suffix}/reaction_path_{element}_{phase}.png"
            dot_bin_path = '/Users/blais.ch/anaconda3/pkgs/graphviz-2.40.1-hefbbd9a_2/bin/dot'
            img_path = os.path.join(os.getcwd(), img_file)
            diagram.write_dot(dot_file)
            #print(diagram.get_data())

            print(f"Wrote graphviz input file to '{os.path.join(os.getcwd(), dot_file)}'.")
            os.system(f'dot {dot_file} -Tpng -o{img_file} -Gdpi=200')
            print(f"Wrote graphviz output file to '{img_path}'.")


def run_reactor(cti_file, t_array=[528], p_array=[75], v_array=[0.00424], h2_array=[0.75], co2_array=[0.5], 
                rtol=1.0e-11, atol=1.0e-22, reactor_type=0, energy='off'):
    
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
    try:
        array_i = int(os.getenv('SLURM_ARRAY_TASK_ID'))
    except TypeError:
        array_i = 0
    

    # this should probably be outside of function
    settings  = list(itertools.product(t_array,
                                p_array,
                                v_array,
                                h2_array,
                                co2_array
                                ))

    #constants
    pi = math.pi

    # set initial temps, pressures, concentrations
    temp = settings[array_i][0] # kelvin
    pressure = settings[array_i][1]*ct.one_atm # Pascals


    X_h2 = settings[array_i][3]

    # Per Grabow experiments, add in H2O for X=0.75 H2 run
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
    gas= ct.Solution(cti_file,'gas')

    # surf_grab = ct.Interface(cti_file,'surface1_grab', [gas_grab])
    surf = ct.Interface(cti_file,'surface1', [gas])

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
    if reactor_type == 0:
        r = ct.Reactor(gas, energy=energy)
        reactor_type_str = 'Reactor'
    elif reactor_type == 1:
        r = ct.IdealGasReactor(gas, energy=energy)
        reactor_type_str = 'IdealGasReactor'
    elif reactor_type == 2:
        r = ct.ConstPressureReactor(gas, energy=energy)
        reactor_type_str = 'ConstPressureReactor'
    elif reactor_type == 3:
        r = ct.IdealGasConstPressureReactor(gas, energy=energy)
        reactor_type_str = 'IdealGasConstPressureReactor'
    
    rsurf = ct.ReactorSurface(surf, r, A=cat_area)
    r.volume = rvol
    surf.coverages = 'X(1):1.0'

    # flow controllers (Graaf measured flow at 293.15 and 1 atm)
    one_atm = ct.one_atm
    FC_temp = 293.15
    volume_flow = settings[array_i][2]                            # [m^3/s]
    molar_flow = volume_flow*one_atm/(8.3145*FC_temp)               # [mol/s]
    mass_flow = molar_flow*(X_co*mw_co+X_co2*mw_co2+X_h2*mw_h2+X_h2o*mw_h2o)   # [kg/s]
    mfc = ct.MassFlowController(inlet, r, mdot=mass_flow)

    # initialize reactor network
    sim = ct.ReactorNet([r])
    # set relative and absolute tolerances on the simulation
    sim.rtol = 1.0e-11
    sim.atol = 1.0e-22


    #################################################
    # Run single reactor 
    #################################################

    # round numbers so they're easier to read
    # 
    # temp_str = '%s' % '%.3g' % tempn

    cat_area_str = '%s' % '%.3g' % cat_area
    results_path = os.path.dirname(os.path.abspath(__file__)) + '/results'
    flux_path = os.path.dirname(os.path.abspath(__file__)) + '/flux_diagrams'
    try:  
        os.mkdir(results_path)  
        os.mkdir(flux_path)
    except OSError as error:  
        print(error) 
    
    output_filename = results_path + f'/Spinning_basket_area_{cat_area_str}_{reactor_type_str}_energy_{energy}.csv'
    outfile = open(output_filename,'w')
    writer = csv.writer(outfile)
    writer.writerow(['T (C)', 'P (atm)', 'V (M^3/s)', 'X_co initial','X_co2 initial','X_h2 initial','X_h2o initial',
                    'CO2/(CO2+CO)','(CO+CO2/H2)', 'T (C) final', 'Rtol', 'Atol'] + gas.species_names + surf.species_names,)

    t = 0.0
    dt = 0.01

    # run the simulation

    while t < 1000.0:
    # while t < 1000000.0:
        t += dt
        sim.advance(t)
    
    writer.writerow([temp, pressure, volume_flow, X_co, X_co2, X_h2, X_h2o, co2_ratio, h2_ratio, gas.T, sim.rtol, sim.atol] +
                    list(gas.X) + list(surf.X),)
    outfile.close()

    # save flux diagrams
    save_flux_diagrams(gas,suffix=flux_path)
    save_flux_diagrams(surf,suffix=flux_path)
    return 


#######################################################################
# Input Parameters for combustor
#######################################################################


# filepath for writing files
cti_file = os.path.dirname(os.path.abspath(__file__)) +'/chem_annotated.cti'
print(cti_file)
# Reactor settings arrays for run
Temps = [483.7-273.25, 499.3-273.25, 516.7-273.25]
Pressures = [15,30,50,75]
volume_flows = [0.00424,0.0106,0.02544]

# Mole fractions for run (min max and mid from Graaf data)
# X_cos = [0.053, 0.1365, 0.220] 
# X_co2s = [0.261,0.1205, 0.261]
# X_h2s = [0.625, 0.7625, 0.900]

# CO+CO2/H2
H2_fraction = [0.8,0.5,0.95,0.75]

#CO2/CO
CO_CO2_ratio = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9] 

run_reactor(cti_file=cti_file)