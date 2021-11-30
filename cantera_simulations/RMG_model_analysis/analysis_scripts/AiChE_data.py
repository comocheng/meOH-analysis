import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn
import time
from pathlib import Path
import os
import csv
import glob
from scipy.stats import gaussian_kde

import scipy.io
import math



###################################
# load in Graaf experimental data
###################################
# exclude feed 8 because it is a monolith reactor
path_str = "../../Graaf_data/combined_experimental_runs.xlsx"
df_graaf = pd.read_excel(path_str , engine='openpyxl')

###################################
# load in Grabow model data
###################################
# load in initial conditions from experiment from matlab
ml_results_file = "matlab_grabow_results.csv"
translation = {39: None, 91: None, 93: None} # remove ', [, and ]
count = 0
first_row = True

for path in Path("../../Grabow_matlab_data/original_runs").rglob('*.mat'):
 
    path_str = str(path)
    
    mat = scipy.io.loadmat(path_str)
    conditions = mat['ccondition']
    reactions = mat['reaction']
    
    # get relevant TOFs, catalyst weights, # sites
    cat_weight = float(conditions["CatalystWeight"])
    num_sites = float(conditions["Sites"])
    flow_site_basis = float(conditions["Fin"]) # molar flow rate divided by moles sites
    flow_cm_3_min = float(conditions["Flow"]) 
    P_atm = float(conditions["P"]) 
    T_grabow_K = float(conditions["T"])
    inlet_CO = float(conditions["molfrac"][0][0][0][0][1][0][0])
    inlet_CO2 = float(conditions["molfrac"][0][0][0][0][1][0][1])
    inlet_h2 = float(conditions["molfrac"][0][0][0][0][1][0][2])
    
    #read in species names to use as columns in dataframe
    species_names = []

    for spec in range(len(mat['species'][0,:])):
        species_string = str(mat['species'][0,spec][0])
        species_string = species_string.translate(translation)
        species_names.append(species_string)

    # read in results (surface coverages, partial pressures)
    results = mat['Y']

    if first_row:
        # create data frame with species names as column headers
        df_grabow = pd.DataFrame(data=results,columns=species_names)
        df_grabow = df_grabow.tail(1)
        df_grabow["Catalyst Weight (g)"] = cat_weight
        df_grabow["Number of Sites (mol)"] = num_sites 
        df_grabow["Fin site basis (mol/site/s)"] = flow_site_basis
        df_grabow["flow (cm^3/min)"] = flow_cm_3_min
        df_grabow["P (atm)"] = P_atm
        df_grabow["T grab (K)"] = T_grabow_K
        df_grabow["inlet mol frac CO"] = inlet_CO
        df_grabow["inlet mol frac CO2"] = inlet_CO2
        df_grabow["inlet mol frac h2"] = inlet_h2
        
        
        # get total pressure by adding each species partial pressure
        total_pressure = 0
        for column in df_grabow:
            if "g" in column[-1].strip():
                total_pressure += float(df_grabow[column])
        df_grabow["total pressure (bar)"] = total_pressure
        
        first_row = False
    else: 
        # add tail from new dataframe
        new_df = pd.DataFrame(data=results,columns=species_names)
        new_df = new_df.tail(1)
        new_df["Catalyst Weight (g)"] = cat_weight
        new_df["Number of Sites (mol)"] = num_sites 
        new_df["Fin site basis (mol/site/s)"] = flow_site_basis
        new_df["flow (cm^3/min)"] = flow_cm_3_min
        new_df["P (atm)"] = P_atm
        new_df["T grab (K)"] = T_grabow_K
        new_df["inlet mol frac CO"] = inlet_CO
        new_df["inlet mol frac CO2"] = inlet_CO2
        new_df["inlet mol frac h2"] = inlet_h2
        
        # get total pressure by adding each species partial pressure
        total_pressure = 0
        for column in new_df:
            if "g" in column[-1].strip():
                total_pressure += float(new_df[column])
                
        new_df["total pressure (bar)"] = total_pressure
        
        df_grabow = df_grabow.append(new_df, ignore_index=True)


# convert all partial pressures to mole fractions 
# any value with a "g" (gas) after the species name
for column in df_grabow:
    if "g" in column[-1].strip():
        df_grabow[column] = df_grabow[column]/df_grabow["total pressure (bar)"]

# calculate grabow turn over frequencies
df_grabow["Grabow MeOH TOF (1/s)"] = df_grabow["CH3OHg"]*df_grabow["Fin site basis (mol/site/s)"]
df_grabow["H2O TOF (1/s)"] = df_grabow["H2Og"]*df_grabow["Fin site basis (mol/site/s)"]


# match grabow runs to their corresponding Graaf run. 
count = 0
found = False
for indexgrabow, row_grabow in df_grabow.iterrows():
    for index_graaf, row_graaf in df_graaf.iterrows():
        
        # get the values in the Graaf experiment that match the grabow values to 1e-3
        if (math.isclose(row_graaf["10^6 * V (M^3/s)"], row_grabow["flow (cm^3/min)"]/60,rel_tol=1e-3) and
            math.isclose(row_graaf["T(K)"], row_grabow["T grab (K)"],rel_tol=1e-3) and
            math.isclose(row_graaf["p (bar)"], row_grabow["P (atm)"]*1.01325,rel_tol=1e-3) and
            math.isclose(row_graaf["feed Yco"], row_grabow["inlet mol frac CO"],rel_tol=1e-3) and
            math.isclose(row_graaf["feed Yco2"], row_grabow["inlet mol frac CO2"],rel_tol=1e-3) and
            math.isclose(row_graaf["feed Yh2"], row_grabow["inlet mol frac h2"],rel_tol=1e-3)):
            
            new_series_data = row_grabow.append(row_graaf)
            found = True
    if found:
        if count == 0:
            df_combined = pd.DataFrame(new_series_data, index = new_series_data.index).transpose()
            count+=1
            
        if count > 0:
            df_combined = df_combined.append(pd.DataFrame(new_series_data).transpose())
            count+=1
            
        found = False
        
            
###################################
# load in RMG model data
###################################
rmg_runs_dir = "/scratch/westgroup/methanol/perturb_5000/"

# load uncertainty analysis cantera results
csv_files = glob.glob(os.path.join(rmg_runs_dir, 'run_****', 'cantera', 'ct_analysis_graaf.csv'))

# load list of csv files
csv_list = []

# get mse of model: 
mse_meoh = []
mse_h2o = []

# save all TOFs to 
RMG_meoh_tof = []
Graaf_meoh_tof = []

for i, csv_file in enumerate(csv_files):
    df_rmg = pd.read_csv(csv_file)
    count = 0
    found = False
    
    for indexrmg, row_rmg in df_rmg.iterrows():
        for index_graaf, row_graaf in df_combined.iterrows():
            if (math.isclose(row_graaf["10^6 * V (M^3/s)"], row_rmg["V (m^3/s)"]*1e6,rel_tol=1e-2) and
                math.isclose(row_graaf["T(K)"], row_rmg["T (K)"],rel_tol=1e-2) and
                math.isclose(row_graaf["p (bar)"], row_rmg["P (Pa)"]/1e5,rel_tol=1e-2) and
                math.isclose(row_graaf["feed Yco"], row_rmg["x_CO initial"],rel_tol=1e-2) and
                math.isclose(row_graaf["feed Yco2"], row_rmg["x_CO2 initial"],rel_tol=1e-2) and
                math.isclose(row_graaf["feed Yh2"], row_rmg["x_H2 initial"],rel_tol=1e-2)):

                new_series_data = row_graaf.append(row_rmg)
                found = True
        if found:
            if count == 0:
                df_combined_rmg = pd.DataFrame(new_series_data, index = new_series_data.index).transpose()
                count+=1

            if count > 0:
                df_combined_rmg = df_combined_rmg.append(pd.DataFrame(new_series_data).transpose())
                count+=1
                
                RMG_meoh_tof = np.append
                Graaf_meoh_tof = []
                
            found = False
            
        
    try:
        mse_meth = df_combined_rmg["error squared MeOH TOF"].sum()/len(df_combined_rmg["error squared MeOH TOF"])
        mse_wat = df_combined_rmg["error squared H2O TOF"].sum()/len(df_combined_rmg["error squared H2O TOF"])
        csv_list.append(csv_file)
        mse_meoh.append(mse_meth)
        mse_h2o.append(mse_wat)
        
        # append the tof values for density plot
        if i == 0: 
            RMG_meoh_tof = df_combined_rmg["RMG MeOH TOF 1/s"].values
            Graaf_meoh_tof = df_combined_rmg["MeOH TOF (mol/site/s)"].values
            RMG_h2o_tof = df_combined_rmg["RMG H2O TOF 1/s"].values
            Graaf_h2o_tof = df_combined_rmg["H2O TOF (mol/site/s)"].values
        else:
            RMG_meoh_tof = np.append(RMG_meoh_tof, df_combined_rmg["RMG MeOH TOF 1/s"].values)
            Graaf_meoh_tof = np.append(Graaf_meoh_tof, df_combined_rmg["MeOH TOF (mol/site/s)"].values)
            RMG_h2o_tof = np.append(RMG_h2o_tof, df_combined_rmg["RMG H2O TOF 1/s"].values)
            Graaf_h2o_tof = np.append(Graaf_h2o_tof, df_combined_rmg["H2O TOF (mol/site/s)"].values )          
        
    except:
        print(f'no column in {csv_file}')

    

# get the model number that is correllated to 
meoh_csv = csv_list[mse_meoh.index(min(mse_meoh))]
h2o_csv = csv_list[mse_h2o.index(min(mse_h2o))]

###################################
# generate density plot for RMG TOFs
################################### 
# Calculate the MeOH TOF point density
xy = np.vstack([Graaf_meoh_tof,RMG_meoh_tof])
z = gaussian_kde(xy)(xy)        
fig, ax = plt.subplots()
ax.scatter(Graaf_meoh_tof, RMG_meoh_tof, c=z, s=300)

# plot the MeOH TOF density graph
plt.xlabel('Graaf TOF (1/s)')
plt.ylabel("RMG MeOH TOF (1/s)" + field)
plt.savefig('TOF_density_methanol.png')
plt.clf()

# Calculate the H2O TOF point density
xy = np.vstack([Graaf_h20_tof,RMG_h2o_tof])
z = gaussian_kde(xy)(xy)        
fig, ax = plt.subplots()
ax.scatter(Graaf_h20_tof,RMG_h2o_tof, c=z, s=300)

# plot the H2O TOF density graph
plt.xlabel('Graaf TOF (1/s)')
plt.ylabel("RMG H2O TOF (1/s)" + field)
plt.savefig('TOF_density_water.png')
plt.clf()


###################################
# generate density plot for RMG TOFs
################################### 
