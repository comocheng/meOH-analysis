import os
import sys
import numpy as np
import pandas as pd
import itertools
from multiprocessing import Pool
from min_sbr import MinSBR
# import time


def load_graaf_data():
    """
    Julia Treese
    Edited part to get Graaf conditions
    get Graaf conditions into a list of lists to run
    """
    file_name_overall = "/scratch/westgroup/methanol/meOH-analysis/cantera_simulations/Graaf_data/combined_experimental_runs.xlsx"
    file_name_feed1 = "/scratch/westgroup/methanol/meOH-analysis/cantera_simulations/Graaf_data/Feed_1.xlsx"
    file_name_feed2 = "/scratch/westgroup/methanol/meOH-analysis/cantera_simulations/Graaf_data/Feed_2.xlsx"
    file_name_feed3 = "/scratch/westgroup/methanol/meOH-analysis/cantera_simulations/Graaf_data/Feed_3.xlsx"
    file_name_feed4 = "/scratch/westgroup/methanol/meOH-analysis/cantera_simulations/Graaf_data/Feed_4.xlsx"
    file_name_feed5 = "/scratch/westgroup/methanol/meOH-analysis/cantera_simulations/Graaf_data/Feed_5.xlsx"
    file_name_feed6a = "/scratch/westgroup/methanol/meOH-analysis/cantera_simulations/Graaf_data/Feed_6a.xlsx"
    file_name_feed6b = "/scratch/westgroup/methanol/meOH-analysis/cantera_simulations/Graaf_data/Feed_6b.xlsx"
    file_name_feed7a = "/scratch/westgroup/methanol/meOH-analysis/cantera_simulations/Graaf_data/Feed_7a.xlsx"
    file_name_feed7b = "/scratch/westgroup/methanol/meOH-analysis/cantera_simulations/Graaf_data/Feed_7b.xlsx"
    
    
    df_overall = pd.read_excel(file_name_overall, engine='openpyxl')
    # Needed: [T, P, V, YH2, YCO2, wcat, meoh TOF, water TOF] -- Create a list of lists
    # Should be columns 2 (T), 1 (P), 3 (V), 6 (YH2), 5 (YCO2) 6 (cat weight)
    # Each list is the conditions of one experimental Graaf run

    # List of dataframes with feed conditions
    # df_list = [df_1, df_2, df_3, df_4, df_5, df_6a, df_6b, df_7a, df_7b]
    df_list = [df_overall]
    # Loop through dataframes and create a list of conditions based on Graaf runs
    # Loop through each row in the dataframes and add that row's conditions to the list of lists

    settings = []

    for i in range(len(df_list)):
        df = df_list[i]
        for row in range(len(df)):
            if not np.isnan(df.iloc[row, df.columns.get_loc('T(K)')]):
                row_conditions = [df.iloc[row, df.columns.get_loc('T(K)')], 
                                df.iloc[row, df.columns.get_loc('p (bar)')]/1.01325, #convert to atm for consistency
                                df.iloc[row,df.columns.get_loc('10^6 * V (M^3/s)')]*1e-6, # convert from graaf format
                                df.iloc[row,df.columns.get_loc('feed Yh2')], 
                                df.iloc[row,df.columns.get_loc('feed Yco')],
                                df.iloc[row,df.columns.get_loc('wcat (g)')]*1e-3, # convert from g to kg
                                df.iloc[row,df.columns.get_loc('MeOH TOF (mol/site/s)')],
                                df.iloc[row,df.columns.get_loc('H2O TOF (mol/site/s)')],]        
                settings.append(row_conditions)
    return settings


if len(sys.argv) < 2:
    raise ValueError("Incorrect usage. Must pass the cantera model file as an argument to this analysis script")

if not os.path.exists(sys.argv[1]):
    raise OSError(f"Path to the cantera model file does not exist: {sys.argv[1]}")

# start = time.time()

# rmg_model_folder = "/home/moon/methanol/perturb_5000/run_0000/"
# rmg_model_folder = "/home/sevy/methanol/perturb_5000/run_0000/"
# rmg_model_folder = "/scratch/westgroup/methanol/perturb_5000/run_0000/"

# cti_file_path = "/home/moon/methanol/perturb_5000/run_0000/cantera/chem_annotated.cti"
# cti_file_path = "/home/sevy/methanol/perturb_5000/run_0000/cantera/chem_annotated.cti"
# cti_file_path = "/scratch/westgroup/methanol/perturb_5000/run_0000/cantera/chem_annotated.cti"

cti_file_path = sys.argv[1]
rmg_model_folder = os.path.dirname(cti_file_path)
csv_path = os.path.join(rmg_model_folder, "ct_analysis_graaf.csv")


# pressures = [75.0]
volume_flows = [3.32416e-5] # updated to duplicate grabow's space velocity of 7.84e-3 m^3/kg/s


# generate settings array
settings = load_graaf_data()


def run_reactor(setting):

    # initialize reactor
    sbr_ss = MinSBR(
        cti_file_path,
        rmg_model_folder,
        temperature=setting[0],
        pressure=setting[1],
        volume_flow=setting[2],
        x_H2=setting[3],
        x_CO2=(1-(setting[3]+setting[4])),
        x_CO=setting[4],
        catalyst_weight=setting[5],
        rtol=1.0e-11,
        atol=1.0e-22,
        reactor_type=1,
        energy="off",
        reactime=1e5,
        meoh_tof=setting[6],
        h2o_tof=setting[7],
    )

    results = sbr_ss.run_reactor_ss_memory()
    return results


# Too much memory? is that why it's slow?
with Pool() as p:
    result = p.map(run_reactor, settings)

df = pd.DataFrame(result)
df.to_csv(csv_path)

# end = time.time()
# print(f"Completed {len(settings)} processes in {end-start} seconds")
