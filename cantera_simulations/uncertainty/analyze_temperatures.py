import os
import sys
import numpy as np
import pandas as pd
import itertools
from multiprocessing import Pool
from min_sbr import MinSBR
# import time


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
csv_path = os.path.join(rmg_model_folder, "ct_analysis.csv")

temperatures = np.linspace(400.0, 700.0, 20)
pressures = np.linspace(30.0, 75.0, 1)
# pressures = [75.0]
volume_flows = [3.32416e-5] # updated to duplicate grabow's space velocity of 7.84e-3 m^3/kg/s


# generate settings array
settings = list(
    itertools.product(
        temperatures,
        pressures,
        volume_flows,
    )
)


def run_reactor(setting):

    # initialize reactor
    sbr_ss = MinSBR(
        cti_file_path,
        rmg_model_folder,
        temperature=setting[0],
        pressure=setting[1],
        volume_flow=setting[2],
        x_H2=0.5,
        x_CO2=0.25,
        x_CO=0.25,
        rtol=1.0e-11,
        atol=1.0e-22,
        reactor_type=1,
        energy="off",
        reactime=1e5,
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
