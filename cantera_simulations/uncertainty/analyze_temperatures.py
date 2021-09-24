import numpy as np
import itertools
from multiprocessing import Pool
from min_sbr import MinSBR
import time


N_proc = 0

start = time.time()

# rmg_model_folder = "/home/moon/methanol/perturb_5000/run_0000/"
# rmg_model_folder = "/home/sevy/methanol/perturb_5000/run_0000/"
rmg_model_folder = "/scratch/westgroup/methanol/perturb_5000/run_0000/"

# cti_file_path = "/home/moon/methanol/perturb_5000/run_0000/cantera/chem_annotated.cti"
# cti_file_path = "/home/sevy/methanol/perturb_5000/run_0000/cantera/chem_annotated.cti"
cti_file_path = "/scratch/westgroup/methanol/perturb_5000/run_0000/cantera/chem_annotated.cti"


temperatures = np.linspace(400.0, 700.0, 40)
pressures = np.linspace(30.0, 75.0, 5)
# pressures = [75.0]
volume_flows = [4.24e-6]


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
        x_H2=0.75,
        x_CO2=0.25,
        x_CO=0,
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

# results = run_reactor(settings[0])

end = time.time()
print(f"Completed {len(settings)} processes in {end-start} seconds using {N_proc} workers")

# completed 40 processes in 309.3567500114441 seconds
