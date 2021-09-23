import numpy as np
from min_sbr import MinSBR

rmg_model_folder = "/home/moon/methanol/perturb_5000/run_0000/"
cti_file_path = "/home/moon/methanol/perturb_5000/run_0000/cantera/chem_annotated.cti"


# initialize reactor
sbr_ss = MinSBR(
    cti_file_path,
    rmg_model_folder,
    temperatures=np.linspace(400.0, 700.0, 20),
    pressures=[75],
    volumes=[4.24e-6],
    H2_fractions=[0.75],
    CO2_fractions=[0.5],
    rtol=1.0e-11,
    atol=1.0e-22,
    reactor_type=0,
    energy="off",
    reactime=1e5,
)

results = sbr_ss.run_reactor_ss_memory()
print("I made a minimal SBR. Now what's in it?")

# run to SS

