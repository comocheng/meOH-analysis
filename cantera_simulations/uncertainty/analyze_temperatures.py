import numpy as np
from min_sbr import MinSBR

rmg_model_folder = "/home/moon/methanol/perturb_5000/run_0000/"
cti_file_path = "/home/moon/methanol/perturb_5000/run_0000/cantera/chem_annotated.cti"


# initialize reactor
sbr_ss = MinSBR(
    cti_file_path,
    rmg_model_folder,
    temperature=400,
    pressure=75,
    volume_flow=4.24e-6,
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
print("I made a minimal SBR. Now what's in it?")

# run to SS
# # generate settings array. - this goes outside the init
#         self.settings = list(
#             itertools.product(
#                 self.temperatures,
#                 self.pressures,
#                 self.volumes,
#                 self.H2_fractions,
#                 self.CO2_fractions,
#                 self.catalyst_weights,
#             )
#         )
