
from surface_reactor_test import sensitivity_test

#######################################################################
# Input Parameters for combustor
#######################################################################

# use example yaml
yaml_file = "methane_pox_on_pt.yaml"

# Reactor settings arrays for run
Temps = [400, 500, 600]
Pressures = [15, 30, 50, 75]
volume_flows = [0.00424, 0.0106, 0.02544]

# Mole fractions from runs for reference (min max and mid from Graaf data)
# X_cos = [0.053, 0.1365, 0.220]
# X_co2s = [0.261,0.1205, 0.261]
# X_h2s = [0.625, 0.7625, 0.900]

# CO+CO2/H2
# H2_fraction = [0.8,0.5,0.95,0.75]
H2_fraction = [0.8, 0.5, 0.95, 0.75]

# CO2/CO
CO_CO2_ratio = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

# reaction time
reactime = 1e3

run_reactor(
    yaml_file=yaml_file,
)