# Compute the local uncertainty for methanol formation on copper
#
# Use this script with vscode to debug, then put result in the Jupyter Notebook

import os

from IPython.display import display, Image

from rmgpy.tools.uncertainty import Uncertainty, process_local_results
from rmgpy.tools.canteramodel import get_rmg_species_from_user_species
from rmgpy.species import Species

# Need to find a way to combine surface and gas chemkin files

# Must use annotated chemkin file
rmg_repo = '/home/moon/methanol/meOH-synthesis'
chemkin_file = os.path.join(rmg_repo, 'base/chemkin/chem_annotated-surface.inp')
dict_file = os.path.join(rmg_repo, 'base/chemkin/species_dictionary.txt')

# Initialize the Uncertainty class instance and load the model
uncertainty = Uncertainty(output_directory='./temp/uncertainty')
uncertainty.load_model(chemkin_file, dict_file)

# Map the species to the objects within the Uncertainty class
methanol = Species().from_smiles('C[OH]')

h_ads = Species().from_adjacency_list("""
1 H u0 p0 c0 {2,S}
2 X u0 p0 c0 {1,S}
""")

co_ads = Species().from_adjacency_list("""
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,D}
3 X u0 p0 c0 {2,D}
""")

co2_ads = Species().from_adjacency_list("""
1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
4 X u0 p0 c0
""")

co2_ads = Species().from_adjacency_list("""
1 O u0 p2 c0 {3,D}
2 O u0 p2 c0 {3,D}
3 C u0 p0 c0 {1,D} {2,D}
4 X u0 p0 c0
""")

methanol_ads = Species().from_adjacency_list("""
CH3OH*(23)
1 O u0 p2 c0 {2,S} {6,S}
2 C u0 p0 c0 {1,S} {3,S} {4,S} {5,S}
3 H u0 p0 c0 {2,S}
4 H u0 p0 c0 {2,S}
5 H u0 p0 c0 {2,S}
6 H u0 p0 c0 {1,S}
7 X u0 p0 c0
""")

# also try formate, other species from Grabow 2011
mapping = get_rmg_species_from_user_species([methanol_ads, co2_ads, co_ads, h_ads], uncertainty.species_list)

# Define the reaction conditions
# pick a H2, CO, CO2 concentration
initial_mole_fractions = {mapping[h_ads]: 0.8, mapping[co2_ads]: 0.18, mapping[co2_ads]: 0.02}
T = (400, 'K')
P = (15, 'bar')
termination_time = (0.5, 'ms')
sensitive_species = [mapping[methanol_ads]]

####################################################
# Run sensitivity analysis
# Perform the sensitivity analysis
uncertainty.sensitivity_analysis(initial_mole_fractions, sensitive_species, T, P, termination_time, number=5, fileformat='.png')

# Show the sensitivity plots
for species in sensitive_species:
    print('{}: Reaction Sensitivities'.format(species))
    index = species.index
    display(Image(filename=os.path.join(uncertainty.output_directory, 'solver', 'sensitivity_1_SPC_{}_reactions.png'.format(index))))

    print('{}: Thermo Sensitivities'.format(species))
    display(Image(filename=os.path.join(uncertainty.output_directory, 'solver', 'sensitivity_1_SPC_{}_thermo.png'.format(index))))
