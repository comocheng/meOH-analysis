###############################################################################
# make park model
###############################################################################
import cantera as ct
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt


import rmgpy
from rmgpy.data.thermo import ThermoDatabase
from rmgpy.data.kinetics import KineticsDatabase
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.reaction import Reaction
import inspect
import copy
from rmgpy.kinetics.surface import SurfaceArrhenius
from rmgpy.kinetics.surface import StickingCoefficient
from rmgpy.quantity import ScalarQuantity

import rmgpy.chemkin as Chemkin
from cantera import ck2cti

###############################################################################
# useful functions
###############################################################################
def get_thermo(spec_str):
    '''
    takes a string input and returns a species object with complete thermo
    this may already exist in RMG.
    '''
    spec = Species()
    spec.from_smiles(spec_str)
    est_thermo = thermo_database.get_thermo_data(spec,metal_to_scale_to="Cu111")
    spec.thermo = est_thermo
    return spec

def get_gas_phase_precurs(spec):
    ''' 
    adapted from ThermoDatabase method: 
    get_thermo_data_for_surface_species()
    gets a Species object corresponding to the gas phase precursor for 
    a given surface species
    does NOT apply adsorption correction!
    '''
    dummy_molecules = spec.molecule[0].get_desorbed_molecules()
    
    for mol in dummy_molecules:
        mol.clear_labeled_atoms()
    if len(dummy_molecules) == 0:
        raise RuntimeError(f"Cannot get thermo for gas-phase molecule")

    # if len(molecule) > 1, it will assume all resonance structures have already been 
    #generated when it tries to generate them, so evaluate each configuration separately 
    # and pick the lowest energy one by H298 value
    
    gas_phase_species_from_libraries = []
    gas_phase_species_estimates = []
    for dummy_molecule in dummy_molecules:
        dummy_species = Species()
        dummy_species.molecule = [dummy_molecule]
        dummy_species.generate_resonance_structures()
        dummy_species.thermo = thermo_database.get_thermo_data(dummy_species)
        if dummy_species.thermo.label:
            gas_phase_species_from_libraries.append(dummy_species)
        else:
            gas_phase_species_estimates.append(dummy_species)

    # define the comparison function to find the lowest energy
    def lowest_energy(species):
        if hasattr(species.thermo, 'H298'):
            print(species.thermo.H298.value_si)
            return species.thermo.H298.value_si
        else:
            print(species.thermo.get_enthalpy(298.0))
            return species.thermo.get_enthalpy(298.0)
        

    if gas_phase_species_from_libraries:
        species = min(gas_phase_species_from_libraries, key=lowest_energy)
    else:
        species = min(gas_phase_species_estimates, key=lowest_energy)

    thermo = species.thermo
    
    return species

def update_thermo(spec, name, be1, be2):
    '''
    updates species thermo given an input for binding energy.
    input species object (spec)
    park name as string (name)
    two floats for the original binding
    energy (be1) and the "correct" binding energy (be2)
    '''
    spec_new = copy.deepcopy(spec)
    ev_2_kj = 9.6e4
    be_diff = (be_dict[name] - be_dict_park[name])*9.6e4
    new_h298 = spec.thermo.H298.value_si - be_diff
    spec_new.thermo.H298.value_si = new_h298
    print(name, id(spec_new.thermo.H298.value_si), id(spec.thermo.H298.value_si))
    print(name, spec_new.thermo.H298.value_si, spec.thermo.H298.value_si, be_diff)
    return spec_new

def make_reaction(reactants, products, rxn_str, A, Ea, stick = False,):
    '''
    make a rmgpy reaction object. 
    takes a list of the species objects for 
    the reactants and products.
    takes a string for the reaction string
    if Stick is true, A-factor is the sticking coefficient
    '''
    if stick: 
        kinetics = StickingCoefficient(
            A=A, 
            n=0.0, 
            Ea=Ea, 
            T0=(1.0, "K"), 
            Tmin=None, 
            Tmax=None, 
            Pmin=None, 
            Pmax=None,
            coverage_dependence=None, 
            comment=''
        )
    else: 
        kinetics = SurfaceArrhenius(
            A=A, 
            n=0.0, 
            Ea=Ea, 
            T0=(1.0, "K"), 
            Tmin=None, 
            Tmax=None, 
            Pmin=None, 
            Pmax=None,
            coverage_dependence=None, 
            comment=''
        ) 
        
    
    # use the rmgpy reaction object 
    rxn = Reaction(
        index=-1,
        label=rxn_str,
        reactants=reactants,
        products=products,
        specific_collider=None,
        kinetics=kinetics,
        network_kinetics=None,
        reversible=True,
        transition_state=None,
        duplicate=False,
        degeneracy=1,
        pairs=None,
        allow_pdep_route=False,
        elementary_high_p=False,
        allow_max_rate_violation=False,
        rank=None,
        comment='',
        is_forward=None,
        )
    
    return rxn

def convert_to_nasa(spec):
    thermo = spec.thermo
    thermo_nasa = thermo.to_nasa(298, 1500, 1000)
    spec.thermo = thermo_nasa

###############################################################################
# initiialize things
###############################################################################

# quick check that we are using the correct rmgpy and version

print('using rmgpy at: ',inspect.getfile(rmgpy))
print('using rmgpy version: ', rmgpy.__version__)

# save rmgpy and db directory. db is assumed to be in the same
# folder as RMG-Py
rmg_py_path = inspect.getfile(rmgpy).split("rmgpy")[0]
rmg_db_path = rmg_py_path.split("RMG-Py")[0] + "RMG-database/"

# import data 
# set absolute location, using './' in jupyter performs differently
# in vscode
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
park_xl_file =os.path.join(__location__,'park_thermo_and_rates.xlsx')
BE_sheet='Binding Energies'
rxn_sheet = 'reactions'
be_df = pd.read_excel(park_xl_file, sheet_name=BE_sheet, engine='openpyxl')
rxn_df = pd.read_excel(park_xl_file, sheet_name=rxn_sheet, engine='openpyxl')

# output files 
chemkin_gas_file = os.path.join(__location__, 'park_gas.inp')
chemkin_surface_file = os.path.join(__location__ + '/park_surf.inp') # why do we need a / for surface?
cantera_file = os.path.join(__location__,'park_mech.cti')

###############################################################################
# Constants/values
###############################################################################
site_density_mol_cm = 2.943e-09
site_density_si = site_density_mol_cm * 1e4
site_density_object = ScalarQuantity(site_density_si, 'mol/m^2')

###############################################################################
# get thermo for all species in RMG model. adjust BEs per the sheet values
###############################################################################
db_input_path = rmg_db_path + 'input/'

# load the thermo database
library_path = db_input_path + 'thermo/'
thermo_libraries = [
    'surfaceThermoPt111',
    
]
thermo_database = ThermoDatabase()
thermo_database.load(
    library_path,
    libraries=thermo_libraries,
    depository=False,
    surface=True
    )

# load the kinetics database
kin_libraries_dir = db_input_path + "kinetics/libraries/Surface/"
kin_fam_dir = db_input_path + "kinetics/families/"

kinetics_libraries = [
    'CPOX_Pt/Deutschmann2006_adjusted',   
]

kinetics_families = ['surface']

kinetics_database = KineticsDatabase()
kinetics_database.load_recommended_families(kin_fam_dir  + 'recommended.py')
kinetics_database.load_families(
    path=kin_fam_dir,
    families=kinetics_families,
)

kinetics_database.load_libraries(
    kin_libraries_dir,
    libraries=kinetics_libraries
)

# get binding energies 
# need a dictionary translating species names to smiles
# need a dictionary translating species names to smiles
spec_smiles_dict = {
    'CO*':'O=C=[*]', 
    'CO2*':'O=C=O.[*]', 
    'H*':'[H]*',
    'H2O*':'O.[*]',
    'CH3OH*':'CO.[*]',
    'O*':'O=[*]',
    'OH*':'O[*]',
    'HCO*':'O=C*',
#     'HCOO**':'O=CO[*][*]', #formate, bidentate
    'HCOO**':'O=CO[*].[*]', # formate, bidentate, plus extra X
    'H2CO2*':'[*]OCO[*]',
    'COOH*':'O=C(O)[*]',
    'CH2O*':'C=O.[*]',
    'CH3O*':'CO[*]',
    'CH3O2*':'OCO[*]',
    '*':'[*]',
}

# also need a dict of gas phase species to get be's from
# key is surface species, value is Gas phase precursor 
# either from RMGs estimate or if it's explicitly known, 
# just the gas phase version (e.g. 'CO2*': 'CO2') 
gas_pre_dict = {
    'CO*':'[C-]#[O+]', 
    'CO2*':'O=C=O', 
    'H*':'[H]',
    'H2O*':'O',
    'CH3OH*':'CO',
    'O*':'[O]',
    'OH*':'[OH]',
    'HCO*':'[CH]=O',
    'HCOO**':'[O]C=O', #formate, bidentate
    'H2CO2*':'[O]C[O]',
    'COOH*':'O=[C]O',
    'CH2O*':'C=O',
    'CH3O*':'C[O]',
    'CH3O2*':'[O]CO',
    '*':'[*]',
}

# all of the gas phase species in the model
gas_smiles_dict = {
    'CO':'[C-]#[O+]', 
    'CO2':'O=C=O', 
    'H2O':'O',
    'CH3OH':'CO',
    'CH2O':'C=O',
    'H2':'[H][H]',
}

# construct a dictionary of binding energies
be_dict = {}
for label in spec_smiles_dict.keys():
    surf_spec = get_thermo(spec_smiles_dict[label])
    gas_spec = get_thermo(gas_pre_dict[label])
    
    surf_h298 = surf_spec.thermo.get_enthalpy(298)
    gas_h298 = gas_spec.thermo.get_enthalpy(298)
    
    be_dict[label] = (surf_h298 - gas_h298)/9.6e4

species_dict = {}
for spec_name in be_df['Species']:
    smiles = spec_smiles_dict[spec_name.strip()]
    spec = get_thermo(smiles)
    spec.label = spec_name
    species_dict[spec_name.strip()] = spec
    
# # manually add surface site to species_dict
# species_dict['*'] = get_thermo(spec_smiles_dict['*'])

gas_species_dict = {}
for spec_name in gas_smiles_dict.keys():
    smiles = gas_smiles_dict[spec_name.strip()]
    spec = get_thermo(smiles)
    spec.label = spec_name
    gas_species_dict[spec_name.strip()] = spec

# make binding energy dictionary from park data
be_dict_park = {}
for i in range(len(be_df)):
    species = be_df['Species'][i].strip()
    be_park = be_df["BE"][i]
    be_dict_park[species] = be_park

# update thermo to be closer to bark BE values 
new_thermo_spec_dict = {}
for name, spec in species_dict.items():
    spec_new = update_thermo(
        spec, 
        name,
        be_dict[name], 
        be_dict_park[name],   
    )
    new_thermo_spec_dict[name] = spec_new

# combine gas and surface species dicts
combined_species_dict = {**new_thermo_spec_dict, **gas_species_dict}

# now that we've solidified the thermo, convert to nasa so chemkin conversion
# is a little easier
for spec in combined_species_dict.values():
    convert_to_nasa(spec)


# pull the information for rea ctants, products, 
# and arrhenius prefactors for the equations below
rxn_spec_dict = {}
rxn_dict = {}
rxn_dict_coeff = {}
rxn_list = {}
for index, row in rxn_df.iterrows():
    rxn_raw = row['eqtn']
    rxn = rxn_raw.strip()
    reactants, products = rxn.split("<=>")
    reac_spl = reactants.split("+")
    prod_spl = products.split("+")
    
    # retain to list with stoichiometric coeff
    # just in case we need it
    reac_spl_coeff = reac_spl
    prod_spl_coeff = prod_spl
    
    # expand split reactant/product string so 
    # reactants with "2" as prefix become two 
    # separate strings 
    # e.g. 2OH --> OH, OH
    for reac in reac_spl:
        if reac.startswith("2"):
            reac_dup = reac.replace("2","")
            reac_spl.remove(reac)
            reac_spl.extend([reac_dup]*2)
            
    for prod in prod_spl:
        if prod.startswith("2"):
            prod_dup = prod.replace("2","")
            prod_spl.remove(prod)
            prod_spl.extend([prod_dup]*2) 
    
    rxn_dict[rxn] = [reac_spl, prod_spl] 
    rxn_dict_coeff[rxn] = [reac_spl_coeff, prod_spl_coeff] 
    
    if row['Af'] == 'N/A' and row['stick']:
        # if no rate info and sticking coefficient
        A = 1.0 # units of mol/m^2/s
    elif row['Af'] != 'N/A' and row['stick']:
        # if we supply a sticking coefficient         
        A = float(row['Af'])
    else:
        # we are making a concession here. rates that do 
        # not have an A-factor or Ea specified are quasi-
        # equilibrated, so I am setting the A-factor to the 
        # highest value (1e22 1/s) in the mechanism, and
        # making it barrierless (Ea=0 eV)

        if len(reac_spl) > 1:
            A = (float(row['Af'] / site_density_si), 'm^2/(mol*s)') # units of mol/m^2/s
        else:
            A = (float(row['Af'] / site_density_si), 's^-1') # units of mol/m^2/s
        
    Ea = (float(row['Ef (eV)'] * 9.6e4), 'J/mol') # units of J/mol
        
    
    rxn_spec_dict[rxn] = [
        [combined_species_dict[reac] for reac in reac_spl], 
        [combined_species_dict[prod] for prod in prod_spl],
    ]
    
    rxn_obj = make_reaction(
        rxn_spec_dict[rxn][0], 
        rxn_spec_dict[rxn][1], 
        rxn, 
        A, 
        Ea, 
        stick = row['stick'],
    )
    rxn_list[rxn] = rxn_obj

# finally, make inputs into lists for chemkin file write
chemkin_specs = []
for spec in combined_species_dict.values():
    chemkin_specs.append(spec)

chemkin_rxns = []
for rxn in rxn_list.values():
    chemkin_rxns.append(rxn)

# write chemkin file
# make inputs into lists for chemkin file write
chemkin_specs = []
for spec in gas_species_dict.values():
    chemkin_specs.append(spec)

chemkin_rxns = []

Chemkin.save_chemkin_file(
    chemkin_gas_file, 
    chemkin_specs, 
    chemkin_rxns, 
    verbose=True, 
    check_for_duplicates=True,
    )

# make inputs into lists for chemkin file write
chemkin_specs = []
for spec in new_thermo_spec_dict.values():
    chemkin_specs.append(spec)

chemkin_rxns = []
for rxn in rxn_list.values():
    chemkin_rxns.append(rxn)
    
Chemkin.save_chemkin_surface_file(
    chemkin_surface_file, 
    chemkin_specs, 
    chemkin_rxns, 
    verbose=True, 
    check_for_duplicates=True,
    surface_site_density=site_density_object,
    )

parser = ck2cti.Parser()
parser.convertMech(
    chemkin_gas_file, 
    outName=cantera_file, 
    quiet=True, 
    permissive=True, 
    surfaceFile=chemkin_surface_file
)


# test that model works by attempting to load it
gas = ct.Solution(cantera_file, "gas")
surf = ct.Interface(cantera_file,"surface1", [gas])