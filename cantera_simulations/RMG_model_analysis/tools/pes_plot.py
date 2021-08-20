import cantera as ct
import os
import sys
sys.path.append(f'{os.getcwd()}/tools/PyEnergyDiagram')
print(os.getcwd())
from energydiagram import ED
import logging

############################################
#
#   Plots a potential energy surface 
#   (enthalpy vs rxn coordinate) for a 
#   given cti file mechanism
#   
#   uses https://github.com/giacomomarchioro/PyEnergyDiagrams
#   To install: 
#   1. activate desired anaconda environment
#   2. pip install git+https://github.com/giacomomarchioro/PyEnergyDiagrams
#
############################################

class pes_reaction_combine():
    """
    feed in a cantera reaction, get an object out of it that we can use for making a chart

    arguments: 
    reaction - a ct reaction object
    phase_gas - gas phase in mechanism (for looking up species)
    phase_surf - solid phase in mechanism (for looking up species)

    properties:
    reactants - dict, reactant string as key (e.g. "CO2+H2O"), combined Hf as value
    products - dict, reactant string as key (e.g. "CO2+H2O"), combined Hf as value
    barrier - float, Ea for reaction as value
    equation - string, cantera reaction equation
    links - list of ids for connecting the diagram
            [reactant id, Ea id, product id]
    """
    def __init__(
        self,
        reaction,
        phase_gas,
        phase_surf,
        ):

        self.equation = reaction.equation
        self.reactants  = {}
        self.products  = {}
        self.links = [-1,-1,-1]

        # lookup each reactant, put in dictionary as 
        # {species_name : enthalpy at phase temp (in Ev) * stoich coeff}
        total_reac_enth = 0
        reac_str = ""
        for i in reaction.reactants: 
            if i in phase_gas.species_names:
                phase = phase_gas
            elif i in phase_surf.species_names:
                phase = phase_surf
            else: 
                logging.error(f"species {i} cannot be found")

            reac_str += f"{i} "
            total_reac_enth += reaction.reactants[i] * (phase.species(i).thermo.h(phase.T)/1000**2)/96
        
        self.links[0] = -1

        reac_str = reac_str.strip()
        self.reactants[reac_str] = total_reac_enth 
        self.reactants[reac_str] = round(self.reactants[reac_str],3)
        self.links[1] = -1

        total_prod_enth = 0
        prod_str = ""
        for i in reaction.products: 
            if i in phase_gas.species_names:
                phase = phase_gas
            elif i in phase_surf.species_names:
                phase = phase_surf
            else: 
                logging.error(f"species {i} cannot be found")

            prod_str += f"{i} "
            total_prod_enth += reaction.products[i] * (phase.species(i).thermo.h(phase.T)/1000**2)/96
        
        self.links[2] = -1
        prod_str = prod_str.strip()
        self.products[prod_str] = total_prod_enth 
        self.products[prod_str] = round(self.products[prod_str],3)

        # reaction activation energy. need to add to 
        # reactant enthalpy to get barrier
        self.barrier = (reaction.rate.activation_energy/1000**2)/96 + total_reac_enth
        self.barrier = round(self.barrier, 3)

class pes_reaction():
    """
    feed in a cantera reaction, get an object out of it that we can use for making a chart

    arguments: 
    reaction - a ct reaction object
    phase_gas - gas phase in mechanism (for looking up species)
    phase_surf - solid phase in mechanism (for looking up species)

    properties:
    reactants - dict, species names as keys, Hf as value
    products - dict, species names as keys, Hf as value
    barrier - float, Ea for reaction as value
    equation - string, cantera reaction equation
    links - list of ids for connecting the diagram
            [[reactant ids], Ea id, [product ids]]
    """
    def __init__(
        self,
        reaction,
        phase_gas,
        phase_surf,
        ):

        self.equation = reaction.equation
        self.reactants  = {}
        self.products  = {}

        # lookup each reactant, put in dictionary as 
        # {species_name : enthalpy at phase temp (in Ev) * stoich coeff}
        total_reac_enth = 0
        for i in reaction.reactants: 
            if i in phase_gas.species_names:
                phase = phase_gas
            elif i in phase_surf.species_names:
                phase = phase_surf
            else: 
                logging.error(f"species {i} cannot be found")
            
            self.reactants[i] = reaction.reactants[i] * (phase.species(i).thermo.h(phase.T)/1000**2)/96
            total_reac_enth += self.reactants[i]
            self.reactants[i] = round(self.reactants[i],3)

        for i in reaction.products: 
            if i in phase_gas.species_names:
                phase = phase_gas
            elif i in phase_surf.species_names:
                phase = phase_surf
            else: 
                logging.error(f"species {i} cannot be found")

            self.products[i] = reaction.products[i] * (phase.species(i).thermo.h(phase.T)/1000**2)/96
            self.products[i] = round(self.products[i],3)

        # reaction activation energy. need to add to 
        # reactant enthalpy to get barrier
        self.barrier = (reaction.rate.activation_energy/1000**2)/96 + total_reac_enth
        self.barrier = round(self.barrier, 3)

class pes_plot():
    """
    Plots a potential energy surface 
    (enthalpy vs rxn coordinate) for a 
    given cti file mechanism
    """
    def __init__(
        self,
        yaml_file,
        temp=528,
        press=75,
        ):
        """
        initialize model
        yaml_file = cti or yaml file for mechanism
        temp = temperature (K)
        press = pressure (atm)
        """

        # set initial temps & pressures
        self.temp = temp # kelvin
        self.pressure =  press * ct.one_atm  # Pascals

        # create thermo phases
        self.yaml_file = yaml_file
        self.gas = ct.Solution(yaml_file, "gas")
        self.surf = ct.Interface(yaml_file,"surface1", [self.gas])

        # initialize T and P for each phase
        self.gas.TP = self.temp, self.pressure
        self.surf.TP = self.temp, self.pressure

    def get_h_ev(self, species, temp):
        """
        gets species enthalpy in eV. 
        species is a cantera Species object
        """
        h_eV = (species.thermo.h(temp)/1000**2)/96
        print(f'{species.name} enthalpy = {h_eV} eV')
        return h_eV

    def get_ea_ev(self, reaction):
        """
        gets reaction Ea in eV. 
        reaction is a cantera Reaction object
        """
        Ea_eV = (reaction.rate.activation_energy/1000**2)/96
        print(f'{reaction.equation} enthalpy = {Ea_eV} eV')
        return Ea_eV

    def find_reactions(self, species, temp):
        """
        find all reactions that involve a certain species.
        species is a species object
        rxns is a dictionary, reaction equation is the key, reaction object is the value
        """
        rxns = {}
        species_names = [i.name for i in species]
        species_ea = self.get_h_ev(species[1], temp)
        
        for i,j in enumerate(self.gas.reaction_equations()):
            if all(x in j for x in species_names):
                rxns[j] = self.gas.reaction(i)
                Ea = species_ea + (self.gas.reaction(i).rate.activation_energy/1000**2)/96
                print(j, Ea)
                
        for i,j in enumerate(self.surf.reaction_equations()):
            if all(x in j for x in species_names):
                rxns[j] = self.surf.reaction(i)
                Ea = species_ea + (self.surf.reaction(i).rate.activation_energy/1000**2)/96
                print(j, Ea)
                
        return rxns

    def plot_pes_diagram(
        self, 
        species, 
        width, 
        height, 
        offset=None,
        dimension=None,
        space=None,
        combined=True,
        ):
        """
        plots a potential energy surface given an input for species, 

        inputs:
        species - str or a list of strs matching species name in cantera mechanism.
        width - float matplotlib plot width in inches
        height - float matplotlib plot height in inches
        offset - float vertical distance that energy level and upper/lower labels are spaced
        dimension - float width of platform used for energy level 
        space - float distance between bars for energy levels
        combined - bool if true combine all reactants to a single energy level. do the same for products.
        """
        # create reaction diagram object
        self.diagram = ED()

        # get a list of all reactions containing the two species identified
        species_obj = []
        for i in species:
            if i in self.gas.species_names:
                species_obj.append(self.gas.species(i))
            elif i in self.surf.species_names:
                species_obj.append(self.surf.species(i))

        rxns = self.find_reactions(species_obj, self.temp)

        # for each reaction, make a "pes_rxn" object. add to a dict.  
        # contains information for enthalpy and barriers
        pes_rxn_dict = {}

        for i,j in rxns.items():
            # create pes_rxns
            if combined:
                pes_rxn_dict[i] = pes_reaction_combine(j, self.gas, self.surf)
            else: 
                pes_rxn_dict[i] = pes_reaction(j, self.gas, self.surf)

        first = True
        link_num = 0
        for i,j in rxns.items():
            # generate a pes plot for each pes_reaction
            for k,l in pes_rxn_dict[i].reactants.items():
                reac = k
                H_r = l

                # if it is the first one, make a new energy level
                if first:
                    self.diagram.add_level( H_r, k)
                    first = False
                else:
                    self.diagram.add_level(H_r, k,'l')
                    
            pes_rxn_dict[i].links[0] = link_num
            link_num+=1


        first = True
        for i,j in rxns.items():   
            # plot rxn Ea. for it to show up between species, should be here
            rxn_eq = pes_rxn_dict[i].equation
            rxn_Ea = pes_rxn_dict[i].barrier

            # if it is the first one, make a new energy level
            if first: 
                self.diagram.add_level(rxn_Ea, rxn_eq)
                first = False
            else:
                self.diagram.add_level(rxn_Ea, rxn_eq, 'l')
            
            pes_rxn_dict[i].links[1] = link_num
            link_num+=1

        first = True
        for i,j in rxns.items():  
            # generate a pes plot for each pes_reaction 
            for k,l in pes_rxn_dict[i].products.items():
                prod = k
                H_p = l

                # if it is the first one, make a new energy level
                if first:
                    self.diagram.add_level(H_p, prod)
                    first = False
                else:
                    self.diagram.add_level(H_p, prod,'l')

            # add link id for line drawing
            pes_rxn_dict[i].links[2] = link_num
            link_num+=1

        
        for i in pes_rxn_dict.values():
            # get connections between each reac - Ea and each Ea-product
            self.diagram.add_link(i.links[0],i.links[1])
            self.diagram.add_link(i.links[1],i.links[2])

        # optional arguements 
        if space: 
            self.diagram.space = space
        if offset:
            self.diagram.offset = offset
        if dimension:
            self.diagram.dimension = dimension

        self.diagram.plot(show_IDs=True, ylabel="Energy / $eV$", width=width, height=height)




    