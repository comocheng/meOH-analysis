
import pandas as pd
import numpy as np
import time
import cantera as ct
from matplotlib import pyplot as plt
import csv
import math
import os
import sys
import re
import itertools
import logging
from collections import defaultdict
import git
import time

from rmgpy.molecule import Molecule
from rmgpy.data.base import Database

class sbr:

###############################################
# Stirred batch reactor script
# Chris Blais
# Northeastern University
###############################################
    def __init__(
        self,
        yaml_file,
        rmg_model_path, 
        t_array=[528],
        p_array=[75],
        v_array=[0.00424],
        h2_array=[0.75],
        co2_array=[0.5],
        rtol=1.0e-11,
        atol=1.0e-22,
        reactor_type=0,
        energy="off",
        sensitivity=0,
        sensatol=1e-6,
        sensrtol=1e-6,
        reactime=1e5,
        grabow=False,
        sens_species="CH3OH(8)",
        graaf = False,
        ):
        """
        initialize object 
        yaml_file = cti or yaml file for mechanism
        t_array = temperature (K)
        p_array = pressure array (atm)
        v_array = volume flow rate array (m^3/s)
        h2_array=[0.75],
        co2_array=[0.5],
        rtol = relative tolerance
        atol = absolute tolerance
        reactor_type = 0 - Reactor, 1 - IdealGasReactor, 
                       2 - ConstPressureReactor, 3 - IdealGasConstPressureReactor
        energy = str, "on" if non-isothermal, "off" if isothermal
        sensitivity = int (0-3) perform sensitivity analysis on: 
                        0: nothing
                        1: Kinetics
                        2: Thermo
                        3: Kinetics + Thermo
        sensatol = sensitivity atol
        sensrtol = sensitivity rtol
        reactime = time to run the reactor
        grabow = use grabow data (not tied to a commit)
        graaf = bool, use graaf inlet conditions 
        """
        try:
            array_i = int(os.getenv("SLURM_ARRAY_TASK_ID"))
        except TypeError:
            array_i = 0

        self.t_array = t_array
        self.p_array = p_array
        self.v_array = v_array 
        self.h2_array = h2_array
        self.co2_array = co2_array
        # generate settings array. if using Graaf values, load them from csv
        if graaf:
            self.settings = self.load_graaf_data()
        else: 
            self.settings = list(
                itertools.product(
                    self.t_array, 
                    self.p_array, 
                    self.v_array, 
                    self.h2_array, 
                    self.co2_array
                    )
                )
        
        # get information for git repository that model comes from. 
        # alternatively, if the grabow model is used, use a dummy 
        # date, commit and hash
        if grabow: 
            # format grabow model the same as the others
            self.rmg_model_path = rmg_model_path
            repo = git.Repo(rmg_model_path)
            date = time.localtime()
            git_date = "0000_00_00_0000"
            git_sha = '000000'
            git_msg = "Grabow model"
            self.git_file_string = f"{git_date}_{git_sha}_{git_msg}"      
        else:
            # get git commit hash and message
            self.rmg_model_path = rmg_model_path
            repo = git.Repo(rmg_model_path)
            date = time.localtime(repo.head.commit.committed_date)
            git_date = time.strftime("%Y_%m_%d_%H%M", date)
            git_sha = str(repo.head.commit)[0:6]
            git_msg = str(repo.head.commit.message)[0:50].replace(" ", "_").replace("'", "_").replace("\n", "")
            self.git_file_string = f"{git_date}_{git_sha}_{git_msg}"

        # set sensitivity string for file path name
        self.sensitivity = sensitivity
        if self.sensitivity == 1: 
            self.sensitivity_str = "kinetics_sensitivity"
        elif self.sensitivity == 2: 
            self.sensitivity_str = "thermo_sensitivity"
        elif self.sensitivity == 3: 
            self.sensitivity_str = "all_sensitivity"
        else: 
            self.sensitivity_str = "no_sensitivity"

        self.sens_species = sens_species
        # constants
        pi = math.pi

        # set initial temps& pressures
        self.temp = self.settings[array_i][0] # kelvin
        self.temp_str = str(self.temp)[0:3]
        self.pressure =  self.settings[array_i][1] * ct.one_atm  # Pascals

        # set mole fractions of each species
        self.X_h2 = self.settings[array_i][3]
        self.x_h2_str = str(self.X_h2)[0:3].replace(".", "_")
        self.x_CO_CO2_str = str(self.settings[array_i][4])[0:3].replace(".", "_")

        # Per Grabow experiments, add in H2O for X=0.75 H2 run
        if self.X_h2 == 0.75:
            self.X_h2o = 0.05
        else:
            self.X_h2o = 0

        self.X_co = (1 - (self.X_h2 + self.X_h2o)) * (self.settings[array_i][4])
        self.X_co2 = (1 - (self.X_h2 + self.X_h2o)) * (1 - self.settings[array_i][4])

        # molecular weights for mass flow calculations
        mw_co = 28.01e-3  # [kg/mol]
        mw_co2 = 44.01e-3  # [kg/mol]
        mw_h2 = 2.016e-3  # [kg/mol]
        mw_h2o = 18.01528e-3  # [kg/mol]

        # define ratios of reactants for plots
        # CO2/(CO+CO2) and (CO2+CO)/H2
        self.co2_ratio = self.X_co2 / (self.X_co + self.X_co2)
        self.h2_ratio = (self.X_co2 + self.X_co) / self.X_h2

        # CO/CO2/H2/H2: typical is
        if grabow:
            self.concentrations_rmg = {
                "CO": self.X_co, 
                "CO2": self.X_co2, 
                "H2": self.X_h2, 
                "H2O": self.X_h2o,
                }
        else:
            self.concentrations_rmg = {
                "CO(3)": self.X_co, 
                "CO2(4)": self.X_co2, 
                "H2(2)": self.X_h2, 
                "H2O(5)": self.X_h2o,
                }

        # create thermo phases
        self.yaml_file = yaml_file
        self.gas = ct.Solution(yaml_file, "gas")
        self.surf = ct.Interface(yaml_file,"surface1", [self.gas])

        # initialize T and P
        self.gas.TPX = self.temp, self.pressure, self.concentrations_rmg
        self.surf.TP = self.temp, self.pressure

        # if a mistake is made with the input moles, 
        # cantera will normalize the mole fractions. 
        # make sure that we are reporting
        # the normalized values
        if grabow:
            self.X_co = float(self.gas["CO"].X)
            self.X_co2 = float(self.gas["CO2"].X)
            self.X_h2 = float(self.gas["H2"].X)
            self.X_h2o = float(self.gas["H2O"].X)
        else:           
            self.X_co = float(self.gas["CO(3)"].X)
            self.X_co2 = float(self.gas["CO2(4)"].X)
            self.X_h2 = float(self.gas["H2(2)"].X)
            self.X_h2o = float(self.gas["H2O(5)"].X)

        # create gas inlet
        self.inlet = ct.Reservoir(self.gas)

        # create gas outlet
        self.exhaust = ct.Reservoir(self.gas)

        # Reactor volume (divide by 2 per Graaf paper)
        self.rradius = 35e-3
        self.rlength = 70e-3
        self.rvol = (self.rradius ** 2) * math.pi * self.rlength/2

        # Catalyst Surface Area
        self.site_density = (
            self.surf.site_density * 1000
        )  # [mol/m^2]cantera uses kmol/m^2, convert to mol/m^2
        self.cat_weight = 4.24e-3  # [kg]
        self.cat_site_per_wt = (300 * 1e-6) * 1000  # [mol/kg] 1e-6mol/micromole, 1000g/kg
        self.cat_area = self.site_density / (self.cat_weight * self.cat_site_per_wt)  # [m^3]
        self.cat_area_str = "%s" % "%.3g" % self.cat_area

        # reactor initialization
        self.reactor_type = reactor_type
        self.energy = energy
        if self.reactor_type == 0:
            self.r = ct.Reactor(self.gas, energy=self.energy)
            self.reactor_type_str = "Reactor"
        elif reactor_type == 1:
            self.r = ct.IdealGasReactor(self.gas, energy=self.energy)
            self.reactor_type_str = "IdealGasReactor"
        elif reactor_type == 2:
            self.r = ct.ConstPressureReactor(self.gas, energy=self.energy)
            self.reactor_type_str = "ConstPressureReactor"
        elif reactor_type == 3:
            self.r = ct.IdealGasConstPressureReactor(self.gas, energy=self.energy)
            self.reactor_type_str = "IdealGasConstPressureReactor"

        # calculate the available catalyst area in a differential reactor
        self.rsurf = ct.ReactorSurface(self.surf, self.r, A=self.cat_area)
        self.r.volume = self.rvol
        if grabow:
            self.surf.coverages = "X:1.0"
        else:
            self.surf.coverages = "X(1):1.0"

        # flow controllers (Graaf measured flow at 293.15 and 1 atm)
        one_atm = ct.one_atm
        FC_temp = 293.15
        self.volume_flow = self.settings[array_i][2]  # [m^3/s]
        self.molar_flow = self.volume_flow * one_atm / (8.3145 * FC_temp)  # [mol/s]
        self.mass_flow = self.molar_flow * (
            self.X_co * mw_co + self.X_co2 * mw_co2 + self.X_h2 * mw_h2 + self.X_h2o * mw_h2o
        )  # [kg/s]
        self.mfc = ct.MassFlowController(self.inlet, self.r, mdot=self.mass_flow)

        # A PressureController has a baseline mass flow rate matching the 'master'
        # MassFlowController, with an additional pressure-dependent term. By explicitly
        # including the upstream mass flow rate, the pressure is kept constant without
        # needing to use a large value for 'K', which can introduce undesired stiffness.
        self.outlet_mfc = ct.PressureController(self.r, self.exhaust, master=self.mfc, K=0.01)


        # initialize reactor network
        self.sim = ct.ReactorNet([self.r])

        # set relative and absolute tolerances on the simulation
        self.sim.rtol = rtol
        self.sim.atol = atol
        
        self.sensrtol = sensrtol
        self.sensatol = sensatol

        # set reactime for transient reactor
        self.reactime = reactime
   
    def load_graaf_data(self):
        """
        Julia Treese
        Edited part to get Graaf conditions
        get Graaf conditions into a list of lists to run
        """
        file_name_feed1 = "../Graaf_data/Feed_1.xlsx"
        file_name_feed2 = "../Graaf_data/Feed_2.xlsx"
        file_name_feed3 = "../Graaf_data/Feed_3.xlsx"
        file_name_feed4 = "../Graaf_data/Feed_4.xlsx"
        file_name_feed5 = "../Graaf_data/Feed_5.xlsx"

        df_1 = pd.read_excel(file_name_feed1, engine='openpyxl')
        df_2 = pd.read_excel(file_name_feed2, engine='openpyxl')
        df_3 = pd.read_excel(file_name_feed3, engine='openpyxl')
        df_4 = pd.read_excel(file_name_feed4, engine='openpyxl')
        df_5 = pd.read_excel(file_name_feed5, engine='openpyxl')


        # Needed: [T, P, V, YH2, YCO2] -- Create a list of lists
        # Should be columns 2 (T), 1 (P), 3 (V), 6 (YH2), 5 (YCO2)
        # Each list is the conditions of one experimental Graaf run

        # List of dataframes with feed conditions
        df_list = [df_1, df_2, df_3, df_4, df_5]

        # Loop through dataframes and create a list of conditions based on Graaf runs
        # Loop through each row in the dataframes and add that row's conditions to the list of lists

        settings = []

        for i in range(len(df_list)):
            df = df_list[i]
            for row in range(len(df)):
                row_conditions = [df.iloc[row, df.columns.get_loc('T(K)')], 
                                df.iloc[row, df.columns.get_loc('p (bar)')], 
                                df.iloc[row,df.columns.get_loc('10^6 * V (M^3/s)')], 
                                df.iloc[row,df.columns.get_loc('feed Yco2')], 
                                df.iloc[row,df.columns.get_loc('CO2/(CO+CO2)')],]
                settings.append(row_conditions)
        return settings

    def save_pictures(self, git_path="", species_path="", overwrite=False):
        """
        Save a folder full of molecule pictures, needed for the pretty dot files.

        Saves them in the results directory, in a subfolder "species_pictures".
        Unless you set overwrite=True, it'll leave alone files that are
        already there.
        """
        dictionary_filename = git_path + "/base/chemkin/species_dictionary.txt"
        specs = Database().get_species(dictionary_filename, resonance=False)

        images_dir = os.path.join(species_path)
        os.makedirs(images_dir, exist_ok=True)
        for name, species in specs.items():
            filepath = os.path.join(images_dir, name + ".png")
            if not overwrite and os.path.exists(filepath):
                continue
            species.molecule[0].draw(filepath)

    def prettydot(self, species_path="", dotfilepath="", strip_line_labels=False):
        """
        Make a prettier version of the dot file (flux diagram)

        Assumes the species pictures are stored in a directory
        called 'species_pictures' alongside the dot file.
        """
        pictures_directory = f'{species_path}/'

        if strip_line_labels:
            print("stripping edge (line) labels")

        reSize = re.compile('size="5,6"\;page="5,6"')
        reNode = re.compile(
            '(?P<node>s\d+)\ \[\ fontname="Helvetica",\ label="(?P<label>[^"]*)"\]\;'
        )

        rePicture = re.compile("(?P<smiles>.+?)\((?P<id>\d+)\)\.png")
        reLabel = re.compile("(?P<name>.+?)\((?P<id>\d+)\)$")

        species_pictures = dict()
        for picturefile in os.listdir(pictures_directory):
            match = rePicture.match(picturefile)
            if match:
                species_pictures[match.group("id")] = picturefile
            else:
                pass
                # print(picturefile, "didn't look like a picture")

        filepath = dotfilepath

        if not open(filepath).readline().startswith("digraph"):
            raise ValueError("{0} - not a digraph".format(filepath))

        infile = open(filepath)
        prettypath = filepath.replace(".dot", "", 1) + "-pretty.dot"
        outfile = open(prettypath, "w")

        for line in infile:
            (line, changed_size) = reSize.subn('size="12,12";page="12,12"', line)
            match = reNode.search(line)
            if match:
                label = match.group("label")
                idmatch = reLabel.match(label)
                if idmatch:
                    idnumber = idmatch.group("id")
                    if idnumber in species_pictures:
                        line = (
                            f'%s [ image="{pictures_directory}%s" label="" width="0.5" height="0.5" imagescale=true fixedsize=false color="none" ];\n'
                            % (match.group("node"), species_pictures[idnumber])
                        )

            # rankdir="LR" to make graph go left>right instead of top>bottom
            if strip_line_labels:
                line = re.sub('label\s*=\s*"\s*[\d.]+"', 'label=""', line)

            # change colours
            line = re.sub('color="0.7,\ (.*?),\ 0.9"', r'color="1.0, \1, 0.7*\1"', line)

            outfile.write(line)

        outfile.close()
        infile.close()
        print(f"Graph saved to: {prettypath}")
        os.system(f'dot {prettypath} -Tpng -o{prettypath.replace(".dot", "", 1) + ".png"} -Gdpi=300')
        return prettypath

    def show_flux_diagrams(self, suffix="", embed=False):
        """
        Shows the flux diagrams in the notebook.
        Loads them from disk.
        Does not embed them, to keep the .ipynb file small,
        unless embed=True. Use embed=True if you might over-write the files,
        eg. you want to show flux at different points.
        """
        import IPython

        for element in "CHONX":
            for phase_object in (self.gas, self.surf):
                phase = phase_object.name
                img_file = (
                    f"reaction_path_{element}_{phase}{'_' if suffix else ''}{suffix}.png"
                )
                display(IPython.display.HTML(f"<hr><h2>{element} {phase}</h2>"))
                if embed:
                    display(IPython.display.Image(filename=img_file, width=400, embed=True))
                else:
                    display(IPython.display.Image(url=img_file, width=400, embed=False))

            # Now do the combined
            img_file = f"reaction_path_mass{'_' if suffix else ''}{suffix}.png"
            display(IPython.display.HTML(f"<hr><h2>Combined mass</h2>"))
            if embed:
                display(IPython.display.Image(filename=img_file, width=400, embed=True))
            else:
                display(IPython.display.Image(url=img_file, width=400, embed=False))

    def save_flux_diagrams(self, *phases, suffix="", timepoint="", species_path=""):
        """
        Saves the flux diagrams. The filenames have a suffix if provided,
        so you can keep them separate and not over-write.
        """
        for element in "CHONX":
            for phase_object in phases:
                phase = phase_object.name

                diagram = ct.ReactionPathDiagram(phase_object, element)
                diagram.threshold = 1e-5
                diagram.title = f"Reaction path diagram following {element} in {phase}"
                # diagram.label_threshold = 0.0001

                dot_file = f"{suffix}/reaction_path_{element}_{phase}_{timepoint}.dot"
                img_file = f"{suffix}/reaction_path_{element}_{phase}_{timepoint}.png"
                dot_bin_path = (
                    "/Users/blais.ch/anaconda3/pkgs/graphviz-2.40.1-hefbbd9a_2/bin/dot"
                )
                img_path = os.path.join(os.getcwd(), img_file)
                diagram.write_dot(dot_file)

                #also make a prettydot file
                self.prettydot(species_path, dot_file, strip_line_labels=False)

                # print(diagram.get_data())

                print(
                    f"Wrote graphviz input file to '{os.path.join(os.getcwd(), dot_file)}'."
                )
                os.system(f"dot {dot_file} -Tpng -o{img_file} -Gdpi=200")
                print(f"Wrote graphviz output file to '{img_path}'.")

    def run_reactor_transient(self):
        """
        Run single reactor to a set time
        """

        # self.git_file_string = self.git_file_string.replace("&", "\&")
        # set paths for saving species pictures, fluxes, and csv files
        self.species_path = (
            os.path.dirname(os.path.abspath(__file__))
            + f"/{self.git_file_string}/species_pictures"
        )
        self.results_path = (
            os.path.dirname(os.path.abspath(__file__))
        + f"/{self.git_file_string}/transient/{self.reactor_type_str}/energy_{self.energy}/{self.sensitivity_str}/{self.temp_str}/results"
        )
        results_path_csp = (
            os.path.dirname(os.path.abspath(__file__))
            + f"/{self.git_file_string}/transient/{self.reactor_type_str}/energy_{self.energy}/{self.sensitivity_str}/{self.temp_str}/csp"
        )

        self.flux_path = (
            os.path.dirname(os.path.abspath(__file__))
            + f"/{self.git_file_string}/transient/{self.reactor_type_str}/energy_{self.energy}/{self.sensitivity_str}/{self.temp_str}/flux_diagrams/{self.x_h2_str}/{self.x_CO_CO2_str}"
        )
        
        # create folders if they don't already exist
        try:
            os.makedirs(self.species_path, exist_ok=True)
            self.save_pictures(git_path=self.rmg_model_path, species_path=self.species_path)
        except OSError as error:
            print(error)

        try:
            os.makedirs(self.results_path, exist_ok=True)
        except OSError as error:
            print(error)

        try:
            os.makedirs(results_path_csp, exist_ok=True)
            print(results_path_csp)
        except OSError as error:
            print(error)
            print("no results for CSP saved")

        try:
            os.makedirs(self.flux_path, exist_ok=True)
        except OSError as error:
            print(error)
        
        gas_ROP_str = [i + " ROP [kmol/m^3 s]" for i in self.gas.species_names]

        # surface ROP reports gas and surface ROP. these values might be redundant, not sure.

        gas_surf_ROP_str = [i + " surface ROP [kmol/m^2 s]" for i in self.gas.species_names]
        surf_ROP_str = [i + " ROP [kmol/m^2 s]" for i in self.surf.species_names]

        gasrxn_ROP_str = [i + " ROP [kmol/m^3 s]" for i in self.gas.reaction_equations()]
        surfrxn_ROP_str = [i + " ROP [kmol/m^2 s]" for i in self.surf.reaction_equations()]
        self.output_filename = (
            self.results_path
            + f"/Spinning_basket_area_{self.cat_area_str}_energy_{self.energy}"
            + f"_temp_{self.temp}_h2_{self.x_h2_str}_COCO2_{self.x_CO_CO2_str}.csv"
        )
        self.output_filename_csp = (
            results_path_csp
            + f"/CSP_Spinning_basket_area_{self.cat_area_str}_energy_{self.energy}"
            + f"_temp_{self.temp}_h2_{self.x_h2_str}_COCO2_{self.x_CO_CO2_str}.dat"
        )
        outfile = open(self.output_filename, "w")
        outfile_csp = open(self.output_filename_csp, "w")
        writer = csv.writer(outfile)
        writer_csp = csv.writer(outfile_csp, delimiter='\t')

        # log the filepath so we know where to go if there's a problem 
        logging.warning(self.results_path+self.output_filename)
                

        # get list  of preconditions for csv writer
        preconditions = [
            "time (s)",
            "T (K)",
            "P (Pa)",
            "V (M^3/s)",
            "X_co initial",
            "X_co2initial",
            "X_h2 initial",
            "X_h2o initial",
            "CO2/(CO2+CO)",
            "(CO+CO2/H2)",
            "T (C) final",
            "Rtol",
            "Atol",
            "reactor type",
            "energy on?"
            ]

        # Sensitivity atol, rtol, and strings for gas and surface reactions if selected
        # slows down script by a lot
        if self.sensitivity ==1: # kinetic sensitivity
            self.sim.rtol_sensitivity = self.sensrtol
            self.sim.atol_sensitivity = self.sensatol

            # turn on sensitive reactions
            for i in range(self.gas.n_reactions):
                self.r.add_sensitivity_reaction(i)

            for i in range(self.surf.n_reactions):
                self.rsurf.add_sensitivity_reaction(i)

            for j in self.sens_species:
                gasrxn_sens_str = [
                    j + " sensitivity to " + i for i in self.gas.reaction_equations()
                ]
                surfrxn_sens_str = [
                    j + " sensitivity to " + i for i in self.surf.reaction_equations()
                ]
                sens_list = gasrxn_sens_str + surfrxn_sens_str  


            writer.writerow(
                preconditions
                + self.gas.species_names
                + self.surf.species_names
                + gas_ROP_str
                + gas_surf_ROP_str
                + surf_ROP_str
                + gasrxn_ROP_str
                + surfrxn_ROP_str
                + sens_list
            )

        if self.sensitivity ==2: # thermo sensitivity
            self.sim.rtol_sensitivity = self.sensrtol
            self.sim.atol_sensitivity = self.sensatol

            # turn on sensitive species
            for i in range(self.gas.n_species):
                self.r.add_sensitivity_species_enthalpy(i)

            for i in range(self.surf.n_species):
                self.rsurf.add_sensitivity_species_enthalpy(i)

            for j in self.sens_species:
                gastherm_sens_str = [
                    j + " thermo sensitivity to " + i for i in self.gas.species_names
                    ]
                surftherm_sens_str = [
                    j + " thermo sensitivity to " + i for i in self.surf.species_names
                    ]
                sens_list = gastherm_sens_str + surftherm_sens_str


            writer.writerow(
                preconditions
                + self.gas.species_names
                + self.surf.species_names
                + gas_ROP_str
                + gas_surf_ROP_str
                + surf_ROP_str
                + gasrxn_ROP_str
                + surfrxn_ROP_str
                + sens_list
            )

        if self.sensitivity ==3: # both sensitivity
            self.sim.rtol_sensitivity = self.sensrtol
            self.sim.atol_sensitivity = self.sensatol

            # turn on sensitive reactions/species
            for i in range(self.gas.n_reactions):
                self.r.add_sensitivity_reaction(i)

            for i in range(self.surf.n_reactions):
                self.rsurf.add_sensitivity_reaction(i)

            for i in range(self.gas.n_species):
                self.r.add_sensitivity_species_enthalpy(i)

            for i in range(self.surf.n_species):
                self.rsurf.add_sensitivity_species_enthalpy(i)

            for j in self.sens_species:
                gasrxn_sens_str = [
                    j + " sensitivity to " + i for i in self.gas.reaction_equations()
                ]
                surfrxn_sens_str = [
                    j + " sensitivity to " + i for i in self.surf.reaction_equations()
                ]
                gastherm_sens_str = [
                    j + " thermo sensitivity to " + i for i in self.gas.species_names
                    ]
                surftherm_sens_str = [
                    j + " thermo sensitivity to " + i for i in self.surf.species_names
                    ]
                sens_list = gasrxn_sens_str + surfrxn_sens_str  + gastherm_sens_str + surftherm_sens_str


            writer.writerow(
                preconditions
                + self.gas.species_names
                + self.surf.species_names
                + gas_ROP_str
                + gas_surf_ROP_str
                + surf_ROP_str
                + gasrxn_ROP_str
                + surfrxn_ROP_str
                + sens_list
            )

        else: # no sensitivity
            writer.writerow(
                preconditions
                + self.gas.species_names
                + self.surf.species_names
                + gas_ROP_str
                + gas_surf_ROP_str
                + surf_ROP_str
                + gasrxn_ROP_str
                + surfrxn_ROP_str
            )

        writer_csp.writerow([
            "iter", 
            "t", 
            "dt", 
            "Density[kg/m3]", 
            "Pressure[Pascal]", 
            "Temperature[K]",
            ]
            + self.gas.species_names
            + self.surf.species_names
        )

        # set initial time, timestep, step number
        t = 0.0
        dt = 0.1
        iter_ct = 0

        # run the simulation
        first_run = True
        while t < self.reactime:
            # save flux diagrams at beginning of run
            if first_run == True:
                self.save_flux_diagrams(
                    self.gas, 
                    suffix=self.flux_path, 
                    timepoint="beginning", 
                    species_path=self.species_path
                    )
                self.save_flux_diagrams(
                    self.surf, 
                    suffix=self.flux_path, 
                    timepoint="beginning", 
                    species_path=self.species_path)
                first_run = False
            t += dt
            self.sim.advance(t)

            # some models have the special case where they do not have any
            # gas phase reactions. if this is true, skip over gas ROPs
            if len(self.gas.reactions()) > 0: 
                gas_net_rates_of_progress = list(self.gas.net_rates_of_progress)
            else: 
                gas_net_rates_of_progress = []
                
            precondition_values = [
                self.sim.time,
                self.temp,
                self.pressure,
                self.volume_flow,
                self.X_co,
                self.X_co2,
                self.X_h2,
                self.X_h2o,
                self.co2_ratio,
                self.h2_ratio,
                self.gas.T,
                self.sim.rtol,
                self.sim.atol,
                self.reactor_type_str,
                self.energy,
            ]
            # if sensitivity, get sensitivity for sensitive 
            # species i (e.g. methanol) in reaction j

            if self.sensitivity == 1: # kinetic sensitivity
                for i in self.sens_species:
                    g_nrxn = self.gas.n_reactions
                    s_nrxn = self.surf.n_reactions

                    gas_sensitivities = [
                        self.sim.sensitivity(i, j) for j in range(g_nrxn)
                        ]
                    surf_sensitivities = [
                        self.sim.sensitivity(i, j) for j in range(g_nrxn, g_nrxn + s_nrxn)
                    ]

                    sensitivities_all = (
                        gas_sensitivities
                        + surf_sensitivities
                    )

                writer.writerow(
                    precondition_values
                    + list(self.gas.X)
                    + list(self.surf.X)
                    + list(self.gas.net_production_rates)
                    + list(self.surf.net_production_rates)
                    + gas_net_rates_of_progress
                    + list(self.surf.net_rates_of_progress)
                    + sensitivities_all
                )


            elif self.sensitivity == 2: # thermo_sensitivity
                for i in self.sens_species:
                    g_nspec = self.gas.n_species
                    s_nspec = self.surf.n_species

                    gas_therm_sensitivities = [
                        self.sim.sensitivity(i,j)
                        for j in range(g_nrxn+s_nrxn,g_nrxn+s_nrxn+g_nspec)
                        ]
                    surf_therm_sensitivities = [
                        self.sim.sensitivity(i,j)
                        for j in range(g_nrxn+s_nrxn+g_nspec,g_nrxn+s_nrxn+g_nspec+s_nspec)
                        ]

                    sensitivities_all = (
                        gas_therm_sensitivities
                        + surf_therm_sensitivities
                    )

                writer.writerow(
                    precondition_values
                    + list(self.gas.X)
                    + list(self.surf.X)
                    + list(self.gas.net_production_rates)
                    + list(self.surf.net_production_rates)
                    + gas_net_rates_of_progress
                    + list(self.surf.net_rates_of_progress)
                    + sensitivities_all
                )               


            elif self.sensitivity == 3: # all_sensitivity
                for i in self.sens_species:
                    g_nrxn = self.gas.n_reactions
                    s_nrxn = self.surf.n_reactions
                    g_nspec = self.gas.n_species
                    s_nspec = self.surf.n_species

                    gas_sensitivities = [
                        self.sim.sensitivity(i, j) for j in range(g_nrxn)
                        ]
                    surf_sensitivities = [
                        self.sim.sensitivity(i, j) for j in range(g_nrxn, g_nrxn + s_nrxn)
                    ]

                    gas_therm_sensitivities = [
                        self.sim.sensitivity(i,j)
                        for j in range(g_nrxn+s_nrxn,g_nrxn+s_nrxn+g_nspec)
                        ]
                    surf_therm_sensitivities = [
                        self.sim.sensitivity(i,j)
                        for j in range(g_nrxn+s_nrxn+g_nspec,g_nrxn+s_nrxn+g_nspec+s_nspec)
                        ]

                    sensitivities_all = (
                        gas_sensitivities
                        + surf_sensitivities
                        + gas_therm_sensitivities
                        + surf_therm_sensitivities
                    )

                writer.writerow(
                    precondition_values
                    + list(self.gas.X)
                    + list(self.surf.X)
                    + list(self.gas.net_production_rates)
                    + list(self.surf.net_production_rates)
                    + gas_net_rates_of_progress
                    + list(self.surf.net_rates_of_progress)
                    + sensitivities_all
                )

            else: # no sensitivity
                writer.writerow(
                    precondition_values
                    + list(self.gas.X)
                    + list(self.surf.X)
                    + list(self.gas.net_production_rates)
                    + list(self.surf.net_production_rates)
                    + gas_net_rates_of_progress
                    + list(self.surf.net_rates_of_progress)
                )
                
            writer_csp.writerow(
                [
                    iter_ct,
                    self.sim.time,
                    dt,
                    self.gas.density,
                    self.gas.P,
                    self.gas.T,
                ]
                + list(self.gas.X)
                + list(self.surf.X)
            )

            iter_ct += 1

        outfile.close()
        outfile_csp.close()
        
        
        # save flux diagrams at the end of the run
        self.save_flux_diagrams(
            self.gas, 
            suffix=self.flux_path, 
            timepoint="end", 
            species_path=self.species_path
            )
        self.save_flux_diagrams(
            self.surf, 
            suffix=self.flux_path, 
            timepoint="end", 
            species_path=self.species_path
            )
        return

    def run_reactor_ss(self):
        """
        Run single reactor to steady state
        """
        # set paths for saving species pictures, fluxes, and csv files
        self.species_path = (
            os.path.dirname(os.path.abspath(__file__))
            + f"/{self.git_file_string}/species_pictures"
        )
        self.results_path = (
            os.path.dirname(os.path.abspath(__file__))
        + f"/{self.git_file_string}/steady_state/{self.reactor_type_str}/energy_{self.energy}/{self.sensitivity_str}/{self.temp_str}/results"
        )
        self.flux_path = (
            os.path.dirname(os.path.abspath(__file__))
            + f"/{self.git_file_string}/steady_state/{self.reactor_type_str}/energy_{self.energy}/{self.sensitivity_str}/{self.temp_str}/flux_diagrams/{self.x_h2_str}/{self.x_CO_CO2_str}"
        )

        # create folders if they don't already exist
        try:
            os.makedirs(self.species_path, exist_ok=True)
            self.save_pictures(git_path=self.rmg_model_path, species_path=self.species_path)
        except OSError as error:
            print(error)

        try:
            os.makedirs(self.results_path, exist_ok=True)
        except OSError as error:
            print(error)

        try:
            os.makedirs(self.flux_path, exist_ok=True)
        except OSError as error:
            print(error)
        
        gas_ROP_str = [i + " ROP [kmol/m^3 s]" for i in self.gas.species_names]

        # surface ROP reports gas and surface ROP. these values might be redundant, not sure.
        gas_surf_ROP_str = [i + " surface ROP [kmol/m^2 s]" for i in self.gas.species_names]
        surf_ROP_str = [i + " ROP [kmol/m^2 s]" for i in self.surf.species_names]

        gasrxn_ROP_str = [i + " ROP [kmol/m^3 s]" for i in self.gas.reaction_equations()]
        surfrxn_ROP_str = [i + " ROP [kmol/m^2 s]" for i in self.surf.reaction_equations()]
        self.output_filename = (
            self.results_path
            + f"/Spinning_basket_area_{self.cat_area_str}_energy_{self.energy}"
            + f"_temp_{self.temp}_h2_{self.x_h2_str}_COCO2_{self.x_CO_CO2_str}.csv"
        )

        outfile = open(self.output_filename, "w")
        writer = csv.writer(outfile)

        # log the filepath so we know where to go if there's a problem 
        logging.warning(self.results_path+self.output_filename)
                
        # get list  of preconditions for csv writer
        preconditions = [
            "time (s)",
            "T (K)",
            "P (Pa)",
            "V (M^3/s)",
            "X_co initial",
            "X_co2initial",
            "X_h2 initial",
            "X_h2o initial",
            "CO2/(CO2+CO)",
            "(CO+CO2/H2)",
            "T (C) final",
            "Rtol",
            "Atol",
            "reactor type",
            "energy on?"
            ]

        # Sensitivity atol, rtol, and strings for gas and surface reactions if selected
        # slows down script by a lot
        if self.sensitivity ==1: # kinetic sensitivity
            self.sim.rtol_sensitivity = self.sensrtol
            self.sim.atol_sensitivity = self.sensatol

            # turn on sensitive reactions
            for i in range(self.gas.n_reactions):
                self.r.add_sensitivity_reaction(i)

            for i in range(self.surf.n_reactions):
                self.rsurf.add_sensitivity_reaction(i)

            for j in self.sens_species:
                gasrxn_sens_str = [
                    j + " sensitivity to " + i for i in self.gas.reaction_equations()
                ]
                surfrxn_sens_str = [
                    j + " sensitivity to " + i for i in self.surf.reaction_equations()
                ]
                sens_list = gasrxn_sens_str + surfrxn_sens_str  


            writer.writerow(
                preconditions
                + self.gas.species_names
                + self.surf.species_names
                + gas_ROP_str
                + gas_surf_ROP_str
                + surf_ROP_str
                + gasrxn_ROP_str
                + surfrxn_ROP_str
                + sens_list
            )

        elif self.sensitivity ==2: # thermo sensitivity
            self.sim.rtol_sensitivity = self.sensrtol
            self.sim.atol_sensitivity = self.sensatol

            # turn on sensitive species
            for i in range(self.gas.n_species):
                self.r.add_sensitivity_species_enthalpy(i)

            for i in range(self.surf.n_species):
                self.rsurf.add_sensitivity_species_enthalpy(i)

            for j in self.sens_species:
                gastherm_sens_str = [
                    j + " thermo sensitivity to " + i for i in self.gas.species_names
                    ]
                surftherm_sens_str = [
                    j + " thermo sensitivity to " + i for i in self.surf.species_names
                    ]
                sens_list = gastherm_sens_str + surftherm_sens_str


            writer.writerow(
                preconditions
                + self.gas.species_names
                + self.surf.species_names
                + gas_ROP_str
                + gas_surf_ROP_str
                + surf_ROP_str
                + gasrxn_ROP_str
                + surfrxn_ROP_str
                + sens_list
            )

        elif self.sensitivity ==3: # both sensitivity
            self.sim.rtol_sensitivity = self.sensrtol
            self.sim.atol_sensitivity = self.sensatol

            # turn on sensitive reactions/species
            for i in range(self.gas.n_reactions):
                self.r.add_sensitivity_reaction(i)

            for i in range(self.surf.n_reactions):
                self.rsurf.add_sensitivity_reaction(i)

            for i in range(self.gas.n_species):
                self.r.add_sensitivity_species_enthalpy(i)

            for i in range(self.surf.n_species):
                self.rsurf.add_sensitivity_species_enthalpy(i)

            for j in self.sens_species:
                gasrxn_sens_str = [
                    j + " sensitivity to " + i for i in self.gas.reaction_equations()
                ]
                surfrxn_sens_str = [
                    j + " sensitivity to " + i for i in self.surf.reaction_equations()
                ]
                gastherm_sens_str = [
                    j + " thermo sensitivity to " + i for i in self.gas.species_names
                    ]
                surftherm_sens_str = [
                    j + " thermo sensitivity to " + i for i in self.surf.species_names
                    ]
                sens_list = gasrxn_sens_str + surfrxn_sens_str  + gastherm_sens_str + surftherm_sens_str


            writer.writerow(
                preconditions
                + self.gas.species_names
                + self.surf.species_names
                + gas_ROP_str
                + gas_surf_ROP_str
                + surf_ROP_str
                + gasrxn_ROP_str
                + surfrxn_ROP_str
                + sens_list
            )

        else: # no sensitivity
            writer.writerow(
                preconditions
                + self.gas.species_names
                + self.surf.species_names
                + gas_ROP_str
                + gas_surf_ROP_str
                + surf_ROP_str
                + gasrxn_ROP_str
                + surfrxn_ROP_str
            )

        # run the simulation
        self.sim.advance_to_steady_state()

        # some models have the special case where they do not have any
        # gas phase reactions. if this is true, skip over gas ROPs
        if len(self.gas.reactions()) > 0: 
            gas_net_rates_of_progress = list(self.gas.net_rates_of_progress)
        else: 
            gas_net_rates_of_progress = []
            
        precondition_values = [
            self.sim.time,
            self.temp,
            self.pressure,
            self.volume_flow,
            self.X_co,
            self.X_co2,
            self.X_h2,
            self.X_h2o,
            self.co2_ratio,
            self.h2_ratio,
            self.gas.T,
            self.sim.rtol,
            self.sim.atol,
            self.reactor_type_str,
            self.energy,
        ]
        # if sensitivity, get sensitivity for sensitive 
        # species i (e.g. methanol) in reaction j
        if self.sensitivity == 1: # kinetic sensitivity
            for i in self.sens_species:
                g_nrxn = self.gas.n_reactions
                s_nrxn = self.surf.n_reactions

                gas_sensitivities = [
                    self.sim.sensitivity(i, j) for j in range(g_nrxn)
                    ]
                surf_sensitivities = [
                    self.sim.sensitivity(i, j) for j in range(g_nrxn, g_nrxn + s_nrxn)
                ]

                sensitivities_all = (
                    gas_sensitivities
                    + surf_sensitivities
                )

            writer.writerow(
                precondition_values
                + list(self.gas.X)
                + list(self.surf.X)
                + list(self.gas.net_production_rates)
                + list(self.surf.net_production_rates)
                + gas_net_rates_of_progress
                + list(self.surf.net_rates_of_progress)
                + sensitivities_all
            )


        elif self.sensitivity == 2: # thermo_sensitivity
            for i in self.sens_species:
                g_nspec = self.gas.n_species
                s_nspec = self.surf.n_species

                gas_therm_sensitivities = [
                    self.sim.sensitivity(i,j)
                    for j in range(g_nspec)
                    ]
                surf_therm_sensitivities = [
                    self.sim.sensitivity(i,j)
                    for j in range(g_nspec,g_nspec+s_nspec)
                    ]

                sensitivities_all = (
                    gas_therm_sensitivities
                    + surf_therm_sensitivities
                )

            writer.writerow(
                precondition_values
                + list(self.gas.X)
                + list(self.surf.X)
                + list(self.gas.net_production_rates)
                + list(self.surf.net_production_rates)
                + gas_net_rates_of_progress
                + list(self.surf.net_rates_of_progress)
                + sensitivities_all
            )               


        elif self.sensitivity == 3: # all_sensitivity
            for i in self.sens_species:
                g_nrxn = self.gas.n_reactions
                s_nrxn = self.surf.n_reactions
                g_nspec = self.gas.n_species
                s_nspec = self.surf.n_species

                gas_sensitivities = [
                    self.sim.sensitivity(i, j) for j in range(g_nrxn)
                    ]
                surf_sensitivities = [
                    self.sim.sensitivity(i, j) for j in range(g_nrxn, g_nrxn + s_nrxn)
                ]

                gas_therm_sensitivities = [
                    self.sim.sensitivity(i,j)
                    for j in range(g_nrxn+s_nrxn,g_nrxn+s_nrxn+g_nspec)
                    ]
                surf_therm_sensitivities = [
                    self.sim.sensitivity(i,j)
                    for j in range(g_nrxn+s_nrxn+g_nspec,g_nrxn+s_nrxn+g_nspec+s_nspec)
                    ]

                sensitivities_all = (
                    gas_sensitivities
                    + surf_sensitivities
                    + gas_therm_sensitivities
                    + surf_therm_sensitivities
                )

            writer.writerow(
                precondition_values
                + list(self.gas.X)
                + list(self.surf.X)
                + list(self.gas.net_production_rates)
                + list(self.surf.net_production_rates)
                + gas_net_rates_of_progress
                + list(self.surf.net_rates_of_progress)
                + sensitivities_all
            )

        else: # no sensitivity
            writer.writerow(
                precondition_values
                + list(self.gas.X)
                + list(self.surf.X)
                + list(self.gas.net_production_rates)
                + list(self.surf.net_production_rates)
                + gas_net_rates_of_progress
                + list(self.surf.net_rates_of_progress)
            )

        outfile.close()
        
        
        # save flux diagrams at the end of the run
        self.save_flux_diagrams(self.gas, suffix=self.flux_path, timepoint="end", species_path=self.species_path)
        self.save_flux_diagrams(self.surf, suffix=self.flux_path, timepoint="end", species_path=self.species_path)
        return


def run_sbr_test():
    """
    mostly used to run in debugger in vscode
    """
    rmg_model_folder = "/Users/blais.ch/_01_code/05_Project_repos_Github/meOH_repos/MeOH_RMG/meOH-synthesis/"
    cti_file_path = rmg_model_folder + "base/cantera/chem_annotated.cti"

    # initialize reactor
    sbr_ss = sbr(
        cti_file_path,
        rmg_model_folder, 
        t_array=[528],
        p_array=[75],
        v_array=[0.00424],
        h2_array=[0.75],
        co2_array=[0.5],
        rtol=1.0e-11,
        atol=1.0e-22,
        reactor_type=0,
        energy="off",
        sensitivity=False,
        sensatol=1e-6,
        sensrtol=1e-6,
        reactime=1e5,
        grabow=False,
    )
    # run to SS
    sbr_ss.run_reactor_ss()


if __name__ == "__main__":
    # execute only if run as a script
    run_sbr_test()
