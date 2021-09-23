# minimal stirred batch reactor to analyze 5000 rmg runs
###############################################
# Stirred batch reactor script
# Chris Blais, Sevy Harris
# Northeastern University
#
# runs a spinning basket reactor to the specifications in:
# Kinetics of low-pressure methanol synthesis
# Gh Graaf; Ej Stamhuis; Aacm Beenackers
# 1988
# 10.1016/0009-2509(88)85127-3
###############################################


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
import time


class MinSBR:
    """
    A minimal stirred batch reactor
    Takes away some of the options from the regular sbr class
    """
    def __init__(
        self,
        yaml_file,
        rmg_model_path,
        temperatures=[528],
        pressures=[75],
        volumes=[4.24e-6],
        H2_fractions=[0.75],
        CO2_fractions=[0.5],
        catalyst_weights=[4.24e-3],
        rtol=1.0e-11,
        atol=1.0e-22,
        sensitivity=0,
        sensatol=1e-6,
        sensrtol=1e-6,
        sensitive_species="CH3OH(8)",
        reactor_type=1,
        energy="off",
        reactime=1e5,
        timestep=0.1,
    ):
        """
        initialize sbr object
        yaml_file = cti or yaml file for mechanism
        temperatures = list of floats, temperature (K)
        pressures = list of float, pressure array (atm)
        volumes = list of float, volume flow rate array (m^3/s)
        H2_fractions = list of floats,
        CO2_fractions = list of float[0.5],
        rtol = float, relative tolerance
        atol = float, absolute tolerance
        sensitivity = int (0-3) perform sensitivity analysis on:
                        0: nothing
                        1: Kinetics
                        2: Thermo
                        3: Kinetics + Thermo
        sensatol = float, sensitivity atol
        sensrtol = float, sensitivity rtol
        sens_species = list of str, species to get sensitivity for
        energy = str, "on" if non-isothermal, "off" if isothermal
        reactime = flaot, total time to run the reactor (for transient simulation)
        timestep = float, step taken for reactor simulation (for transient simulation)
        catalyst_weights = list of float, weight of the catalyst (in kg)
        """
        try:
            self.SLURM_index = int(os.getenv("SLURM_ARRAY_TASK_ID"))
        except TypeError:
            self.SLURM_index = 0

        self.temperatures = temperatures
        self.pressures = pressures
        self.volumes = volumes
        self.H2_fractions = H2_fractions  # TODO is this inlet H2 or what?
        self.CO2_fractions = CO2_fractions
        self.catalyst_weights = catalyst_weights  # [kg]
        self.rmg_model_path = rmg_model_path

        # generate settings array.
        self.settings = list(
            itertools.product(
                self.temperatures,
                self.pressures,
                self.volumes,
                self.H2_fractions,
                self.CO2_fractions,
                self.catalyst_weights,
            )
        )

        # TODO get information for which RMG run this is

        # set initial temperatures and pressures
        self.temperature = self.settings[self.SLURM_index][0]  # kelvin
        self.temperature_str = str(self.temperature)[0:3]
        self.pressure = self.settings[self.SLURM_index][1] * ct.one_atm  # Pascals

        # set mole fractions of each species
        self.x_H2 = self.settings[self.SLURM_index][3]
        self.x_H2_str = str(self.x_H2)[0:5].replace(".", "_")
        self.x_CO_CO2_str = str(self.settings[self.SLURM_index][4])[0:5].replace(".", "_")

        # TODO consolidate this into something else
        # Per Grabow experiments, add in H2O for X=0.75 H2 run
        if self.x_H2 == 0.75:
            self.x_H2O = 0.05
        else:
            self.x_H2O = 0

        self.x_CO = (1 - (self.x_H2 + self.x_H2O)) * (self.settings[self.SLURM_index][4])
        self.x_CO2 = (1 - (self.x_H2 + self.x_H2O)) * (1 - self.settings[self.SLURM_index][4])

        # molecular weights for mass flow calculations
        MW_CO = 28.01e-3  # [kg/mol]
        MW_CO2 = 44.01e-3  # [kg/mol]
        MW_H2 = 2.016e-3  # [kg/mol]
        MW_H2O = 18.01528e-3  # [kg/mol]

        # define ratios of reactants for plots
        # CO2/(CO+CO2) and (CO2+CO)/H2
        self.CO2_ratio = self.x_CO2 / (self.x_CO + self.x_CO2)
        self.H2_ratio = (self.x_CO2 + self.x_CO) / self.x_H2

        # CO/CO2/H2/H2: typical is
        self.concentrations_rmg = {
            "CO(3)": self.x_CO,
            "CO2(4)": self.x_CO2,
            "H2(2)": self.x_H2,
            "H2O(5)": self.x_H2O,
        }

        # create thermo phases
        self.yaml_file = yaml_file
        self.gas = ct.Solution(yaml_file, "gas")
        self.surf = ct.Interface(yaml_file, "surface1", [self.gas])

        # initialize T and P
        self.gas.TPX = self.temperature, self.pressure, self.concentrations_rmg
        self.surf.TP = self.temperature, self.pressure

        # if a mistake is made with the input moles,
        # cantera will normalize the mole fractions.
        # make sure that we are reporting
        # the normalized values
        self.x_CO = float(self.gas["CO(3)"].X)
        self.x_CO2 = float(self.gas["CO2(4)"].X)
        self.x_H2 = float(self.gas["H2(2)"].X)
        self.x_H2O = float(self.gas["H2O(5)"].X)

        # create gas inlet
        self.inlet = ct.Reservoir(self.gas)

        # create gas outlet
        self.exhaust = ct.Reservoir(self.gas)

        # Reactor volume (divide by 2 per Graaf paper)
        self.rradius = 35e-3
        self.rlength = 70e-3
        self.rvol = (self.rradius ** 2) * math.pi * self.rlength / 2.0

        # Catalyst Surface Area
        self.site_density = (
            self.surf.site_density * 1000
        )  # [mol/m^2]cantera uses kmol/m^2, convert to mol/m^2
        self.cat_weight = self.settings[self.SLURM_index][5]
        self.cat_site_per_wt = (300 * 1e-6) * 1000  # [mol/kg] 1e-6mol/micromole, 1000g/kg
        self.cat_area = self.site_density / (self.cat_weight * self.cat_site_per_wt)  # [m^3]
        self.cat_area_str = "%s" % "%.3g" % self.cat_area

        # reactor initialization
        # always use an IdealGasReactor
        # TODO check about reactor type, see if you can get rid of this block
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

        self.surf.coverages = "X(1):1.0"

        # flow controllers (Graaf measured flow at 293.15 and 1 atm)
        one_atm = ct.one_atm
        FC_temp = 293.15
        self.volume_flow = self.settings[self.SLURM_index][2]  # [m^3/s]
        self.molar_flow = self.volume_flow * one_atm / (8.3145 * FC_temp)  # [mol/s]
        self.mass_flow = self.molar_flow * (
            self.x_CO * MW_CO + self.x_CO2 * MW_CO2 + self.x_H2 * MW_H2 + self.x_H2O * MW_H2O
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
        self.timestep = timestep

    def run_reactor_ss(self):
        """
        Run single reactor to steady state
        """
        # set paths for saving species pictures, fluxes, and csv files

        # save csv results regardless
        # TODO create a sensible path
        self.results_path = os.path.join(self.rmg_model_path, 'steady_state')
        try:
            os.makedirs(self.results_path, exist_ok=True)
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
            + f"/run_number_{self.SLURM_index}_Spinning_basket_area_{self.cat_area_str}_energy_{self.energy}"
            + f"_temp_{self.temperatures}_h2_{self.x_H2_str}_COCO2_{self.x_CO_CO2_str}.csv"
        )

        outfile = open(self.output_filename, "w")
        writer = csv.writer(outfile)

        # log the filepath so we know where to go if there's a problem
        logging.warning(self.results_path + self.output_filename)

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
            "catalyst weight (kg)"
        ]

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
            self.temperature,
            self.pressure,
            self.volume_flow,
            self.x_CO,
            self.x_CO2,
            self.x_H2,
            self.x_H2O,
            self.CO2_ratio,
            self.H2_ratio,
            self.gas.T,
            self.sim.rtol,
            self.sim.atol,
            self.reactor_type_str,
            self.energy,
            self.cat_weight,
        ]

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


def run_sbr_test():
    """
    mostly used to run in debugger in vscode
    """
    rmg_model_folder = "/home/moon/methanol/perturb_5000/run_0000/"
    cti_file_path = "/home/moon/methanol/perturb_5000/run_0000/cantera/chem_annotated.cti"

    # initialize reactor
    sbr_ss = MinSBR(
        cti_file_path,
        rmg_model_folder,
        temperatures=[528],
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

    # run to SS
    sbr_ss.run_reactor_ss()


if __name__ == "__main__":
    # execute only if run as a script
    run_sbr_test()
