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


import cantera as ct
import math


class MinSBR:
    """
    A minimal stirred batch reactor
    Takes away some of the options from the regular sbr class
    """
    def __init__(
        self,
        yaml_file,
        temperature=528,
        pressure=75,
        volume_flow=7.84e-06,  # [m^3/s]
        x_H2=0.75,  # responsibility is on caller to get these right
        x_CO2=0.10,
        x_CO=0.10,
        x_H2O=0.05,
        catalyst_weight=4.24e-3,
        rtol=1.0e-11,
        atol=1.0e-22,
        sensatol=1e-6,
        sensrtol=1e-6,
        reactor_type=1,
        energy="off",
        reactime=1e5,
        timestep=0.1,
        meoh_tof=0,
        h2o_tof=0,
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
        reactime = float, total time to run the reactor (for transient simulation)
        timestep = float, step taken for reactor simulation (for transient simulation)
        catalyst_weights = list of float, weight of the catalyst (in kg)
        """

        self.temperature = temperature
        self.pressure = pressure* ct.one_atm # cantera input is in pascals, so convert
        self.volume_flow = volume_flow
        self.x_H2 = x_H2
        self.x_CO2 = x_CO2
        self.x_CO = x_CO
        self.x_H2O = x_H2O
        self.catalyst_weight = catalyst_weight  # [kg]
        
        # load the experimental TOFs
        self.graaf_meoh_tof = meoh_tof
        self.graaf_h2o_tof = h2o_tof

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
        self.cat_site_per_wt = 5*61.67*1e-6*1e3 # [mol/kg] 1e-6mol/micromole, 1000g/kg
        self.cat_area = (self.catalyst_weight * self.cat_site_per_wt)/self.site_density  # [m^3]
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
        FC_temp = 293.15
        self.molar_flow = self.volume_flow * ct.one_atm / (8.3145 * FC_temp)  # [mol/s]
        self.mass_flow = self.molar_flow * (
            self.x_CO * MW_CO + self.x_CO2 * MW_CO2 + self.x_H2 * MW_H2 + self.x_H2O * MW_H2O
        )  # [kg/s]
        self.mfc = ct.MassFlowController(self.inlet, self.r, mdot=self.mass_flow)

        # A PressureController has a baseline mass flow rate matching the 'master'
        # MassFlowController, with an additional pressure-dependent term. By explicitly
        # including the upstream mass flow rate, the pressure is kept constant without
        # needing to use a large value for 'K', which can introduce undesired stiffness.
        self.outlet_mfc = ct.PressureController(self.r, self.exhaust, master=self.mfc, K=0.01)
        
        # turn off all reactions containing species x
        for spec in self.surf.species_names:
            if spec.startswith("H2*(") or spec.startswith("H2X("):
                off_species = spec
        for index,rxn in enumerate(self.surf.reactions()):
            if off_species in rxn.equation:
                self.surf.set_multiplier(0,index)
                
        rxn_strings  = [
            'CO*(14) + O*(11) <=> CO2*(15) + X(1)',
            'CO2*(15) + HCO*(16) <=> CO*(14) + HCOO*(17)',
            'CO2(4) + 2 X(1) <=> CO*(14) + O*(11)',
        ]
        for rxn_string in rxn_strings:
            rxn_index = self.surf.reaction_equations().index(rxn_string)
            self.surf.set_multiplier(0, rxn_index)

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

    def alter_thermo(self, species, dH):
        """
        perturb species NASA by dh. 
        species - name of the species
        dh - enthalpy perturbation for 5th and 13th NASA coeff
        """
        if species in self.gas.species_names:
            phase = self.gas
        else: 
            phase = self.surf
        S = phase.species(species)
        st_orig = S.thermo
        coeffs = st_orig.coeffs
        coeffs[[6, 13]] += dH / ct.gas_constant
        s_new = ct.NasaPoly2(st_orig.min_temp, st_orig.max_temp, st_orig.reference_pressure, coeffs)
        S.thermo = s_new

        phase.modify_species(phase.species_index(species), S)
        print(S.thermo.h(self.gas.T)/1000**2/96)

    def alter_reaction(
        self,
        reaction, 
        updated_rxn,
        ):
        """
        update ct reacion object with new arrhenius coefficients,
        reaction - ct reaction string
        updated_rxn - updated ct reaction object
        """
        if reaction in self.gas.reaction_equations:
            phase = self.gas
        else: 
            phase = self.surf

        rxn_index = phase.reaction_equations().index(reaction)
        phase.reaction()

    def run_reactor_ss_memory(self):
        """
        Run single reactor to steady state and save the results to an ordered dictionary in memory
        """

        gas_ROP_str = [i + " ROP [kmol/m^3 s]" for i in self.gas.species_names]

        # surface ROP reports gas and surface ROP.
        gas_surf_ROP_str = [i + " surface ROP [kmol/m^2 s]" for i in self.gas.species_names]
        surf_surf_ROP_str = [i + " ROP [kmol/m^2 s]" for i in self.surf.species_names]
        surf_ROP_str =  gas_surf_ROP_str + surf_surf_ROP_str

        gasrxn_ROP_str = [i + " ROP [kmol/m^3 s]" for i in self.gas.reaction_equations()]
        surfrxn_ROP_str = [i + " ROP [kmol/m^2 s]" for i in self.surf.reaction_equations()]


        # run the simulation
        self.sim.advance_to_steady_state()

        results = {}
        results['time (s)'] = self.sim.time
        results['T (K)'] = self.temperature
        results['P (Pa)'] = self.gas.P
        results['V (m^3/s)'] = self.volume_flow
        results['x_CO initial'] = self.x_CO
        results['x_CO2 initial'] = self.x_CO2
        results['x_H2 initial'] = self.x_H2
        results['x_H2O initial'] = self.x_H2O
        results['CO2/(CO2+CO)'] = self.CO2_ratio
        results['(CO+CO2/H2)'] = self.H2_ratio
        results['T (K) final'] = self.gas.T
        results['Rtol'] = self.sim.rtol
        results['Atol'] = self.sim.atol
        results['reactor type'] = self.reactor_type_str
        results['energy on?'] = self.energy
        results['catalyst weight (kg)'] = self.catalyst_weight
        results['graaf MeOH TOF 1/s'] = self.graaf_meoh_tof 
        results['graaf H2O TOF 1/s'] = self.graaf_h2o_tof
        results['RMG MeOH TOF 1/s'] = self.surf.net_production_rates[self.gas.species_index("CH3OH(8)")]/self.surf.site_density
        results['RMG H2O TOF 1/s'] = self.surf.net_production_rates[self.gas.species_index("H2O(5)")]/self.surf.site_density
        results['error squared MeOH TOF'] = (results['graaf MeOH TOF 1/s'] - results['RMG MeOH TOF 1/s'])**2
        results['error squared H2O TOF'] = (results['graaf H2O TOF 1/s'] - results['RMG H2O TOF 1/s'])**2
        
        for i in range(0, len(self.gas.X)):
            results[self.gas.species_names[i]] = self.gas.X[i]
        for i in range(0, len(self.surf.X)):
            results[self.surf.species_names[i]] = self.surf.X[i]
        

        # # TODO debug the fact that surf.net_production_rates also includes the gas phase
        # if len(self.gas.species_names) != len(self.gas.net_production_rates):
        #     raise ValueError('Gas species production rates do not match self.gas.net_production_rates')
        
        # # production rate of surface species contains the "surface production"
        # # i.e. desorption/adsorption reactions for gas phase species. 
        # if len(self.surf.species_names+self.gas.species_names) != len(self.surf.net_production_rates):
        #     raise ValueError('Surface species production rates do not match self.surf.net_production_rates')
        # if len(self.gas.reaction_equations()) != len(self.gas.net_rates_of_progress):
        #     raise ValueError('Gas net production rates does not match number of reaction equations')
        # if len(self.surf.reaction_equations()) != len(self.surf.net_rates_of_progress):
        #     raise ValueError('Surface net production rates does not match number of reaction equations')

        # Enter the ROP's
        for i in range(0, len(self.gas.net_production_rates)):
            results[gas_ROP_str[i]] = self.gas.net_production_rates[i]

        for i in range(0, len(self.surf.net_production_rates)):
            results[surf_ROP_str[i]] = self.surf.net_production_rates[i]

        for i in range(0, len(self.surf.net_rates_of_progress)):
            results[surfrxn_ROP_str[i]] = self.surf.net_rates_of_progress[i]

        for i in range(0, len(self.gas.net_rates_of_progress)):
            results[gasrxn_ROP_str[i]] = self.gas.net_rates_of_progress[i]

        return results

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
