import cantera as ct
import numpy as np
import copy

class sensitivity_test:

    def __init__(
        self, 
        yaml_file,
        rtol=1.0e-11,
        atol=1.0e-22,
        energy="on",
        sensatol=1e-6,
        sensrtol=1e-6,
        reactime=1e5,
        ):

        """
        initialize object
        """
        self.yaml_file = yaml_file
        self.gas = ct.Solution(yaml_file)
        self.surf = ct.Interface(yaml_file,"Pt_surf", [self.gas])

        self.r = ct.IdealGasReactor(self.gas, energy=energy)
        # set initial temps, pressures, concentrations
        self.temp = 400              # kelvin
        self.pressure =  ct.one_atm  # Pascals

        # mole fractions: 
        self.concentrations = {"CH4": 0.095, "O2": 0.21, "AR": 0.79}

        # initialize T and P
        self.gas.TPX = self.temp, self.pressure, self.concentrations
        self.surf.TP = self.temp, self.pressure

        # create gas inlet
        self.inlet = ct.Reservoir(self.gas)

        # create gas outlet
        self.exhaust = ct.Reservoir(self.gas)

        # Reactor volume
        self.rradius = 35e-3
        self.rlength = 70e-3
        self.rvol = (self.rradius ** 2) * np.pi * self.rlength/2

        # Catalyst Surface Area
        self.site_density = (
            self.surf.site_density * 1000
        )  # [mol/m^2]cantera uses kmol/m^2, convert to mol/m^2
        self.cat_weight = 4.24e-3  # [kg]
        self.cat_site_per_wt = (300 * 1e-6) * 1000  # [mol/kg] 1e-6mol/micromole, 1000g/kg
        self.cat_area = self.site_density / (self.cat_weight * self.cat_site_per_wt)  # [m^3]

        # reactor initialization
        self.r = ct.IdealGasReactor(self.gas, energy=energy)

        # calculate the available catalyst area in a differential reactor
        self.rsurf = ct.ReactorSurface(self.surf, self.r, A=self.cat_area)
        self.r.volume = self.rvol
        self.surf.coverages = "PT(S):1.0"

        # initialize reactor network
        self.sim = ct.ReactorNet([self.r])

        # set relative and absolute tolerances on the simulation
        self.sim.rtol = rtol
        self.sim.atol = atol

    def brute_unperturbed(self, species):
        S = self.surf.species(species)
        index = self.surf.species_index(species)
        st = S.thermo
        coeffs = st.coeffs

        # get initial H298
        Ho = S.thermo.h(298)
        print(f"initial H298: {Ho}")

        # run simulation to a short time and get concentration
        self.sim.advance(1e-8)

        # second advance. this will have to change to reflect
        # the step from the "normal_thermo" function
        self.sim.step()
        self.sim.advance(1e-8 + 2.2349856286302647e-10)

        # # get unperturbed value for Concentration
        to = self.sim.time
        Co = self.surf.concentrations[index]
        print(f"initial conc: {Co}")
        return [Co,Ho,to]

    def brute_perturbed(self, species, dH):
        S = self.surf.species(species)
        index = self.surf.species_index(species)
        st = S.thermo
        coeffs = st.coeffs

        # perturb coefficients for new thermo
        coeffs[[6, 13]] += dH / ct.gas_constant
        snew = ct.NasaPoly2(st.min_temp, st.max_temp, st.reference_pressure, coeffs)
        S.thermo = snew
        self.surf.modify_species(self.surf.species_index(species), S)

        # get new H298
        Hf = S.thermo.h(298)
        print(f"final H298: {Hf}")
        self.sim.advance(1e-8)
        self.sim.step()
        # self.sim.advance(1e-8 + 2.2349856286302647e-10)

        tf = self.sim.time
        Cf = self.surf.concentrations[index]
        print(f"final conc: {Cf}")

        # sensitivity1 = (Ho/Co)*(Co-Cf)/(Ho-Hf)
        # sensitivity2 = (Hf/Cf)*(Co-Cf)/(Ho-Hf)
        # print(sensitivity1)
        # print(sensitivity2)
        return [Cf,Hf,tf]

    def normal_thermo(self, species):
        # get index of species in solution. 
        index = self.surf.species_index(species)
        S = self.surf.species(species)

        # print initial Conc for troubleshooting (should be 0)
        Co = self.surf.concentrations[index]
        print(f"initial conc: {Co}")

        # get initial H298, make sure it is properly reset
        Ho = S.thermo.h(298)
        print(f"initial H298: {Ho}")

        # add sensitivity using our new function
        self.rsurf.add_sensitivity_species_enthalpy(index)
        
        # advance the simulation forward so we have a non zero start 
        # record initial concentration for troubleshooting
        self.sim.advance(1e-8)
        to = self.sim.time
        print(f"simulation time before step:{to}")
        Ci = self.surf.concentrations[index]
        print(f"n-1 conc: {Ci}")
        
        # step the simulation forward once
        self.sim.step()
        tf = self.sim.time
        print(f"simulation time after step:{tf}")
        print(f"amount of time for step: {tf-to}")

        # record final concentration for troubleshooting
        Cf = self.surf.concentrations[index]
        print(f"final conc: {Cf}")

        # get final sensitivity value
        return self.sim.sensitivities()[index]



