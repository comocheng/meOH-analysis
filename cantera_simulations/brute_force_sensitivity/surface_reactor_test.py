import cantera as ct
import numpy as np
import copy
import sys, getopt

class sensitivity_test:

    def __init__(
        self,
        yaml_file,
        rtol=1.0e-9,
        atol=1.0e-18,
        sensatol=1e-6,
        sensrtol=1e-6,
        reactime=1e5,
        ):

        """
        initialize object
        """
        self.yaml_file = yaml_file
        self.gas = ct.Solution(yaml_file)
        self.surf = ct.Interface(yaml_file,"surface1", [self.gas])

        # set initial temps, pressures, concentrations
        self.temp = 700              # kelvin
        self.pressure =  ct.one_atm  # Pascals

        # mole fractions:
        self.concentrations = {"CH4(2)": 0.148, "O2(3)": 0.296, "AR": 0.556}

        # Reactor volume
        self.rradius = 8.25e-3
        self.rlength = 70e-3
        self.porosity = 0.81
        self.area = (self.rradius ** 2) * np.pi
        # Catalyst Surface Area
        self.site_density = (
            self.surf.site_density * 1000
        )  # [mol/m^2]cantera uses kmol/m^2, convert to mol/m^2
        # self.cat_weight = 4.24e-3  # [kg]
        # self.cat_site_per_wt = (300 * 1e-6) * 1000  # [mol/kg] 1e-6mol/micromole, 1000g/kg
        # self.cat_area = self.site_density / (self.cat_weight * self.cat_site_per_wt)  # [m^3]
        self.cat_area_per_vol = 1.6e4  # m2/m3, which is 160 cm2/cm3, as used in Horn 2006
        self.tot_flow = 0.208  # constant inlet flow rate in mol/min, equivalent to 4.7 slpm
        self.flow_rate = 4.7 * .001 / 60  # m^3/s, as seen in as seen in Horn 2007
        self.velocity = self.flow_rate / self.area  # m/s

        self.on_catalyst = 1000  # catalyst length 10mm, from Ref 17
        self.off_catalyst = 2000
        self.dt = 1.0

        # new sensitivities
        # length = 510 * mm  # Reactor length - m
        # N_reactors = 51001
        # on_catalyst = 1000  # catalyst length 10mm, from Ref 17
        # off_catalyst = 51000
        self.n_reactors = 7001
        self.reactor_len = self.rlength/(self.n_reactors-1)
        self.rvol = self.area * self.reactor_len * self.porosity

        # catalyst area in one reactor
        self.cat_area = self.cat_area_per_vol * self.rvol

        # set relative and absolute tolerances on the simulation
        self.rtol = rtol
        self.atol = atol
        self.sens_atol = sensatol
        self.sens_rtol = sensrtol

    def set_simulation(self, species):
        self.gas.TPX = 273.15, self.pressure, self.concentrations
        mass_flow_rate = self.flow_rate * self.gas.density_mass
        self.gas.TPX = self.temp, self.pressure, self.concentrations
        self.surf.TP = self.temp, self.pressure

        self.surf.coverages = "X(1):1.0"

        idx = self.surf.species_index(species)

        # create a reactor
        r = ct.IdealGasReactor(self.gas)

        # set the volume of the reactor
        r.volume = self.rvol

        # create a reservoir to represent the reactor immediately upstream. Note
        # that the gas object is set already to the state of the upstream reactor
        upstream = ct.Reservoir(self.gas, name='upstream')

        # create a reservoir for the reactor to exhaust into. The composition of
        # this reservoir is irrelevant.
        downstream = ct.Reservoir(self.gas, name='downstream')

        rsurf = ct.ReactorSurface(self.surf, r, A=self.cat_area)

        # The mass flow rate into the reactor will be fixed by using a
        # MassFlowController object.
        m = ct.MassFlowController(upstream, r, mdot=mass_flow_rate)

        # We need an outlet to the downstream reservoir. This will determine the
        # pressure in the reactor. The value of K will only affect the transient
        # pressure difference.
        v = ct.PressureController(r, downstream, master=m, K=1e-5)

        sim = ct.ReactorNet([r])
        sim.max_err_test_fails = 12

        rsurf.add_sensitivity_species_enthalpy(idx)
        sim.rtol = self.rtol
        sim.atol = self.atol
        sim.rtol_sensitivity = self.sens_rtol
        sim.atol_sensitivity = self.sens_atol

        gas_out = []
        surf_out = []
        dist_array = []
        T_array = []
        sens_array = []

        for n in range(self.n_reactors):
            self.gas.TDY = r.thermo.TDY
            upstream.syncState()
            if n == self.on_catalyst:
                self.surf.set_multiplier(1.0)
            if n == self.off_catalyst:
                self.surf.set_multiplier(0.0)
            sim.reinitialize()
            sim.advance_to_steady_state()
            dist = n * self.reactor_len * 1e3 #distance in mm
            dist_array.append(dist)
            T_array.append(self.surf.T)
            sens_array.append(sim.sensitivities()[idx])
            kmole_flow_rate = mass_flow_rate / self.gas.mean_molecular_weight  # kmol/s
            gas_out.append(1000 * 60 * kmole_flow_rate * self.gas.X.copy())  # molar flow rate in moles/minute
            surf_out.append(self.surf.X.copy())

            if n >= 1001:
                if np.max(abs(np.subtract(gas_out[-2], gas_out[-1]))) < 1e-15:
                    break

        gas_out = np.array(gas_out)
        surf_out = np.array(surf_out)
        gas_names = np.array(self.gas.species_names)
        surf_names = np.array(self.surf.species_names)
        data_out = gas_out, surf_out, gas_names, surf_names, dist_array, T_array, sens_array
        print(len(dist_array))
        return data_out

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

    def brute_perturbed(self, species):
        S = self.surf.species(species)
        index = self.surf.species_index(species)
        st = S.thermo
        coeffs = st.coeffs
        dH = st.h(self.temp) * 0.01

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
        # get index of species in solution
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

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print('surface_reactor_test.py -i <inputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('surface_reactor_test.py -i <inputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
            return inputfile

if __name__ == "__main__":
    yaml_file = main(sys.argv[1:])
    sens = sensitivity_test(yaml_file)
    sens_value = sens.set_simulation("HX(21)")
    print(sens_value)
    # sens_value_bf = sens.brute_perturbed("COX(23)") - sens.brute_unperturbed("COX(23)") / 0.01
    # print(sens_value, sens_value_bf)
