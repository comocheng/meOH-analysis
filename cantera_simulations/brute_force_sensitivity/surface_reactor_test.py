import cantera as ct
import numpy as np
import copy
import sys, getopt

class sensitivity_test:

    def __init__(
        self,
        yaml_file,
        rtol=1.0e-10,
        atol=1.0e-20,
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
        self.surf = ct.Interface(yaml_file,"surface1", [self.gas])

        self.r = ct.IdealGasReactor(self.gas, energy=energy)
        # set initial temps, pressures, concentrations
        self.temp = 700              # kelvin
        self.pressure =  ct.one_atm  # Pascals

        # mole fractions:
        self.concentrations = {"CH4(2)": 0.148, "O2(3)": 0.296, "AR": 0.556}

        # initialize T and P
        self.gas.TPX = self.temp, self.pressure, self.concentrations
        self.surf.TP = self.temp, self.pressure

        # create gas inlet
        self.inlet = ct.Reservoir(self.gas)

        # create gas outlet
        self.exhaust = ct.Reservoir(self.gas)

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

        # reactor initialization
        self.r = ct.IdealGasReactor(self.gas, energy=energy)

        # calculate the available catalyst area in a differential reactor
        self.rsurf = ct.ReactorSurface(self.surf, self.r, A=self.cat_area)
        self.r.volume = self.rvol
        self.surf.coverages = "X(1):1.0"

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
    sens_value_bf = sens.brute_perturbed("COX(23)") - sens.brute_unperturbed("COX(23)") / 0.01
    sens_value = sens.normal_thermo("COX(23)")
    print(sens_value, sens_value_bf)
