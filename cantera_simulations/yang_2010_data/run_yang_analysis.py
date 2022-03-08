import cantera as ct
import numpy as np
import pandas as pd
import math 
import matplotlib.pyplot as plt


class yang_reactor():
    def __init__(self, *yamls): 
        '''
        create yang batch reactor. if more than 1 input file is specified, 
        load the first file to as a gas phase, and the second file as 
        the surface phase 
        '''

        # Load in values from plot above
        # that we already know 
        temp = [525, 550, 575, 600]

        meoh_ln_rate = [
            -6.691144708,
            -5.978401728,
            -4.48812095,
            -3.894168467,
        ]

        rwgs_ln_rate = [
            -0.578342066,
            0.572607525,
            1.171517945,
            2.072487534,
        ]

        # convert to molecules/cm^2/sec
        meoh_rates_cm = np.exp(meoh_ln_rate)*10**15
        meoh_rates = dict(zip(temp, meoh_rates_cm))

        rwgs_rates_cm = np.exp(rwgs_ln_rate)*10**15
        rwgs_rates = dict(zip(temp, rwgs_rates_cm))
        "{:2e}".format(meoh_rates_cm[525])

        cat_area = 1e2

        # convert to pascals
        p_co2 = 0.5 * ct.one_atm
        p_h2 =  4.5 * ct.one_atm
        p_total = p_co2+p_h2

        # get total pressure at temp using ig law pv = nrt
        initial_temp = 300 #[k]
        p_total_at_temp = np.array(temp)*p_total/initial_temp

        x_co2 = p_co2/p_total
        x_h2 = p_h2/p_total

        # create thermo phases
        if len(yamls) == 1:
            self.gas = ct.Solution(yamls[0], "gas")
            self.surf = ct.Interface(yamls[0],"surface1", [gas])

        elif len(yamls) == 2: 
            self.gas = ct.Solution(yamls[0], "gas")
            self.surf = ct.Interface(yamls[1],"surface1", [gas])
            
            for name in self.gas.species_names:
                if name.startswith("CO2("):
                    co2_str = name
                if name.startswith("H2("):
                    h2_str = name
                if name.startswith("CH3OH("):
                    meoh_str = name

                    
        mole_fractions = {co2_str:x_co2,h2_str:x_h2}

        # molecular weights for mass flow calculations
        mw_co = 28.01e-3  # [kg/mol]
        mw_co2 = 44.01e-3  # [kg/mol]
        mw_h2 = 2.016e-3  # [kg/mol]
        mw_h2o = 18.01528e-3  # [kg/mol]

        # convert to pascals
        p_co2 = 0.5 * ct.one_atm
        p_h2 =  4.5 * ct.one_atm
        p_total = p_co2+p_h2

        # get total pressure at temp using ig law pv = nrt
        initial_temp = 300 #[k]
        p_total_at_temp = np.array(temp)*p_total/initial_temp

        x_co2 = p_co2/p_total
        x_h2 = p_h2/p_total

        print(x_h2, x_co2, p_total)
        print(temp, p_total_at_temp)

        self.gas.TPX = temp[0], p_total_at_temp[0], mole_fractions
        self.surf.TP = temp[0], p_total_at_temp[0]

        self.r = ct.IdealGasReactor(self.gas, energy="off")
        self.rsurf = ct.ReactorSurface(surf, r, A=cat_area)
        # initialize reactor network
        self.sim = ct.ReactorNet([r])
        self.sim.atol = 1e-16
        self.sim.rtol = 1e-14   
    
    def run():
        '''
        performs a run at all temperatures specified by "temps"
        '''

        reactime = 200*60
        t = 0
        dt = 5.0
        conversion = 0 

        nspecies = len(self.gas.species_names)
        meoh_gas_index = self.gas.kinetics_species_index(meoh_str)
        meoh_surf_index = self.surf.kinetics_species_index(meoh_str)

        meoh_rop = []  
        meoh_moles_norm = []
        temps = []
        press = []
        times = []
        conversions = []
        self.r.volume = 0.1

        # see if the Rop is constant
        while t < reactime and conversion < 0.05:
            t += dt
            self.sim.advance(t)
            meoh_total_rop = self.gas.net_production_rates[meoh_gas_index] + self.surf.net_production_rates[meoh_surf_index]
            meoh_total_rop = meoh_total_rop/self.surf.site_density
            meoh_rop.append(meoh_total_rop)
            times.append(sim.time) # time in minutes
            temps.append(self.gas.T)
            press.append(self.gas.P)
            
            moles_co2_initial = p_co2*self.r.volume/(ct.gas_constant*initial_temp) # Pco2*v/RT
            moles_meoh = (self.gas[meoh_str].X*self.gas.P*self.r.volume)/(ct.gas_constant*initial_temp) # mole frac*total moles
            
            # calculate conversion:
            # (moles meoh possible (starting moles co2)-moles meoh current)/moles co2 initial
            conversion = 1 - (moles_co2_initial - (moles_meoh))/moles_co2_initial
            conversions.append(conversion)
            
            # calculate the moles methanol normalized
            # to the number of surface sites
            meoh_moles_norm.append(float(moles_meoh/(self.surf.site_density*cat_area))) 
        