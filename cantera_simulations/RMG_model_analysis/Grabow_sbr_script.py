import cantera as ct
import logging
import sys
import sbr

# filepath for writing files
if len(sys.argv) == 3:
    # Pass the cti file and rmg file as an argument to the script
    cti_file = sys.argv[1]
    rmg_model_folder  = sys.argv[2]
else:
    rmg_model_folder  = "../../../meOH-synthesis/"
    cti_file = rmg_model_folder  + "base/cantera/chem_annotated.cti"

# Reactor settings arrays for run
Temps = [400, 500, 528, 600]

# pressure and volume flow are used in Graaf, but to limit 
# the number of runs we will use the ones used for the plots
Pressures = [75]
# volume_flows = [0.00424]
volume_flows = [7.84e-06] # from grabow CO2/CO ratio study

# H2 mole fraction
H2_fraction = [0.5, 0.8, 0.75, 0.95]

# CO2/(CO+CO2)
step = 0.05
CO_CO2_ratio = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
# CO_CO2_ratio = np.arange(0,1.0+step, step)


energy = "off"

# reactor conditions
reactime = 1e4 # time to run reactor, not residence time
reactor_type = 1 # ideal gas reactor
energy = "off"

# rtol=1.0e-11
# atol=1.0e-22
rtol=1.0e-13
atol=1.0e-24


# sensitivity settings
sensitivity = 0 # 1 = kinetic, 2 = thermo, 3 = both, 0 = none
sensatol = 1e-6
sensrtol = 1e-6
# sens_species = ["H2O(5)", "CH3OH(8)", "CO(3)", "CO2(4)", "H2(2)"]
sens_species = ["H2O(5)"]
# use grabow model
grabow = False

solved = False

while not solved:  
    try: 
        # run to SS. if solver fails, change rtol and atol. 
        # initialize reactor
        sbr_ss = sbr.sbr(
            cti_file,
            rmg_model_folder, 
            t_array=Temps,
            p_array=Pressures,
            v_array=volume_flows,
            h2_array=H2_fraction,
            co2_array=CO_CO2_ratio,
            rtol=rtol,
            atol=atol,
            reactor_type=reactor_type,
            energy=energy,
            sensitivity=sensitivity,
            sensatol=sensatol,
            sensrtol=sensrtol,
            reactime=reactime,
            grabow=grabow,
            sens_species=sens_species,
            fluxes = True,
            csp = True,
        )
        sbr_ss.run_reactor_ss()
        print('sim completed')
        solved = True
    except ct.CanteraError as e:
        logging.error(f"{e}\n",
                    "Steady state run failed, conditions:"
                    f"atol: {atol}", 
                    f"rtol: {rtol}", 
                    f"sensitivity atol: {sensatol}",
                    f"sensitivity rtol: {sensrtol}",
                    f"Temp: {sbr_ss.gas.T}",
                    f"Pressure: {sbr_ss.gas.P}",
                    f"time: {sbr_ss.sim.time}",
                    f"H2: {sbr_ss.X_h2}"
                    f"CO2: {sbr_ss.X_co2}",
                    f"CO: {sbr_ss.X_co}",
                    f"H2O: {sbr_ss.X_h2o}",
                    )
        if sensatol > atol:       
            sensatol = sensatol/10
        else: 
            atol = atol/10

        if sensrtol > rtol:       
            sensrtol = sensatol/10
        else: 
            rtol = rtol/10

        if atol < 1E-30:
            raise Exception('giving up, atol too small')
        pass

