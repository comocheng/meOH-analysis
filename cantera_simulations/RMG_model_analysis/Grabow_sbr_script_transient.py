import logging
import sys
import sbr
import cantera as ct

# filepath for writing files
if len(sys.argv) == 3:
    # Pass the cti file and rmg file as an argument to the script
    cti_file = sys.argv[1]
    rmg_model_folder  = sys.argv[2]
else:
    rmg_model_folder  = "../../../meOH-synthesis/"
    cti_file = rmg_model_folder  + "base/cantera/chem_annotated.cti"

# Reactor temp settings arrays for run
Temps = [600, 500, 400]

# pressure and volume flow are used in Graaf, but to limit 
# the number of runs we will use the ones used for the plots
Pressures = [75]
volume_flows = [0.00424]

# H2 mole fraction
H2_fraction = [0.8,0.5,0.95,0.75]

# CO2/(CO+CO2)
CO_CO2_ratio = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

# reactor conditions
reactime = 1e2 # time to run reactor, not residence time
reactor_type = 1 # ideal gas reactor
energy = "off"
rtol=1.0e-11
atol=1.0e-22
timestep = 1

# sensitivity settings
sensitivity = 2
sensatol = 1e-6
sensrtol = 1e-6

# use grabow model
grabow = False

solved = False

while not solved:  
    try: 
    # initialize reactor
        sbr_tran = sbr.sbr(
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
            graaf=False,
            timestep=timestep,
        )
        # run to time specified by reactime
        sbr_tran.run_reactor_transient()
        print("run completed. Good Job!!!")
        solved = True

    # if run fails, set to a little before where it failed. 
    # proceed to divide reactime and timestep by 10 if this fails again
    except ct.CanteraError as e:
            logging.error(f"Steady state run failed, conditions:\n",
                        f"atol: {atol}\n", 
                        f"rtol: {rtol}\n", 
                        f"sensitivity atol: {sensatol}\n",
                        f"sensitivity rtol: {sensrtol}\n",
                        f"Temp: {sbr_tran.gas.T}\n",
                        f"Pressure: {sbr_tran.gas.P}\n",
                        f"time: {sbr_tran.sim.time}\n",
                        f"H2: {sbr_tran.X_h2}\n"
                        f"CO2: {sbr_tran.X_co2}\n",
                        f"CO: {sbr_tran.X_co}\n",
                        f"H2O: {sbr_tran.X_h2o}\n",
                        )

            if reactime > sbr_tran.sim.time:
                reactime = sbr_tran.sim.time - sbr_tran.sim.time*1e-1
            else: 
                reactime = reactime * 1e-1

            timestep = reactime*1e-2
            if reactime < 1e-20:
                raise Exception('giving up, reactime too small')
            pass


