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

# Reactor settings arrays for run
Temps = [400]

# pressure and volume flow are used in Graaf, but to limit 
# the number of runs we will use the ones used for the plots
# Pressures = [15, 30, 50, 75]
# volume_flows = [0.00424, 0.0106, 0.02544]
Pressures = [75]
volume_flows = [0.00424]

# Mole fractions from runs for reference (min max and mid from Graaf data)
# X_cos = [0.053, 0.1365, 0.220]
# X_co2s = [0.261,0.1205, 0.261]
# X_h2s = [0.625, 0.7625, 0.900]

# CO+CO2/H2
# H2_fraction = [0.8,0.5,0.95,0.75]
H2_fraction = [0.5]

# CO2/CO
CO_CO2_ratio = [0.1]



# reactor conditions
reactime = 1e4 # time to run reactor, not residence time
reactor_type = 1 # ideal gas reactor
energy = "off"
rtol=1.0e-11
atol=1.0e-22


# sensitivity settings
sensitivity = False
sensatol = 1e-6
sensrtol = 1e-6

# use grabow model
grabow = False


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
)


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
    )
    # run to time specified by reactime
    sbr_tran.run_reactor_transient()

except ct.CanteraError as e:
    print("transient run failed",
          f"atol: {atol}", 
          f"rtol: {rtol}", 
          f"sensitivity atol: {sensatol}",
          f"sensitivity rtol: {sensrtol}",
          f"Temp: {sbr_tran.gas.T}",
          f"Pressure: {sbr_tran.gas.P}",
          f"time: {sbr_tran.sim.time}",
          f"H2: {sbr_tran.X_h2}"
          f"CO2: {sbr_tran.X_co2}",
          f"CO: {sbr_tran.X_co}",
          f"H2O: {sbr_tran.X_h2o}",
         )