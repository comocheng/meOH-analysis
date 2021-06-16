exec='/work/westgroup/ChrisB/tchem_system/install/csplib/example/index_class/run_index_ODE_TCSTR_TChem.exe'

inputs="/work/westgroup/ChrisB/meoh-synthesis_RMG/meOH-synthesis/base/chemkin/"
chemfile=$inputs"chem-gas.inp"
thermfile=$inputs"therm.dat"
chemSurffile=$inputs"chem-surface_no_sites.inp"
thermSurffile=$inputs"thermSurf.dat"
inputfile="/work/westgroup/ChrisB/meOH-synthesis_cantera/meOH-synthesis/cantera_simulations/RMG_wDeut/2021_5_5_1342_199260_fixed_typo_in_Cu111_library_for_a5_nasa_parameter/IdealGasReactor/energy_off/sensitivity_off/400/csp/CSP_Spinning_basket_area_0.0231_energy_off_temp_400_h2_0_5_COCO2_0_1.dat"

rtol=1e-6
atol=1e-13
Acat=2.77191e3
Vol=2.77191e-2
mdotIn=5e-3
isoThermic=true


$exec --isoThermic=$isoThermic --rtol=$rtol --atol=$atol --mdotIn=$mdotIn --Acat=$Acat --Vol=$Vol --inputfile=$inputfile --chemfile=$chemfile --thermfile=$thermfile --chemSurffile=$chemSurffile --thermSurffile=$thermSurffile

