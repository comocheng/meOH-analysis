exec="/work/westgroup/ChrisB/tchem_system/install/csplib/example/index_class/run_index_ODE_TChem.exe"

inputs="/work/westgroup/ChrisB/meoh-synthesis_RMG/meOH-synthesis/base/chemkin/"
chemfile=$inputs"chem-gas.inp"
thermfile=$inputs"therm.dat"
chemSurffile=$inputs"chem-surface.inp"
thermSurffile=$inputs"thermSurf.dat"
inputfile="../RMG_wDeut/6fccd2_restrict_CO2XX_and_H/IdealGasConstPressureReactor/500/results/csp/Spinning_basket_area_0.0231_energy_off_temp_500_h2_0_5_COCO2_0_5.csv"
rtol=1e-3
atol=1e-13
Area=0.00053
Pcat=0.1

$exec --Area=$Area --Pcat=$Pcat --inputfile=$inputfile --rtol=$rtol --atol=$atol --chemfile=$chemfile --thermfile=$thermfile --chemSurffile=$chemSurffile --thermSurffile=$thermSurffile
