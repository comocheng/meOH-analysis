# Data sources
database(
    thermoLibraries=['surfaceThermoPt111', 'primaryThermoLibrary'], #'thermo_DFT_CCSDTF12_BAC','DFT_QCI_thermo'], # 'surfaceThermoPt' is the default. Thermo data is derived using bindingEnergies for other metals 
    reactionLibraries = [], # ('Surface/CPOX_Pt/Deutschmann2006', False) when Ni is used change the library to Surface/Deutschmann_Ni 
    seedMechanisms = [],
    kineticsDepositories = ['training'],
    kineticsFamilies = ['surface','default'],
    kineticsEstimator = 'rate rules',
)

catalystProperties(
    bindingEnergies = {  # for Cu(111)    
                          'H':(-2.584, 'eV/molecule'),
                          'C':(-4.960, 'eV/molecule'),
                          'N':(-3.584, 'eV/molecule'),
                          'O':(-4.208, 'eV/molecule')
                      },
    surfaceSiteDensity=(2.72e-9, 'mol/cm^2'), # Default for Pt(111)
)
    
species(
    label='CH3OH',
    reactive=True,
    structure=adjacencyList(
       """
1 C u0 p0 c0 {2,S} {3,S} {4,S} {5,S}
2 O u0 p2 c0 {1,S} {6,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 H u0 p0 c0 {1,S}
6 H u0 p0 c0 {2,S}
"""),
)

species(
   label='CO2',
   reactive=True,
   structure=adjacencyList(
       """
1 O u0 p2 c0 {2,D}
2 C u0 p0 c0 {1,D} {3,D}
3 O u0 p2 c0 {2,D}
"""),
)

species(
   label='CO',
   reactive=True,
   structure=adjacencyList(
       """
1 C u0 p1 c-1 {2,T}
2 O u0 p1 c+1 {1,T}
"""),
)

species(
   label='H2',
   reactive=True,
   structure=adjacencyList(
       """
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
"""),
)

species(
   label='H2O',
   reactive=True,
   structure=adjacencyList(
       """
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
"""),
)

species(
    label='N2',
    reactive=False,
    structure=SMILES("N#N"),
)

species(
    label='vacantX',
    reactive=True,
    structure=adjacencyList("1 X u0"),
)

#----------
# Reaction systems
surfaceReactor(
    temperature=(553,'K'),
    initialPressure=(40, 'bar'),
    initialGasMoleFractions={
        "CO": 60,
        "H2": 30,
        "N2": 5,
    },
    initialSurfaceCoverages={
        "vacantX": 1.0,
    },
    surfaceVolumeRatio=(1.e5, 'm^-1'),
    #terminationTime=(42300000, 's'),
    terminationConversion = { "CO":0.90,},
)


simulator(
    atol=1e-9,#18,
    rtol=1e-6,#12,
)

model(
    toleranceKeepInEdge=0.0,
    toleranceMoveToCore=1e-5,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100,
)

options(
    units='si',
    generateOutputHTML=True,
    generatePlots=False, # Enable to make plots of core and edge size etc. But takes a lot of the total runtime!
    saveEdgeSpecies=False,
    saveSimulationProfiles=True,
)
