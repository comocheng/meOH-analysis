###############################################
# Grabow reactor batch script
# Chris Blais
# Northeastern University
# runs through all reactor conditions
###############################################

import pandas as pd
import numpy as np
import time
import cantera as ct
from matplotlib import pyplot as plt
import csv
import math
import os
import sys
import re
import itertools
import logging
from collections import defaultdict
import git
import time
from IPython.display import display

from rmgpy.molecule import Molecule
from rmgpy.data.base import Database


def save_pictures(git_path="", species_path="", overwrite=False):
    """
    Save a folder full of molecule pictures, needed for the pretty dot files.

    Saves them in the results directory, in a subfolder "species_pictures".
    Unless you set overwrite=True, it'll leave alone files that are
    already there.
    """
    dictionary_filename = git_path + "/base/chemkin/species_dictionary.txt"
    specs = Database().get_species(dictionary_filename, resonance=False)

    images_dir = os.path.join(species_path)
    os.makedirs(images_dir, exist_ok=True)
    for name, species in specs.items():
        filepath = os.path.join(images_dir, name + ".png")
        if not overwrite and os.path.exists(filepath):
            continue
        species.molecule[0].draw(filepath)


def prettydot(species_path="", dotfilepath="", strip_line_labels=False):
    """
    Make a prettier version of the dot file (flux diagram)

    Assumes the species pictures are stored in a directory
    called 'species_pictures' alongside the dot file.
    """
    pictures_directory = f'{species_path}/'

    if strip_line_labels:
        print("stripping edge (line) labels")

    reSize = re.compile(r'size="5,6"\;page="5,6"')
    reNode = re.compile(
        r'(?P<node>s\d+)\ \[\ fontname="Helvetica",\ label="(?P<label>[^"]*)"\]\;'
    )

    rePicture = re.compile(r"(?P<smiles>.+?)\((?P<id>\d+)\)\.png")
    reLabel = re.compile(r"(?P<name>.+?)\((?P<id>\d+)\)$")

    species_pictures = dict()
    for picturefile in os.listdir(pictures_directory):
        match = rePicture.match(picturefile)
        if match:
            species_pictures[match.group("id")] = picturefile
        else:
            pass
            # print(picturefile, "didn't look like a picture")

    filepath = dotfilepath

    if not open(filepath).readline().startswith("digraph"):
        raise ValueError("{0} - not a digraph".format(filepath))

    infile = open(filepath)
    prettypath = filepath.replace(".dot", "", 1) + "-pretty.dot"
    outfile = open(prettypath, "w")

    for line in infile:
        (line, changed_size) = reSize.subn('size="12,12";page="12,12"', line)
        match = reNode.search(line)
        if match:
            label = match.group("label")
            idmatch = reLabel.match(label)
            if idmatch:
                idnumber = idmatch.group("id")
                if idnumber in species_pictures:
                    line = (
                        f'%s [ image="{pictures_directory}%s" label="" width="0.5" height="0.5" imagescale=true fixedsize=false color="none" ];\n'
                        % (match.group("node"), species_pictures[idnumber])
                    )

        # rankdir="LR" to make graph go left>right instead of top>bottom
        if strip_line_labels:
            line = re.sub(r'label\s*=\s*"\s*[\d.]+"', 'label=""', line)

        # change colours
        line = re.sub(r'color="0.7,\ (.*?),\ 0.9"', r'color="1.0, \1, 0.7*\1"', line)

        outfile.write(line)

    outfile.close()
    infile.close()
    print(f"Graph saved to: {prettypath}")
    os.system(f'dot {prettypath} -Tpng -o{prettypath.replace(".dot", "", 1) + ".png"} -Gdpi=300')
    return prettypath


def show_flux_diagrams(self, suffix="", embed=False):
    """
    Shows the flux diagrams in the notebook.
    Loads them from disk.
    Does not embed them, to keep the .ipynb file small,
    unless embed=True. Use embed=True if you might over-write the files,
    eg. you want to show flux at different points.
    """
    import IPython

    for element in "CHONX":
        for phase_object in (self.gas, self.surf):
            phase = phase_object.name
            img_file = (
                f"reaction_path_{element}_{phase}{'_' if suffix else ''}{suffix}.png"
            )
            display(IPython.display.HTML(f"<hr><h2>{element} {phase}</h2>"))
            if embed:
                display(IPython.display.Image(filename=img_file, width=400, embed=True))
            else:
                display(IPython.display.Image(url=img_file, width=400, embed=False))

        # Now do the combined
        img_file = f"reaction_path_mass{'_' if suffix else ''}{suffix}.png"
        display(IPython.display.HTML(f"<hr><h2>Combined mass</h2>"))
        if embed:
            display(IPython.display.Image(filename=img_file, width=400, embed=True))
        else:
            display(IPython.display.Image(url=img_file, width=400, embed=False))


def save_flux_diagrams(*phases, suffix="", timepoint="", species_path=""):
    """
    Saves the flux diagrams. The filenames have a suffix if provided,
    so you can keep them separate and not over-write.
    """
    for element in "CHONX":
        for phase_object in phases:
            phase = phase_object.name

            diagram = ct.ReactionPathDiagram(phase_object, element)
            diagram.title = f"Reaction path diagram following {element} in {phase}"
            diagram.label_threshold = 0.001

            dot_file = f"{suffix}/reaction_path_{element}_{phase}_{timepoint}.dot"
            img_file = f"{suffix}/reaction_path_{element}_{phase}_{timepoint}.png"
            dot_bin_path = (
                "/Users/blais.ch/anaconda3/pkgs/graphviz-2.40.1-hefbbd9a_2/bin/dot"
            )
            img_path = os.path.join(os.getcwd(), img_file)
            diagram.write_dot(dot_file)

            # also make a prettydot file
            prettydot(species_path, dot_file, strip_line_labels=False)

            # print(diagram.get_data())

            print(
                f"Wrote graphviz input file to '{os.path.join(os.getcwd(), dot_file)}'."
            )
            os.system(f"dot {dot_file} -Tpng -o{img_file} -Gdpi=200")
            print(f"Wrote graphviz output file to '{img_path}'.")


def run_reactor(
    cti_file,
    rmg_model_path,
    t_array=[528],
    p_array=[75],
    v_array=[0.00424],
    h2_array=[0.75],
    co2_array=[0.5],
    rtol=1.0e-11,
    atol=1.0e-22,
    reactor_type=0,
    energy="off",
    sensitivity=False,
    sensatol=1e-6,
    sensrtol=1e-6,
    reactime=1e5,
    grabow=False,
):
    try:
        array_i = int(os.getenv("SLURM_ARRAY_TASK_ID"))
    except TypeError:
        array_i = 0

    if grabow:
        # format grabow model the same as the others
        repo = git.Repo(rmg_model_path)
        date = time.localtime()
        git_date = "0000_00_00_0000"
        git_sha = '000000'
        git_msg = "Grabow model"
        git_file_string = f"{git_date}_{git_sha}_{git_msg}"
    else:
        # get git commit hash and message
        repo = git.Repo(rmg_model_path)
        date = time.localtime(repo.head.commit.committed_date)
        git_date = time.strftime("%Y_%m_%d_%H%M", date)
        git_sha = str(repo.head.commit)[0:6]
        git_msg = str(repo.head.commit.message)[0:50].replace(" ", "_").replace("'", "_").replace("\n", "")
        git_file_string = f"{git_date}_{git_sha}_{git_msg}"

    # set sensitivity string for file path name
    if sensitivity:
        sensitivity_str = "on"
    else:
        sensitivity_str = "off"

    # this should probably be outside of function
    settings = list(itertools.product(t_array, p_array, v_array, h2_array, co2_array))

    # constants
    pi = math.pi

    # set initial temps, pressures, concentrations
    temp = settings[array_i][0]  # kelvin
    temp_str = str(temp)[0:3]
    pressure = settings[array_i][1] * ct.one_atm  # Pascals

    X_h2 = settings[array_i][3]
    x_h2_str = str(X_h2)[0:3].replace(".", "_")
    x_CO_CO2_str = str(settings[array_i][4])[0:3].replace(".", "_")

    # Per Grabow experiments, add in H2O for X=0.75 H2 run
    if X_h2 == 0.75:
        X_h2o = 0.05
    else:
        X_h2o = 0

    X_co = (1 - (X_h2 + X_h2o)) * (settings[array_i][4])
    X_co2 = (1 - (X_h2 + X_h2o)) * (1 - settings[array_i][4])

    # normalize mole fractions just in case
    # X_co = X_co/(X_co+X_co2+X_h2)
    # X_co2= X_co2/(X_co+X_co2+X_h2)
    # X_h2 = X_h2/(X_co+X_co2+X_h2)

    mw_co = 28.01e-3  # [kg/mol]
    mw_co2 = 44.01e-3  # [kg/mol]
    mw_h2 = 2.016e-3  # [kg/mol]
    mw_h2o = 18.01528e-3  # [kg/mol]

    co2_ratio = X_co2 / (X_co + X_co2)
    h2_ratio = (X_co2 + X_co) / X_h2

    # CO/CO2/H2/H2: typical is
    if grabow:
        concentrations_rmg = {"CO": X_co, "CO2": X_co2, "H2": X_h2, "H2O": X_h2o, }
    else:
        concentrations_rmg = {"CO(3)": X_co, "CO2(4)": X_co2, "H2(2)": X_h2, "H2O(5)": X_h2o, }

    # initialize cantera gas and surface
    gas = ct.Solution(cti_file, "gas")
    surf = ct.Interface(cti_file, "surface1", [gas])

    # initialize T and P
    gas.TPX = temp, pressure, concentrations_rmg
    surf.TP = temp, pressure

    # if a mistake is made with the input,
    # cantera will normalize the mole fractions.
    # make sure that we are reporting/using
    # the normalized values
    if grabow:
        X_co = float(gas["CO"].X)
        X_co2 = float(gas["CO2"].X)
        X_h2 = float(gas["H2"].X)
        X_h2o = float(gas["H2O"].X)
    else:
        X_co = float(gas["CO(3)"].X)
        X_co2 = float(gas["CO2(4)"].X)
        X_h2 = float(gas["H2(2)"].X)
        X_h2o = float(gas["H2O(5)"].X)

    # create gas inlet
    inlet = ct.Reservoir(gas)

    # create gas outlet
    exhaust = ct.Reservoir(gas)

    # Reactor volume (divide by 2 per Graaf paper)
    rradius = 35e-3
    rlength = 70e-3
    rvol = (rradius ** 2) * pi * rlength / 2

    # Catalyst Surface Area
    site_density = (
        surf.site_density * 1000
    )  # [mol/m^2]cantera uses kmol/m^2, convert to mol/m^2
    cat_weight = 4.24e-3  # [kg]
    cat_site_per_wt = (300 * 1e-6) * 1000  # [mol/kg] 1e-6mol/micromole, 1000g/kg
    cat_area = (cat_weight * cat_site_per_wt) / site_density  # [m^2]

    # reactor initialization
    if reactor_type == 0:
        r = ct.Reactor(gas, energy=energy)
        reactor_type_str = "Reactor"
    elif reactor_type == 1:
        r = ct.IdealGasReactor(gas, energy=energy)
        reactor_type_str = "IdealGasReactor"
    elif reactor_type == 2:
        r = ct.ConstPressureReactor(gas, energy=energy)
        reactor_type_str = "ConstPressureReactor"
    elif reactor_type == 3:
        r = ct.IdealGasConstPressureReactor(gas, energy=energy)
        reactor_type_str = "IdealGasConstPressureReactor"

    # calculate the available catalyst area in a differential reactor
    rsurf = ct.ReactorSurface(surf, r, A=cat_area)
    r.volume = rvol
    if grabow:
        surf.coverages = "X:1.0"
    else:
        surf.coverages = "X(1):1.0"
    # flow controllers (Graaf measured flow at 293.15 and 1 atm)
    one_atm = ct.one_atm
    FC_temp = 293.15
    volume_flow = settings[array_i][2]  # [m^3/s]
    molar_flow = volume_flow * one_atm / (8.3145 * FC_temp)  # [mol/s]
    mass_flow = molar_flow * (
        X_co * mw_co + X_co2 * mw_co2 + X_h2 * mw_h2 + X_h2o * mw_h2o
    )  # [kg/s]
    mfc = ct.MassFlowController(inlet, r, mdot=mass_flow)

    # A PressureController has a baseline mass flow rate matching the 'master'
    # MassFlowController, with an additional pressure-dependent term. By explicitly
    # including the upstream mass flow rate, the pressure is kept constant without
    # needing to use a large value for 'K', which can introduce undesired stiffness.
    outlet_mfc = ct.PressureController(r, exhaust, master=mfc, K=0.01)

    # initialize reactor network
    sim = ct.ReactorNet([r])

    # set relative and absolute tolerances on the simulation
    sim.rtol = 1.0e-11
    sim.atol = 1.0e-22

    #################################################
    # Run single reactor
    #################################################

    # round numbers so they're easier to read
    # temp_str = '%s' % '%.3g' % tempn
    cat_area_str = "%s" % "%.3g" % cat_area
    species_path = (
        os.path.dirname(os.path.abspath(__file__))
        + f"/{git_file_string}/species_pictures"
    )
    print(species_path)
    results_path = (
        os.path.dirname(os.path.abspath(__file__)) + f"/{git_file_string}/{reactor_type_str}/energy_{energy}/sensitivity_{sensitivity_str}/{temp_str}/results"
    )
    results_path_csp = (
        os.path.dirname(os.path.abspath(__file__))
        + f"/{git_file_string}/{reactor_type_str}/energy_{energy}/sensitivity_{sensitivity_str}/{temp_str}/csp"
    )

    flux_path = (
        os.path.dirname(os.path.abspath(__file__))
        + f"/{git_file_string}/{reactor_type_str}/energy_{energy}/sensitivity_{sensitivity_str}/{temp_str}/flux_diagrams/{x_h2_str}/{x_CO_CO2_str}"
    )
    # create species folder for species pictures if it does not already exist
    try:
        os.makedirs(species_path, exist_ok=True)
        save_pictures(git_path=rmg_model_path, species_path=species_path)
    except OSError as error:
        print(error)

    try:
        os.makedirs(results_path, exist_ok=True)
    except OSError as error:
        print(error)

    try:
        os.makedirs(results_path_csp, exist_ok=True)
        print(results_path_csp)
    except OSError as error:
        print(error)
        print("no results for CSP saved")

    try:
        os.makedirs(flux_path, exist_ok=True)
    except OSError as error:
        print(error)

    gas_ROP_str = [i + " ROP [kmol/m^3 s]" for i in gas.species_names]

    # surface ROP reports gas and surface ROP. these values might be redundant, not sure.

    gas_surf_ROP_str = [i + " surface ROP [kmol/m^2 s]" for i in gas.species_names]
    surf_ROP_str = [i + " ROP [kmol/m^2 s]" for i in surf.species_names]

    gasrxn_ROP_str = [i + " ROP [kmol/m^3 s]" for i in gas.reaction_equations()]
    surfrxn_ROP_str = [i + " ROP [kmol/m^2 s]" for i in surf.reaction_equations()]
    output_filename = (
        results_path
        + f"/Spinning_basket_area_{cat_area_str}_energy_{energy}"
        + f"_temp_{temp}_h2_{x_h2_str}_COCO2_{x_CO_CO2_str}.csv"
    )
    output_filename_csp = (
        results_path_csp
        + f"/CSP_Spinning_basket_area_{cat_area_str}_energy_{energy}"
        + f"_temp_{temp}_h2_{x_h2_str}_COCO2_{x_CO_CO2_str}.dat"
    )
    outfile = open(output_filename, "w")
    outfile_csp = open(output_filename_csp, "w")
    writer = csv.writer(outfile)
    writer_csp = csv.writer(outfile_csp, delimiter='\t')

    logging.warning(results_path + output_filename)

    # Sensitivity atol, rtol, and strings for gas and surface reactions if selected
    # slows down script by a lot
    if sensitivity:
        sim.rtol_sensitivity = sensrtol
        sim.atol_sensitivity = sensatol
        sens_species = ["CH3OH(8)"]

        # turn on sensitive reactions/species
        for i in range(gas.n_reactions):
            r.add_sensitivity_reaction(i)

        for i in range(surf.n_reactions):
            rsurf.add_sensitivity_reaction(i)

        # for i in range(gas.n_species):
        #     r.add_sensitivity_species_enthalpy(i)

        # for i in range(surf.n_species):
        #     rsurf.add_sensitivity_species_enthalpy(i)

        for j in sens_species:
            gasrxn_sens_str = [
                j + " sensitivity to " + i for i in gas.reaction_equations()
            ]
            surfrxn_sens_str = [
                j + " sensitivity to " + i for i in surf.reaction_equations()
            ]
            # gastherm_sens_str = [j + " thermo sensitivity to " + i for i in gas.species_names]
            # surftherm_sens_str = [j + " thermo sensitivity to " + i for i in surf.species_names]
            sens_list = gasrxn_sens_str + surfrxn_sens_str  # + gastherm_sens_str

        writer.writerow(
            [
                "time (s)",
                "T (K)",
                "P (Pa)",
                "V (M^3/s)",
                "X_co initial",
                "X_co2initial",
                "X_h2 initial",
                "X_h2o initial",
                "CO2/(CO2+CO)",
                "(CO+CO2/H2)",
                "T (C) final",
                "Rtol",
                "Atol",
                "reactor type",
                "energy on?"
            ]
            + gas.species_names
            + surf.species_names
            + gas_ROP_str
            + gas_surf_ROP_str
            + surf_ROP_str
            + gasrxn_ROP_str
            + surfrxn_ROP_str
            + sens_list
        )

    else:
        writer.writerow(
            [
                "time (s)",
                "T (K)",
                "P (Pa)",
                "V (M^3/s)",
                "X_co initial",
                "X_co2 initial",
                "X_h2 initial",
                "X_h2o initial",
                "CO2/(CO2+CO)",
                "(CO+CO2/H2)",
                "T (C) final",
                "Rtol",
                "Atol",
                "reactor type",
                "energy on?"
            ]
            + gas.species_names
            + surf.species_names
            + gas_ROP_str
            + gas_surf_ROP_str
            + surf_ROP_str
            + gasrxn_ROP_str
            + surfrxn_ROP_str
        )
    # only run csp writer up to 100 seconds. todo: make built in criteria for steady state (ss achieved flag?)
    if sim.time < 100:
        writer_csp.writerow(
            ["iter", "t", "dt", "Density[kg/m3]", "Pressure[Pascal]", "Temperature[K]", ]
            + gas.species_names
            + surf.species_names
        )

    t = 0.0
    dt = 0.1
    iter_ct = 0
    # run the simulation
    first_run = True

    while t < reactime:
        # save flux diagrams at beginning of run
        if first_run is True:
            save_flux_diagrams(gas, suffix=flux_path, timepoint="beginning", species_path=species_path)
            save_flux_diagrams(surf, suffix=flux_path, timepoint="beginning", species_path=species_path)

            first_run = False
        t += dt
        sim.advance(t)
        # if t % 10 < 0.01:

        # some models have the special case where they do not have any gas
        # gas phase reactions. if this is true, skip over gas ROPs
        if len(gas.reactions()) > 0:
            gas_net_rates_of_progress = list(gas.net_rates_of_progress)
        else:
            gas_net_rates_of_progress = []

        if sensitivity:
            # get sensitivity for sensitive species i (e.g. methanol) in reaction j
            for i in sens_species:
                g_nrxn = gas.n_reactions
                s_nrxn = surf.n_reactions
                # g_nspec = gas.n_species
                # s_nspec = surf.n_species

                gas_sensitivities = [sim.sensitivity(i, j) for j in range(g_nrxn)]
                surf_sensitivities = [
                    sim.sensitivity(i, j) for j in range(g_nrxn, g_nrxn + s_nrxn)
                ]
                # gas_therm_sensitivities = [sim.sensitivity(i,j)
                # for j in range(g_nrxn+s_nrxn,g_nrxn+s_nrxn+g_nspec)]
                # surf_therm_sensitivities = [sim.sensitivity(i,j)
                # for j in range(g_nrxn+s_nrxn+g_nspec,g_nrxn+s_nrxn+g_nspec+s_nspec)]

                sensitivities_all = (
                    gas_sensitivities
                    + surf_sensitivities
                    # + gas_therm_sensitivities
                )

            writer.writerow(
                [
                    sim.time,
                    temp,
                    pressure,
                    volume_flow,
                    X_co,
                    X_co2,
                    X_h2,
                    X_h2o,
                    co2_ratio,
                    h2_ratio,
                    gas.T,
                    sim.rtol,
                    sim.atol,
                    reactor_type_str,
                    energy,
                ]
                + list(gas.X)
                + list(surf.X)
                + list(gas.net_production_rates)
                + list(surf.net_production_rates)
                + gas_net_rates_of_progress
                + list(surf.net_rates_of_progress)
                + sensitivities_all
            )

        else:
            writer.writerow(
                [
                    sim.time,
                    temp,
                    pressure,
                    volume_flow,
                    X_co,
                    X_co2,
                    X_h2,
                    X_h2o,
                    co2_ratio,
                    h2_ratio,
                    gas.T,
                    sim.rtol,
                    sim.atol,
                    reactor_type_str,
                    energy,
                ]
                + list(gas.X)
                + list(surf.X)
                + list(gas.net_production_rates)
                + list(surf.net_production_rates)
                + gas_net_rates_of_progress
                + list(surf.net_rates_of_progress)
            )

        writer_csp.writerow(
            [
                iter_ct,
                sim.time,
                dt,
                gas.density,
                gas.P,
                gas.T,
            ]
            + list(gas.X)
            + list(surf.X)
        )

        iter_ct += 1

    outfile.close()
    outfile_csp.close()

    # save flux diagrams at the end of the run
    save_flux_diagrams(gas, suffix=flux_path, timepoint="end", species_path=species_path)
    save_flux_diagrams(surf, suffix=flux_path, timepoint="end", species_path=species_path)
    return


#######################################################################
# Input Parameters for combustor
#######################################################################

# filepath for writing files
if len(sys.argv) == 3:
    # Pass the cti file and rmg file as an argument to the script
    cti_file = sys.argv[1]
    rmg_file = sys.argv[2]
else:
    # cti_file = os.path.dirname(os.path.abspath(__file__)) +'/chem_annotated.cti'
    cti_file = "/work/westgroup/ChrisB/meoh-synthesis_RMG/meOH-synthesis/base/cantera/chem_annotated.cti"
    # cti_file = "/work/westgroup/ChrisB/meOH-synthesis_cantera/meOH-synthesis/External_data/mech_grabow_new.cti"
    rmg_file = "/work/westgroup/ChrisB/meoh-synthesis_RMG/meOH-synthesis"


# Reactor settings arrays for run
Temps = [400, 500, 528, 600]
Pressures = [15, 30, 50, 75]
volume_flows = [0.00424, 0.0106, 0.02544]

# Mole fractions from runs for reference (min max and mid from Graaf data)
# X_cos = [0.053, 0.1365, 0.220]
# X_co2s = [0.261,0.1205, 0.261]
# X_h2s = [0.625, 0.7625, 0.900]

# CO+CO2/H2
# H2_fraction = [0.8,0.5,0.95,0.75]
H2_fraction = [0.8, 0.5, 0.95, 0.75]

# CO2/CO
CO_CO2_ratio = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

# reaction time
reactime = 1e4

# sensitivity settings
sensitivity = False
sensatol = 1e-6
sensrtol = 1e-6

grabow = False

run_reactor(
    cti_file=cti_file,
    rmg_model_path=rmg_file,
    t_array=Temps,
    reactor_type=1,
    h2_array=H2_fraction,
    co2_array=CO_CO2_ratio,
    energy="off",
    sensitivity=sensitivity,
    sensatol=sensatol,
    sensrtol=sensrtol,
    reactime=reactime,
    grabow=grabow
)
