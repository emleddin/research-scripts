# Get the total energy from a GROMACS log file by frame and plot an energy term

import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
from cycler import cycler

## Type of log file ("minimization" or "md") [printout values differ]
type="md"
## Input log file
gromacs_log="my_gromacs_log_file.log"

## Output data file
outfile="separated_energy_gromacs.dat"
## Save using "kj" or "kcal"
unit = "kJ"

## Energy plot type (choose from log file values)
plot_type = "Total Energy"
## Output plot
save_plot = "total_energy_gromacs.png"

## Convert frames (picoseconds) to nanoseconds
frame2ns = 0.001
## Convert kilojoules/mol to kcal/mol
kj2kcal=1/4.184

## Change plot colors
plt.rcParams.update({'axes.prop_cycle': cycler('color', ['#88CCEE', '#44AA99',
'#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499', '#332288'])})

#--------------------------- Function Definitions ---------------------------#

def read_min_log(gromacs_log, unit):
    """Initially search for energy terms from the minimization GROMACS log file.
    Then, create a Pandas data frame with the different energy terms.
    Bond            Angle       Proper Dih.     Improper Dih.    LJ-14
    Coulomb-14      LJ (SR)     Disper. corr.   Coulomb (SR)     Coul. recip.
    Position Rest.  Potential   Pres. DC (bar)  Pressure (bar)   Constr. rmsd
    """
    ## Make empty lists for each of the energy types
    bond_lines = []
    angle_lines = []
    prop_dihed_lines = []
    imp_dihed_lines = []
    lj_14_lines = []
    coul_14_lines = []
    lj_sr_lines = []
    disper_corr_lines = []
    coul_sr_lines = []
    coul_recip_lines = []
    pos_rest_lines = []
    potential_lines = []
    pres_dc_bar_lines = []
    pres_bar_lines = []
    constr_rmsd_lines = []
    pattern = '(?i)Energies \(kJ/mol\)' # case insensitive, escape parentheses

    s_line = 0
    for line in open(gromacs_log).readlines():
        # re.search is anywhere, re.match is explicitly beginning of line
        if re.search(pattern, line):
            # If it matches "Energies", start counter
            s_line += 1
        elif s_line == 1:
            s_line += 1
        elif s_line == 2:
            sl = line.split()
            bond_lines.append(sl[0])
            angle_lines.append(sl[1])
            prop_dihed_lines.append(sl[2])
            imp_dihed_lines.append(sl[3])
            lj_14_lines.append(sl[4])
            s_line += 1
        elif s_line == 3:
            s_line += 1
        elif s_line == 4:
            sl = line.split()
            coul_14_lines.append(sl[0])
            lj_sr_lines.append(sl[1])
            disper_corr_lines.append(sl[2])
            coul_sr_lines.append(sl[3])
            coul_recip_lines.append(sl[4])
            s_line += 1
        elif s_line == 5:
            s_line += 1
        elif s_line == 6:
            sl = line.split()
            pos_rest_lines.append(sl[0])
            potential_lines.append(sl[1])
            pres_dc_bar_lines.append(sl[2])
            pres_bar_lines.append(sl[3])
            constr_rmsd_lines.append(sl[4])
            s_line = 0

    df = pd.DataFrame(list(zip(
    bond_lines, angle_lines, prop_dihed_lines, imp_dihed_lines, lj_14_lines,
    coul_14_lines, lj_sr_lines, disper_corr_lines, coul_sr_lines,
     coul_recip_lines,
    pos_rest_lines, potential_lines, pres_dc_bar_lines, pres_bar_lines,
     constr_rmsd_lines)),
     columns=["Bond", "Angle", "Proper Dih.", "Improper Dih.", "LJ-14",
      "Coulomb-14", "LJ (SR)", "Disper. corr.", "Coulomb (SR)", "Coul. recip.",
      "Position Rest.", "Potential", "Pres. DC (bar)", "Pressure (bar)",
       "Constr. rmsd"])
    ## Save data file
    df.to_csv(outfile, sep='\t', index=True, encoding='utf8', header=True)
    return df

def read_md_log(gromacs_log, unit):
    """Initially search for energy terms from the MD GROMACS log file.
    Then, create a Pandas data frame with the different energy terms.
    Bond            Angle           Proper Dih.     Improper Dih.  LJ-14
    Coulomb-14      LJ (SR)         Disper. corr.   Coulomb (SR)   Coul. recip.
    Position Rest.  Potential       Kinetic En.     Total Energy   Conserved En.
    Temperature     Pres. DC (bar)  Pressure (bar)  Constr. rmsd
    """
    ## Make empty lists for each of the energy types
    bond_lines = []
    angle_lines = []
    prop_dihed_lines = []
    imp_dihed_lines = []
    lj_14_lines = []
    coul_14_lines = []
    lj_sr_lines = []
    disper_corr_lines = []
    coul_sr_lines = []
    coul_recip_lines = []
    pos_rest_lines = []
    potential_lines = []
    kinetic_lines = []
    tot_en_lines = []
    conserv_en_lines = []
    temp_lines = []
    pres_dc_bar_lines = []
    pres_bar_lines = []
    constr_rmsd_lines = []
    pattern = '(?i)Energies \(kJ/mol\)' # case insensitive, escape parentheses

    s_line = 0
    for line in open(gromacs_log).readlines():
        # re.search is anywhere, re.match is explicitly beginning of line
        if re.search(pattern, line):
            # If it matches "Energies", start counter
            s_line += 1
        elif s_line == 1:
            s_line += 1
        elif s_line == 2:
            sl = line.split()
            bond_lines.append(sl[0])
            angle_lines.append(sl[1])
            prop_dihed_lines.append(sl[2])
            imp_dihed_lines.append(sl[3])
            lj_14_lines.append(sl[4])
            s_line += 1
        elif s_line == 3:
            s_line += 1
        elif s_line == 4:
            sl = line.split()
            coul_14_lines.append(sl[0])
            lj_sr_lines.append(sl[1])
            disper_corr_lines.append(sl[2])
            coul_sr_lines.append(sl[3])
            coul_recip_lines.append(sl[4])
            s_line += 1
        elif s_line == 5:
            s_line += 1
        elif s_line == 6:
            sl = line.split()
            pos_rest_lines.append(sl[0])
            potential_lines.append(sl[1])
            kinetic_lines.append(sl[2])
            tot_en_lines.append(sl[3])
            conserv_en_lines.append(sl[4])
            s_line += 1
        elif s_line == 7:
            s_line += 1
        elif s_line == 8:
            temp_lines.append(sl[0])
            pres_dc_bar_lines.append(sl[1])
            pres_bar_lines.append(sl[2])
            constr_rmsd_lines.append(sl[3])
            s_line = 0

    df = pd.DataFrame(list(zip(
    bond_lines, angle_lines, prop_dihed_lines, imp_dihed_lines, lj_14_lines,
    coul_14_lines, lj_sr_lines, disper_corr_lines, coul_sr_lines,
     coul_recip_lines,
    pos_rest_lines, potential_lines, kinetic_lines, tot_en_lines,
     conserv_en_lines,
    temp_lines, pres_dc_bar_lines, pres_bar_lines, constr_rmsd_lines)),
     columns=["Bond", "Angle", "Proper Dih.", "Improper Dih.","LJ-14",
     "Coulomb-14", "LJ (SR)", "Disper. corr.", "Coulomb (SR)", "Coul. recip.",
     "Position Rest.", "Potential", "Kinetic En.", "Total Energy",
      "Conserved En.",
     "Temperature", "Pres. DC (bar)", "Pressure (bar)", "Constr. rmsd"])
    # Change values to float
    df = df.astype(np.float64)
    ## Save averages elsewhere
    df_avg = df.tail(1).index
    ## Convert if necessary
    if unit == "kcal":
        df *= kj2kcal
    ## Drop the last "AVERAGES" row
    df = df.drop(df.tail(1).index)
    ## Save data file
    df.to_csv(outfile, sep='\t', index=True, encoding='utf8', header=True)
    return df, df_avg

def plot_energy(df, plot_type, save_plot, unit):
    plt.plot(df.index.values*frame2ns, df["{}".format(plot_type)])
    plt.xlabel('Time (ns)')
    ## Plot using correct units
    if unit.lower() == "kj":
        plt.ylabel('Energy (kJ/mol)')
    elif unit.lower() == "kcal":
        plt.ylabel('Energy (kcal/mol)')
    plt.tight_layout()
    plt.savefig(save_plot, dpi=300)#, bbox_inches='tight', pad_inches=0)

#------------------------------- Run the Code -------------------------------#
## Get the info
if type.lower() == "minimization":
    df = read_min_log(gromacs_log, unit)
elif type.lower() == "md":
    df, df_avg = read_md_log(gromacs_log, unit)

## Make a plot
plot_energy(df, plot_type, save_plot, unit)
