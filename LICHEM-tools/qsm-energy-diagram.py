import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

## Name of the LICHEM log file
log_file = "LICHEM_QSM.log"

## Name of the output files for the initial guess and final optimization
in_fig="initial_QSM_RE.png"
fin_fig="final_QSM_RE.png"

#---------- No need to modify past the curtain ---------------#
plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'figure.autolayout': True})

def read_log(log_file):
    '''Read the LICHEM log file
    Parameters
    ----------
    log_file : str
        The file name of the LICHEM log file.
    Returns
    -------
    init_RE_kcmol : list
        The reaction energy values from the inital optimization in the log file.
    fin_RE_kcmol : list
        The reaction energy values from the final optimization in the log file.
    '''
    f = open(log_file, 'r')
    log_lines = f.readlines()
    #
    save_lines_on = False
    save_lines = []
    #
    group = 0
    #
    for line in log_lines:
        ## Turns on the line saving for Bead Energies
        # Bead energies:
        if "Bead " in line.strip().split("|")[0]:
            save_lines_on = True
        ## Check for line after beads
        # | Backward Barrier :
        elif "| TS Bead          " in line.strip().split(":")[0]:
            save_lines_on = False
            group += 1
        ## Save the lines with bead information
        elif save_lines_on:
            save_lines.append(line.strip('\n').strip() + "| {}".format(group))
    #
    f.close()
    #
    ## Create a DataFrame of the barriers at each step
    df = pd.DataFrame([x.split("|") for x in save_lines],
     columns = ["Bead", "Coord", "E (a.u.)", "RE (kcal/mol)", "Group"],
     dtype=float)
    #
    ## Drop the line breaks between groups
    df = df.dropna(axis="index")
    #
    ## Get inital guess and final optimized values
    initial = df[df.Group == df.Group.min()]
    final = df[df.Group == df.Group.max()]
    #
    ## Convert to energy lists
    coord = list(initial['Coord'])
    init_RE_kcmol = list(initial['RE (kcal/mol)'])
    fin_RE_kcmol = list(final['RE (kcal/mol)'])
    #
    return coord, init_RE_kcmol, fin_RE_kcmol

def rxn_energy_plot(outputfile, energies, coord):
    """Create a reaction coordinate diagram.
    Parameters
    ----------
    outputfile : str
        The file name to use for the generated plot.
    energies : list
        Energy values for the beads along the reaction coordinate.
    coord : list
        Location of the beads along the reaction coordinate.
    Returns
    -------
    outputfile : png
        The generated plot.
    """
    ## Set up the plot
    fig=plt.figure(figsize=(12,8),dpi=300)
    ax=fig.add_subplot(111)
    ## Plot the dashed line
    plt.plot(coord, energies, "k--")
    ## Add the plot values
    for i in range(len(energies)):
        plt.annotate(str(energies[i]),xy=(coord[i],energies[i]), ha='center',
        va='center', bbox=dict(boxstyle="round", fc="0.8"),
        annotation_clip=False)
    ## Add line at zero
    plt.axhline(y=0, c='0.55', linestyle='dotted')
    ## Add axis labels
    plt.xlabel('Reaction Coordinate', labelpad=20)
    plt.ylabel('Reaction Energy (kcal/mol)', labelpad=20)
    ## Save the figure
    plt.savefig(outputfile,dpi=300)

## Run the code!!!!
coord, init_RE_kcmol, fin_RE_kcmol = read_log(log_file)
rxn_energy_plot(in_fig, init_RE_kcmol, coord)
rxn_energy_plot(fin_fig, fin_RE_kcmol, coord)
