#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# ----- Define variables ---- #

# Total energies
energies = [-52.3, -38.5, -40.1, -28.6, -43.4, -39.2, -58.4]
# energies = [-1.439, -1.445, -1.423]

# Units for the given energies of each component -- (au, kcal/mol, eV)
e_units = "kcal/mol"
# e_units = "au"

# Output Units -- (au, kcal/mol, eV)
rel_units = "kcal/mol"

# Optional: Specify labels to correspond to each energy
#           Default: Reactant, Intermediate / [I1, I2, ...], Product
# labels = ["Reactant", "TS", "Product"]

# Name of output file
out_fig = "relative-free-energy.png"

# Output dimensions
out_width = 5
out_height = 3

# --------------------------- Behind the Curtain --------------------------- #

# Future idea: pull total energy directly from log files if given paths

# Constants as implemented in LICHEM
avogNum = 6.02214129e23
har2eV = 27.21138505
SI2eV = 1/(1.602176565e-19)

kcal2eV = 4184*SI2eV/avogNum
au2kcal = 627.51

# Change font size
# plt.rcParams.update({'font.size': 12})


def get_relative_energy(energies, e_units, rel_units):
    """
    Convert the energy values to a relative scale starting at 0.

    Parameters
    ----------
    energies : list
        Energy values for the structures.
    e_units : str
        The units of each member of `energies`.
    rel_units : str
        The units to use for the relative energies.

    Returns
    -------
    rel_energies : list
        Energy values converted to a relative scale.
    rel_units : str
        The units of each member the `rel_energies`.
    """
    global au2kcal
    global har2eV
    global kcal2eV
    # List of relative energies
    rel_energies = []
    e_rxt = energies[0]
    for e in energies:
        rel_energies.append(e - e_rxt)
    # Check that given units make sense
    if e_units.lower() not in ("au", "kcal", "kcal/mol", "ev"):
        print(f"\nWARNING: Provided units, {e_units}, not understood.\n"
              f"         Units will not be converted to {rel_units}.")
        rel_units = e_units
    if rel_units.lower() not in ("au", "kcal", "kcal/mol", "ev"):
        print(f"\nWARNING: Requested output units, {rel_units}, "
              "not understood.\n"
              f"         Units will not be converted from {e_units}.")
        rel_units = e_units
    # Convert units
    if e_units.lower() == "au" and rel_units.lower() in ("kcal", "kcal/mol"):
        rel_energies = [rel * au2kcal for rel in rel_energies]
    elif e_units.lower() in ("kcal", "kcal/mol") and rel_units.lower() == "au":
        rel_energies = [rel / au2kcal for rel in rel_energies]
    elif e_units.lower() == "au" and rel_units.lower() == "ev":
        rel_energies = [rel * har2eV for rel in rel_energies]
    elif e_units.lower() == "ev" and rel_units.lower() == "au":
        rel_energies = [rel / har2eV for rel in rel_energies]
    elif e_units.lower() in ("kcal", "kcal/mol") and rel_units.lower() == "ev":
        rel_energies = [rel * kcal2eV for rel in rel_energies]
    elif e_units.lower() == "ev" and rel_units.lower() in ("kcal", "kcal/mol"):
        rel_energies = [rel / kcal2eV for rel in rel_energies]
    return rel_energies, rel_units


def plot_energies(rel_energies, rel_units, out_width, out_height, labels=None):
    """
    Plot and save a relative energy diagram.

    Parameters
    ----------
    rel_energies : list
        Energy values converted to a relative scale.
    rel_units : str
        The units of each member the `rel_energies`.
    out_width : real
        Width of the output figure in inches.
    out_height : real
        Height of the output figure in inches.
    labels : list
        Labels for each respective structure.
        If None, labels are: Reactant, Intermediate / [I1, I2, ...], Product

    Notes
    -----
    Based off Free_Energy_Diagram.py by GitHub user markahix
    https://github.com/markahix/Basic-Scripts/blob/main/Basic_Diagrams/Free_Energy_Diagram.py
    """
    # Figure out current font size
    font_size = plt.rcParams['font.size']
    # print(f"Font size: {font_size}")
    # Figure out how many structures there are
    structs = len(rel_energies)
    # If more than 3 structures, assume multiple intermediates
    if structs > 3:
        mult_int = True
    else:
        mult_int = False
    # Scale width of each structure based on total structures
    sx1 = (out_width) / (structs)
    # Scale vertical text placement based on font size
    sy1 = (font_size) / (out_height)
    sx_75 = sx1 - 0.25
    sx1_25 = sx1 + 0.25
    # Scale output linewidth
    fig_lw = 0.5 * out_height
    # Set outline size for border around reactant/product labels
    structname_ol = 0.25
    # Set color for border around reactant/product labels
    # col = "k"
    col = "none"
    # Check the number of labels given (if any)
    if labels is not None:
        if len(labels) != structs:
            print("ERROR: The number of labels provided does not match the "
                  "number of energies.\n"
                  f"       Expected {structs}, but got "
                  f"{len(labels)} ({labels}).\n"
                  "       Continuing with default labels.")
            labels = None
    # List for all bounding boxes
    bboxes = []
    # Set up figure
    fig = plt.figure(figsize=(out_width, out_height), dpi=300)
    ax = fig.add_subplot(111)
    # Since python plots sandwich-style, start with gray line at zero energy
    plt.axhline(0, color='gray', ls='dotted', lw=0.75*fig_lw, alpha=0.75)
    for i in range(0, structs):
        # A horizontal line of set width for each energy value
        # [x1, x2], [y1, y2], style (k is black)
        plt.plot([i + sx_75,       i + sx1_25],
                 [rel_energies[i], rel_energies[i]],
                 "k", lw=fig_lw)
        if i != structs-1:
            # Slope between each vertical line
            plt.plot([i + sx1_25,      i + 1 + sx_75],
                     [rel_energies[i], rel_energies[i+1]],
                     "k--", lw=0.25*fig_lw)
        # Include the energy centered below the horizontal line
        rel_label = str(round(rel_energies[i], 3)) + " " + rel_units
        ax.annotate(str(round(rel_energies[i], 3)),
                    xy=(i + sx1, rel_energies[i] - sy1), ha='center',
                    # Put the energy in a white box with background
                    #  transparency so you can read it regardless of slope
                    bbox=dict(facecolor='#FFFFFF75', edgecolor='none',
                              boxstyle='round,pad=.05'))
        # Label the structure
        if labels is None:
            # Reactant
            if i == 0:
                b = ax.annotate(
                        "Reactant",
                        xy=(i+sx1, rel_energies[i]-2*sy1),
                        ha='center', weight="bold",
                        bbox=dict(facecolor='#FFFFFF75', edgecolor=col,
                                  boxstyle='round,pad=.1', lw=structname_ol))
                bboxes.append(b)
            # Product
            elif i == structs-1:
                b = ax.annotate(
                        "Product",
                        xy=(i+sx1, rel_energies[i]-2*sy1),
                        ha='center', weight="bold",
                        bbox=dict(facecolor='#FFFFFF75', edgecolor=col,
                                  boxstyle='round,pad=.1', lw=structname_ol))
                bboxes.append(b)
            # Multiple intermediates
            elif mult_int is True:
                b = ax.annotate(
                        "I" + str(i+1),
                        xy=(i+sx1, rel_energies[i]-2*sy1),
                        ha='center', weight="bold",
                        bbox=dict(facecolor='#FFFFFF75', edgecolor=col,
                                  boxstyle='round,pad=.1', lw=structname_ol))
                bboxes.append(b)
            else:
                b = ax.annotate(
                        "Intermediate",
                        xy=(i+sx1, rel_energies[i]-2*sy1),
                        ha='center', weight="bold",
                        bbox=dict(facecolor='#FFFFFF75', edgecolor=col,
                                  boxstyle='round,pad=.1', lw=structname_ol))
                bboxes.append(b)
        else:
            b = ax.annotate(
                    str(labels[i]),
                    xy=(i+sx1, rel_energies[i]-2*sy1),
                    ha='center', weight="bold",
                    bbox=dict(facecolor='#FFFFFF75', edgecolor=col,
                              boxstyle='round,pad=.1', lw=structname_ol))
            bboxes.append(b)
    # ----------------- #
    # Finishing Touches #
    # ----------------- #
    # Get min/max x coordinates based on label text: bboxes[i].xy[0]
    bx_min = min(bboxes, key=lambda item: item.xy[0]).xy[0]
    bx_max = max(bboxes, key=lambda item: item.xy[0]).xy[0]
    # Get minimum y coordinates based on label text: bboxes[i].xy[1]
    by_min = min(bboxes, key=lambda item: item.xy[1]).xy[1]
    # Fix the plot limits to show all the annotations
    # X is leftmost annotation to rightmost annotation +/- structs/out_width
    plt.xlim([bx_min-1/sx1, bx_max+1/sx1])
    # Y is lowest annotation to max energy +/- base y_width
    plt.ylim([by_min-sy1, max(rel_energies)+sy1])
    # Turn of x-axis labels
    plt.tick_params(
        axis='x',           # changes apply to the x-axis
        which='both',       # both major and minor ticks are affected
        bottom=False,       # ticks along the bottom edge are off
        top=False,          # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off
    # Set y-axis label
    plt.ylabel("Relative Energy (" + rel_units + ")")
    # Save the plot
    plt.tight_layout()
    plt.savefig(out_fig, dpi=300)
    #
    return


# ---------------------------- Run the Script ---------------------------- #

rel_energies, rel_units = get_relative_energy(energies, e_units, rel_units)

try:
    plot_energies(rel_energies, rel_units, out_width, out_height, labels)
except NameError:
    plot_energies(rel_energies, rel_units, out_width, out_height)
