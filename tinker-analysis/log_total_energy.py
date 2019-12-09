"""Get the total energy from a Tinker log file (Potential + Kinetic) by frame"""

import numpy as np
import pandas as pd
import re

tinkerlog="my_tinker_log_file.out"
outfile="total_energy.dat"

#energy_lines='pot_kin.dat'

# ## grep "Current Potential" tinkerlog
def read_log(tinkerlog):
    """Initially search (based on grep) for potential and kinetic energy
    lines from the Tinker log file. Then, create a Pandas data frame with the
    potential and kinetic energy values."""
    potential_lines = []
    pattern = '(?i)Potential' # case insensitive
    for line in open(tinkerlog).readlines():
        # re.search is anywhere, re.match is explicitly beginning of line
        if re.search(pattern, line):
            potential_lines.append(line)
    kinetic_lines = []
    pattern = '(?i)Kinetic' # case insensitive
    for line in open(tinkerlog).readlines():
        if re.search(pattern, line):
            kinetic_lines.append(line)
    df = pd.DataFrame(list(zip(potential_lines,kinetic_lines)),
     columns=['P', 'K'])
    ## String split the relevant lines
    df[['P Type', 'PE']] = df["P"].str.split("      ", 1, expand=True)
    df["PE"] = df["PE"].str.strip() # Remove leading & trailing spaces + newlines
    df[['PotentialE', 'P Unit']] = df["PE"].str.split(" ", 1, expand=True)
    df[['K Type', 'KE']] = df["K"].str.split("      ", 1, expand=True)
    df["KE"] = df["KE"].str.strip() # Remove leading & trailing spaces + newlines
    df[['KineticE', 'K Unit']] = df["KE"].str.split(" ", 1, expand=True)
    ## Remove the irrelevant columns
    df = df.drop(columns=['P', 'P Type', 'PE', 'K', 'K Type', 'KE'])
    # ## Write a file of all the energy lines for later use
    # df.to_csv(test_energies_csv, index=False, encoding='utf8')
    return df

def calculate_TE(df, outfile):
    """Add the potential and kinetic energies, making a new total energy column.
    Then write an outfile with total energy by frame."""
    df['PotentialE'] = df['PotentialE'].astype(float)
    df['KineticE'] = df['KineticE'].astype(float)
    df['TotalE'] = df['PotentialE'] + df['KineticE']
    energies = df['TotalE']
    ## Set index to 1 for writing
    energies.index = energies.index + 1
    energies.index.name = "#Frame"
    energies.to_csv(outfile, sep='\t', index=True, encoding='utf8', header=True)
    return energies


df = read_log(tinkerlog)
energies = calculate_TE(df, outfile)
# print(df)
