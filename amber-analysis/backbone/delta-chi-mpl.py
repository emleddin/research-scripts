import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cycler import cycler

save_name = "WT-MUT-dihedrals-X.png"

## Read in data
## header = 0 reads header in first row because Python starts at 0
## -------- System 1
d1_a = pd.read_csv('rep1/MUT-protein-system-X-dihedral.dat', \
 delim_whitespace=True, header=0)

d1_b = pd.read_csv('rep2/MUT-protein-system-X-dihedral.dat', \
 delim_whitespace=True, header=0)

d1_c = pd.read_csv('rep3/MUT-protein-system-X-dihedral.dat', \
 delim_whitespace=True, header=0)

## Concatenate together into 1 dataframe
d1 = pd.concat([d1_a, d1_b, d1_c])

sys_a = "MUT"

## -------- System 2
d2_a = pd.read_csv('rep1/WT-protein-system-X-dihedral.dat', \
 delim_whitespace=True, header=0)

d2_b = pd.read_csv('rep2/WT-protein-system-X-dihedral.dat', \
 delim_whitespace=True, header=0)

d2_c = pd.read_csv('rep3/WT-protein-system-X-dihedral.dat', \
 delim_whitespace=True, header=0)

d2 = pd.concat([d2_a, d2_b, d2_c])

sys_b = "WT"

residue_name = "X"

#---------- No need to modify past the curtain ---------------#

## Define the column names to use throughout
## This renames the columns in case something was weird
d1.columns = ['Frame', 'Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon', 'Zeta', 'Chi']
d2.columns = ['Frame', 'Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon', 'Zeta', 'Chi']

## Change font size
# plt.rcParams.update({'font.size': 24})
plt.rcParams.update({'figure.autolayout': True})

## Set Up Graphs
## Set up different color axes based on https://personal.sron.nl/~pault/
# plt.rcParams.update({'axes.prop_cycle': cycler('color', ['#CC6677', '#88CCEE',
# '#44AA99', '#117733', '#999933', '#DDCC77', '#882255', '#AA4499', '#332288'])})

plt.xlim(0,360)
plt.ylim(0,360)

## Plot the First Data
plt.scatter(d1['Chi'], d1['Delta'], alpha=0.3, marker='.', \
 color='#CC6677', label=(sys_a + " " + residue_name))

plt.scatter(d2['Chi'], d2['Delta'], alpha=0.3, marker='.', \
 color='#88CCEE', label=(sys_b + " " + residue_name))

## Flip the legend so WT is listed first...
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(reversed(handles), reversed(labels))

## Add right and top ticks
plt.tick_params(right=True, top=True, direction='in')

## Set the axis labels
plt.xlabel('$\chi$ Degrees ($^\circ$)', labelpad=5)
plt.ylabel('$\delta$ Degrees ($^\circ$)', labelpad=10)

plt.savefig(save_name, dpi=300)
plt.close(save_name)
