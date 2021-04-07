import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cycler import cycler

save_name_hist = "WT-protein-system-pucker-histogram.png"

system = "WT Protein System"

## Read in data
## header = 0 reads header in first row because Python starts at 0

data = {}
i = 1
## Pulls residues 431-450 (select the nucleic acids)
for x in range(431,450+1):
        d_a = pd.read_csv(f'r1/WT-protein-system-pucker-{x}.dat', \
         delim_whitespace=True, header=0)
        d_b = pd.read_csv(f'r2/WT-protein-system-pucker-{x}.dat', \
         delim_whitespace=True, header=0)
        d_c = pd.read_csv(f'r3/WT-protein-system-pucker-{x}.dat', \
         delim_whitespace=True, header=0)

        data[i] = pd.concat([d_a, d_b, d_c])
        ## Do funky stuff with counter because of how DNA residues are labelled
        if i < 10:
            i += 1
        elif i == 10:
            i = 20
        else:
            i -= 1

#--------- Past the curtain, change the residue labels in a1-a20! ---------#

## Define the column names to use throughout
## This renames the columns in case something was weird

for x in data:
    data[x].columns = ['Frame', 'Pucker']

plt.rcParams.update({'axes.prop_cycle': cycler('color', ['#88CCEE', '#44AA99',
'#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499', '#332288'])})

## Set up shape of figure (5 rows, 2 columns)
fig = plt.figure(figsize=(15,12))
fig.tight_layout(rect=[0, 0.03, 1, 0.95])
a1 = fig.add_subplot(541) #1
a11 = fig.add_subplot(542) #2
a6 = fig.add_subplot(543) #3
a16 = fig.add_subplot(544) #4
# New row
a2 = fig.add_subplot(545) #5
a12 = fig.add_subplot(546) #6
a7 = fig.add_subplot(547) #7
a17 = fig.add_subplot(548) #8
# New row
a3 = fig.add_subplot(549) #9
a13 = fig.add_subplot(5,4,10)
a8 = fig.add_subplot(5,4,11)
a18 = fig.add_subplot(5,4,12)
# New row
a4 = fig.add_subplot(5,4,13)
a14 = fig.add_subplot(5,4,14)
a9 = fig.add_subplot(5,4,15)
a19 = fig.add_subplot(5,4,16)
# New row
a5 = fig.add_subplot(5,4,17)
a15 = fig.add_subplot(5,4,18)
a10 = fig.add_subplot(5,4,19)
a20 = fig.add_subplot(5,4,20)

# Set all x- and y- axes
# Use 0-360 for angles
c_xlim = (0,360)
# Set probability density max. Might need to play with this
c_ylim = (0.00, 0.025)

## Add data to each subplot and add label
## Label needs the space because you're adding strings
a1.hist(data[1]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="A1")
a1.set_xlim(c_xlim)
a1.set_ylim(c_ylim)
a1.axvline(data[1]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[1]["Pucker"].mean()))
handles, labels = a1.get_legend_handles_labels()
a1.legend(reversed(handles), reversed(labels))

a2.hist(data[2]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="C2")
a2.set_xlim(c_xlim)
a2.set_ylim(c_ylim)
a2.axvline(data[2]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[2]["Pucker"].mean()))
handles, labels = a2.get_legend_handles_labels()
a2.legend(reversed(handles), reversed(labels))

a3.hist(data[3]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="C3")
a3.set_xlim(c_xlim)
a3.set_ylim(c_ylim)
a3.axvline(data[3]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[3]["Pucker"].mean()))
handles, labels = a3.get_legend_handles_labels()
a3.legend(reversed(handles), reversed(labels))

a4.hist(data[4]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="A4")
a4.set_xlim(c_xlim)
a4.set_ylim(c_ylim)
a4.axvline(data[4]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[4]["Pucker"].mean()))
handles, labels = a4.get_legend_handles_labels()
a4.legend(reversed(handles), reversed(labels))

a5.hist(data[5]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="C5")
a5.set_xlim(c_xlim)
a5.set_ylim(c_ylim)
a5.axvline(data[5]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[5]["Pucker"].mean()))
handles, labels = a5.get_legend_handles_labels()
a5.legend(reversed(handles), reversed(labels))

a6.hist(data[6]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="C6")
a6.set_xlim(c_xlim)
a6.set_ylim(c_ylim)
a6.axvline(data[6]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[6]["Pucker"].mean()))
handles, labels = a6.get_legend_handles_labels()
a6.legend(reversed(handles), reversed(labels))

a7.hist(data[7]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="G7")
a7.set_xlim(c_xlim)
a7.set_ylim(c_ylim)
a7.axvline(data[7]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[7]["Pucker"].mean()))
handles, labels = a7.get_legend_handles_labels()
a7.legend(reversed(handles), reversed(labels))

a8.hist(data[8]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="G8")
a8.set_xlim(c_xlim)
a8.set_ylim(c_ylim)
a8.axvline(data[8]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[8]["Pucker"].mean()))
handles, labels = a8.get_legend_handles_labels()
a8.legend(reversed(handles), reversed(labels))

a9.hist(data[9]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="T9")
a9.set_xlim(c_xlim)
a9.set_ylim(c_ylim)
a9.axvline(data[9]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[9]["Pucker"].mean()))
handles, labels = a9.get_legend_handles_labels()
a9.legend(reversed(handles), reversed(labels))

a10.hist(data[10]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="G10")
a10.set_xlim(c_xlim)
a10.set_ylim(c_ylim)
a10.axvline(data[10]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[10]["Pucker"].mean()))
handles, labels = a10.get_legend_handles_labels()
a10.legend(reversed(handles), reversed(labels))

a11.hist(data[11]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="T1'", color='#CC6677')
a11.set_xlim(c_xlim)
a11.set_ylim(c_ylim)
a11.axvline(data[11]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[11]["Pucker"].mean()))
handles, labels = a11.get_legend_handles_labels()
a11.legend(reversed(handles), reversed(labels))

a12.hist(data[12]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="G2'", color='#CC6677')
a12.set_xlim(c_xlim)
a12.set_ylim(c_ylim)
a12.axvline(data[12]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[12]["Pucker"].mean()))
handles, labels = a12.get_legend_handles_labels()
a12.legend(reversed(handles), reversed(labels))

a13.hist(data[13]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="G3'", color='#CC6677')
a13.set_xlim(c_xlim)
a13.set_ylim(c_ylim)
a13.axvline(data[13]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[13]["Pucker"].mean()))
handles, labels = a13.get_legend_handles_labels()
a13.legend(reversed(handles), reversed(labels))

a14.hist(data[14]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="T4'", color='#CC6677')
a14.set_xlim(c_xlim)
a14.set_ylim(c_ylim)
a14.axvline(data[14]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[14]["Pucker"].mean()))
handles, labels = a14.get_legend_handles_labels()
a14.legend(reversed(handles), reversed(labels))

a15.hist(data[15]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="G5'", color='#CC6677')
a15.set_xlim(c_xlim)
a15.set_ylim(c_ylim)
a15.axvline(data[15]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[15]["Pucker"].mean()))
handles, labels = a15.get_legend_handles_labels()
a15.legend(reversed(handles), reversed(labels))

a16.hist(data[16]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="G6'", color='#CC6677')
a16.set_xlim(c_xlim)
a16.set_ylim(c_ylim)
a16.axvline(data[16]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[16]["Pucker"].mean()))
handles, labels = a16.get_legend_handles_labels()
a16.legend(reversed(handles), reversed(labels))

a17.hist(data[17]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="C7'", color='#CC6677')
a17.set_xlim(c_xlim)
a17.set_ylim(c_ylim)
a17.axvline(data[17]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[17]["Pucker"].mean()))
handles, labels = a17.get_legend_handles_labels()
a17.legend(reversed(handles), reversed(labels))

a18.hist(data[18]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="C8'", color='#CC6677')
a18.set_xlim(c_xlim)
a18.set_ylim(c_ylim)
a18.axvline(data[18]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[18]["Pucker"].mean()))
handles, labels = a18.get_legend_handles_labels()
a18.legend(reversed(handles), reversed(labels))

a19.hist(data[19]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="A9'", color='#CC6677')
a19.set_xlim(c_xlim)
a19.set_ylim(c_ylim)
a19.axvline(data[19]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[19]["Pucker"].mean()))
handles, labels = a19.get_legend_handles_labels()
a19.legend(reversed(handles), reversed(labels))

a20.hist(data[20]["Pucker"], bins=99, edgecolor='black', linewidth=1.2, \
 density=True,  alpha=0.7, label="C10'", color='#CC6677')
a20.set_xlim(c_xlim)
a20.set_ylim(c_ylim)
a20.axvline(data[20]["Pucker"].mean(), color='k', linestyle='-.', linewidth=1, \
 label='{:.2f}'.format(data[20]["Pucker"].mean()))
handles, labels = a20.get_legend_handles_labels()
a20.legend(reversed(handles), reversed(labels))

## Set up shared axis labels
fig.text(0.5, 0.04, 'Pucker ($^\circ$)', ha='center', size="20")
fig.text(0.04, 0.5, 'Probability Density', va='center', rotation='vertical', \
 size="20")
fig.suptitle(system, fontsize=22)

plt.savefig(save_name_hist, dpi=300)
plt.close(save_name_hist)
