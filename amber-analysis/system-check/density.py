#! /usr/bin/python3
## This script works with data generated from process_mdout.perl (from AMBER)

import matplotlib.pyplot as plt
import pandas as pd

save_name_dens = 'test-density.png'

## Labels on plot
sys_lab = "WT"

## Color for the plots -- use Matplotlib color names or a hex code
color = '#EE3377'

## Read in Density Data
## header = 0 reads header in first row because Python starts at 0
d1_dens = pd.read_csv('summary.DENSITY', \
 delim_whitespace=True, header=None)

## Fraction to multiply to convert picoseconds to nanoseconds
ps2ns = 0.001

#---------- No need to modify past the curtain ---------------#

## This renames the columns in case something was weird
d1_dens.columns = ['Time', 'Density']

plt.plot(d1_dens['Time']*ps2ns, d1_dens['Density'], label=sys_lab, color=color)

plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2)

plt.xlabel('Time (ns)', labelpad=5)
plt.ylabel('Density (g cm$^-3$)', labelpad=10)

#plt.savefig(save_name, dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.savefig(save_name_dens, dpi=300)
plt.close(save_name_dens)
plt.clf()
