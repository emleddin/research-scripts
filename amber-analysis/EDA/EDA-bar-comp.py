#! /usr/bin/python3
## This script works with data generated from running `rmagic-EDA-avg.r` or
## `rmagic-EDA-single-run.r` on 2 different systems

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

save_name = 'sysA_sysB_total_barchart.png'
sys_a_tag = 'Sys A'
sys_b_tag = 'Sys B'

## Read in data
## header = 0 reads header in first row because Python starts at 0
d1 = pd.read_csv('SysA_protein_system_res222_tot_avg.dat', \
 delim_whitespace=True, header=0)
d2 = pd.read_csv('SysB_protein_system_res222_tot_avg.dat', \
 delim_whitespace=True, header=0)

#---------- No need to modify past the curtain ---------------#

## This renames the columns in case something was weird
d1.columns = ['ResidueA', 'ResidueB', 'AvgIntTot', 'AvgStdDev']
d2.columns = ['ResidueA', 'ResidueB', 'AvgIntTot', 'AvgStdDev']

## Change font size
plt.rcParams.update({'font.size': 24})
plt.rcParams.update({'figure.autolayout': True})

## X-axis labels (because you need to match the biologists)
labelsx2 = [100, 150, 200, 250, 300, 350, 450, '     BR', '', 1000, 1050]
## X-axis locations
placesx2 = [0, 50, 100, 150, 200, 250, 300, 325, 335, 350, 400]

n_residues = 450
index = np.arange(n_residues)
bar_width = 0.5

fig = plt.gcf()
#fig.set_size_inches(11,8.5)
fig.set_size_inches(22,17)

#111 = 1 row x 1 column x 1 index
ax = plt.subplot(111)
ax.axes.get_xaxis()
ax.set_xticks(index + bar_width / 2)
ax.set_xticks(placesx2)
ax.set_xticklabels(labelsx2, fontdict=None, minor=False)
## 0 to (1 + n_residues)
#ax.set_xlim(0,451)
ax.set_xlim(0,456)
ax.tick_params(axis='both', which='major', pad=10, length=10)
## #0000FA is blue, #FA7D00 is orange
sys_a = ax.bar(index, d1['AvgIntTot'], bar_width, color='#0000FA', \
 align='center', label=sys_a_tag)
sys_b = ax.bar(index + bar_width, d2['AvgIntTot'], bar_width, color='#FA7D00', \
 align='center',label=sys_b_tag)
plt.ylabel('Energy (kcal/mol)')
plt.xlabel('Residue Number')
plt.legend(handles=[sys_a,sys_b])
plt.savefig(save_name)
plt.close(save_name)
