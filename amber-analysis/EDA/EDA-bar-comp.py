#! /usr/bin/python3
## This script works with data generated from running `rmagic-EDA-avg.r` or
## `rmagic-EDA-single-run.r` on 2 different systems

## Needed on MacOS for centering error bars correctly
## See https://github.com/matplotlib/matplotlib/issues/3400
#import matplotlib as mpl
#mpl.use( "cairo" )

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

## Figure file name
save_name = 'sysA_sysB_total_barchart.png'

## Labels on plot
sys_a_tag = 'Sys A'
sys_b_tag = 'Sys B'

## Plot with error bars? True = yes
error_bars=True

## Color for the barplots -- use Matplotlib color names or a hex code
## #0000FA is blue, #FA7D00 is orange
color_a = '#0000FA'
color_b = '#FA7D00'

## Threshold value -- places a horizontal line on the plot at +/- thresh
thresh=1

## Total Residues (set x-axis from 0 to n_res)
n_residues = 450

## Read in data
## header = 0 reads header in first row because Python starts at 0
d1 = pd.read_csv('SysA_protein_system_res220_tot_avg.dat', \
 delim_whitespace=True, header=0)
d2 = pd.read_csv('SysB_protein_system_res220_tot_avg.dat', \
 delim_whitespace=True, header=0)

## X-axis labels (because you need to match the biologists)
labelsx2 = [100, 150, 200, 250, 300, 350, 450, '     BR', '', 1000, 1050]
## X-axis locations
placesx2 = [0, 50, 100, 150, 200, 250, 300, 325, 335, 350, 400]

#---------- No need to modify past the curtain ---------------#

neg_thresh=np.negative(thresh)

## This renames the columns in case something was weird
d1.columns = ['ResidueA', 'ResidueB', 'AvgIntTot', 'AvgStdDev']
d2.columns = ['ResidueA', 'ResidueB', 'AvgIntTot', 'AvgStdDev']

## Change font size
plt.rcParams.update({'font.size': 24})
plt.rcParams.update({'figure.autolayout': True})

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

ax.set_xlim(0,n_residues+1)
ax.tick_params(axis='both', which='major', pad=10, length=10)

## No error bars
if error_bars == False:
    sys_a = ax.bar(index, d1['AvgIntTot'], bar_width, color='#0000FA', \
     align='center', label=sys_a_tag)
    sys_b = ax.bar(index + bar_width, d2['AvgIntTot'], bar_width, color='#FA7D00', \
     align='center',label=sys_b_tag)
else:
    e_bar = {'ecolor': 'black',
             'capsize': 2,
             'capthick': 0.5,
             'elinewidth': 0.5,
             'label': 'Avg. St. Dev.'}
    sys_a = ax.bar(index, d1['AvgIntTot'], bar_width, yerr=d1['AvgStdDev'], \
     color='#0000FA', align='center', label=sys_a_tag, error_kw=e_bar)
    sys_b = ax.bar(index + bar_width, d2['AvgIntTot'], bar_width, \
     yerr=d1['AvgStdDev'], color='#FA7D00', align='center',label=sys_b_tag, \
     error_kw=e_bar)

## Add threshold lines
ax.axhline(y=thresh, color='gray', linestyle='--', dashes=(5, 5))
ax.axhline(y=neg_thresh, color='gray', linestyle='--', dashes=(5, 5))

plt.ylabel('Energy (kcal/mol)')
plt.xlabel('Residue Number')

if error_bars == False:
    plt.legend(handles=[sys_a,sys_b])
else:
    ## Have plot label above error bars in legend
    handles,labels = ax.get_legend_handles_labels()
    handles = [handles[1], handles[3], handles[0]]
    labels = [labels[1], labels[3], labels[0]]
    plt.legend(handles, labels)

#plt.savefig(save_name, dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.savefig(save_name, dpi=300)
plt.close(save_name)
