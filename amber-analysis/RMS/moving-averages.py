#! /usr/bin/python3
## This script works with data generated from cpptraj with `cpptraj_analysis.in`

## Needed for centering error bars correctly
## See https://github.com/matplotlib/matplotlib/issues/3400
import matplotlib as mpl
mpl.use( "cairo" )

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import make_interp_spline, BSpline

## Figure file names
## Note: RMSF is *not* actually a moving average!!!
save_name_RMSD = 'test-moving-avgs-rmsd.png'
save_name_HB = 'test-moving-avgs-hbond.png'
save_name_RMSF = 'test-moving-avgs-rmsf.png'

## Labels on plot
sys1_lab = 'WT'
sys2_lab = 'Mutation A'
sys3_lab = 'Mutation B'
sys4_lab = 'Mutation C'

## Color for the barplots -- use Matplotlib color names or a hex code
## These were chosen by https://personal.sron.nl/~pault/
color_1 = '#EE3377'
color_2 = '#EE7733'
color_3 = '#009988'
color_4 = '#33BBEE'

## Read in RMSD Data
## header = 0 reads header in first row because Python starts at 0
d1rmsd = pd.read_csv('WT-system/WT_protein_system_total_bb_rms.dat', \
 delim_whitespace=True, header=0)
d2rmsd = pd.read_csv('MUT-A-system/MUT_A_system_total_bb_rms.dat', \
 delim_whitespace=True, header=0)
d3rmsd = pd.read_csv('MUT-B-system/MUT_B_system_RNA_total_bb_rms.dat', \
 delim_whitespace=True, header=0)
d4rmsd = pd.read_csv('MUT-C-system/MUT_C_system_total_bb_rms.dat', \
 delim_whitespace=True, header=0)

## Read in HB Data
## header = 0 reads header in first row because Python starts at 0
d1hb = pd.read_csv('WT-system/WT_protein_system_hbond.dat', \
 delim_whitespace=True, header=0)
d2hb = pd.read_csv('MUT-A-system/MUT_A_system_hbond.dat', \
 delim_whitespace=True, header=0)
d3hb = pd.read_csv('MUT-B-system/MUT_B_system_hbond.dat', \
 delim_whitespace=True, header=0)
d4hb = pd.read_csv('MUT-C-system/MUT_C_system_hbond.dat', \
 delim_whitespace=True, header=0)

## Read in RMSF Data
## header = 0 reads header in first row because Python starts at 0
d1rmsf = pd.read_csv('WT-system/WT_protein_system_rmsf_byres.dat', \
 delim_whitespace=True, header=0)
d2rmsf = pd.read_csv('MUT-A-system/MUT_A_system_rmsf_byres.dat', \
 delim_whitespace=True, header=0)
d3rmsf = pd.read_csv('MUT-B-system/MUT_B_system_rmsf_byres.dat', \
 delim_whitespace=True, header=0)
d4rmsf = pd.read_csv('MUT-C-system/MUT_C_system_rmsf_byres.dat', \
 delim_whitespace=True, header=0)

## Summarize the datasets, labels, and colors into arrays
datasets_rmsd = [d1rmsd, d2rmsd, d3rmsd, d4rmsd]
datasets_hb = [d1hb, d2hb, d3hb, d4hb]
datasets_rmsf = [d1rmsf, d2rmsf, d3rmsf, d4rmsf]
labels = [sys1_lab, sys2_lab, sys3_lab, sys4_lab]
colors = [color_1, color_2, color_3, color_4]

## Fraction to multiply to convert frames to nanoseconds
## Berendsen inputs are 1/500, Langevin inputs are 1/100
frame2ns = 1/500

#---------- No need to modify past the curtain ---------------#

## This renames the columns in case something was weird
d1rmsd.columns = ['Frame', 'RMSD', 'AvgStdDev']
d2rmsd.columns = ['Frame', 'RMSD', 'AvgStdDev']
d3rmsd.columns = ['Frame', 'RMSD', 'AvgStdDev']
d4rmsd.columns = ['Frame', 'RMSD', 'AvgStdDev']

## Plot at end because of sandwich effect
i = 0
for ds in datasets_rmsd:
    ## Change Frame to ns
    ds['Time'] = ds['Frame'] * frame2ns
    rolling_mean = ds['RMSD'].rolling(window=int(1/frame2ns)).mean()
    plt.plot(ds['Time'], rolling_mean, label=(labels[i] + ', 1 ns SMA'),
     color=colors[i])
    i += 1

## Add two-column legend
# plt.legend(loc='best', ncol=2)
plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2)

plt.xlabel('Time (ns)', labelpad=5)
plt.ylabel('RMSD ($\AA$)', labelpad=10)

#plt.savefig(save_name, dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.savefig(save_name_RMSD, dpi=300)
plt.close(save_name_RMSD)
plt.clf()

#----------------- HB
## This renames the columns in case something was weird
d1hb.columns = ['Frame', 'HB', 'AvgStdDev']
d2hb.columns = ['Frame', 'HB', 'AvgStdDev']
d3hb.columns = ['Frame', 'HB', 'AvgStdDev']
d4hb.columns = ['Frame', 'HB', 'AvgStdDev']

## Plot at end because of sandwich effect
i = 0
for ds in datasets_hb:
    ## Change Frame to ns
    ds['Time'] = ds['Frame'] * frame2ns
    rolling_mean = ds['HB'].rolling(window=int(1/frame2ns)).mean()
    plt.plot(ds['Time'], rolling_mean, label=(labels[i] + ', 1 ns SMA'),
     color=colors[i])
    i += 1

## Add two-column legend
# plt.legend(loc='best', ncol=2)
plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2)

plt.xlabel('Time (ns)', labelpad=5)
plt.ylabel('Number of hydrogen bonds', labelpad=10)

#plt.savefig(save_name, dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.savefig(save_name_HB, dpi=300)
plt.close(save_name_HB)
plt.clf()

#----------------- RMSF
## This renames the columns in case something was weird
d1rmsf.columns = ['Residue', 'RMSF', 'AvgStdDev']
d2rmsf.columns = ['Residue', 'RMSF', 'AvgStdDev']
d3rmsf.columns = ['Residue', 'RMSF', 'AvgStdDev']
d4rmsf.columns = ['Residue', 'RMSF', 'AvgStdDev']

## Don't plot as moving averages because not time-based
i = 0
for ds in datasets_rmsf:
    plt.plot(ds['Residue'], ds['RMSF'], label=(labels[i]), color=colors[i])
    i += 1

## Add two-column legend
# plt.legend(loc='best', ncol=2)
plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2)

plt.xlabel('Time (ns)', labelpad=5)
plt.ylabel('RMSF ($\AA$)', labelpad=10)

#plt.savefig(save_name, dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.savefig(save_name_RMSF, dpi=300)
plt.close(save_name_RMSF)
