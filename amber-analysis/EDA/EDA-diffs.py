#! /usr/bin/python3
## This script works with data generated from multiple comparisons
## of `rmagic-EDA-avg-diffs.r`

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

## Figure name
save_name = 'WT-SNP_EDA_diffs.png'

## Read in data
## header = 0 reads header in first row because Python starts at 0
d1 = pd.read_csv('WT-MUTA_protein_system_total_interaction_res220_avg.dat', \
 delim_whitespace=True, header=0)
d2 = pd.read_csv('WT-MUTB_protein_system_total_interaction_res220_avg.dat', \
 delim_whitespace=True, header=0)
d3 = pd.read_csv('WT-MUTC_protein_system_total_interaction_res220_avg.dat', \
 delim_whitespace=True, header=0)
d4 = pd.read_csv('WT-MUTD_protein_system_total_interaction_res220_avg.dat', \
 delim_whitespace=True, header=0)
d5 = pd.read_csv('WT-MUTE_protein_system_total_interaction_res220_avg.dat', \
 delim_whitespace=True, header=0)

## Labels for each set of data
sys1_lab = "WT - MUTA"
sys2_lab = "WT - MUTB"
sys3_lab = "WT - MUTC"
sys4_lab = "WT - MUTD"
sys5_lab = "WT - MUTE"

## Set the residues of interest to highlight with a gray bar for each plot
ROI_1=130
ROI_2=130
ROI_3=145
ROI_4=275
ROI_5=375

## Total Residues (set x-axis from 0 to tot_res)
tot_res=450

## X-axis labels (because you need to match the biologists)
labelsx2 = [100, 150, 200, 250, 300, 350, 450, '     BR', '', 1000, 1050]
## X-axis locations
placesx2 = [0, 50, 100, 150, 200, 250, 300, 325, 335, 350, 400]

#---------- No need to modify past the curtain ---------------#

## Define the column names to use throughout
## This renames the columns in case something was weird
d1.columns = ['Residue', 'DiffE', 'AvgSTDEV']
d2.columns = ['Residue', 'DiffE', 'AvgSTDEV']
d3.columns = ['Residue', 'DiffE', 'AvgSTDEV']
d4.columns = ['Residue', 'DiffE', 'AvgSTDEV']
d5.columns = ['Residue', 'DiffE', 'AvgSTDEV']

## Resent legend and tick font
plt.rc('legend', fontsize=5)
plt.rc('xtick', labelsize=5)    # fontsize of the tick labels
plt.rc('ytick', labelsize=5)    # fontsize of the tick labels

## Size in inches WxH
f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=True, sharey=True, \
 figsize=(7.25, 10))

## Add the subplots
## Span goes from res before MUT to res after MUT; alpha is transparency
## r is red to matplotlib
ax1.axvspan(ROI_A-1, ROI_A+1, alpha=0.2, color='gray')
ax1.bar(d1['Residue'], d1['DiffE'], color = 'r', label=sys1_lab)
ax1.legend(loc="lower left")

ax2.axvspan(ROI_B-1, ROI_B+1, alpha=0.2, color='gray')
ax2.bar(d2['Residue'], d2['DiffE'], color = 'orange', label=sys2_lab)
ax2.legend(loc="lower left")

## g is green to matplotlib
ax3.axvspan(ROI_C-1, ROI_C+1, alpha=0.2, color='gray')
ax3.bar(d3['Residue'], d3['DiffE'], color = 'g', label=sys3_lab)
ax3.legend(loc="lower left")

## b is blue to matplotlib
ax4.axvspan(ROI_D-1, ROI_D+1, alpha=0.2, color='gray')
ax4.bar(d4['Residue'], d4['DiffE'], color = 'b', label=sys4_lab)
ax4.legend(loc="lower left")

ax5.axvspan(ROI_E-1, ROI_E+1, alpha=0.2, color='gray')
ax5.bar(d5['Residue'], d5['DiffE'], color = 'purple', label=sys5_lab)
ax5.legend(loc="lower left")

## Make subplots close to each other and hide x ticks for all but bottom plot
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

## Set x-axis
plt.xlim(0, tot_res)

## Put those xticks on the plot
ax1.set_xticks(placesx2)
ax1.set_xticklabels(labelsx2, fontdict=None, minor=False)

# ## Add a ytick
# extratick=[25]
# plt.yticks(list(plt.yticks()[0]) + extratick)

plt.savefig(save_name, dpi=300, bbox_inches = 'tight', pad_inches = 0)
