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
d1 = pd.read_csv('WT-MUTA_total_interaction_res222_avg.dat', \
 delim_whitespace=True, header=0)
d2 = pd.read_csv('WT-MUTB_protein_system_total_interaction_res222_avg.dat', \
 delim_whitespace=True, header=0)
d3 = pd.read_csv('WT-MUTC_protein_system_total_interaction_res222_avg.dat', \
 delim_whitespace=True, header=0)
d4 = pd.read_csv('WT-MUTD_protein_system_total_interaction_res222_avg.dat', \
 delim_whitespace=True, header=0)
d5 = pd.read_csv('WT-MUTE_protein_system_total_interaction_res222_avg.dat', \
 delim_whitespace=True, header=0)

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
## Span goes from residue before SNP to residue after SNP -- 130
## r is red to matplotlib
ax1.axvspan(129, 131, alpha=0.2, color='gray')
ax1.bar(d1['Residue'], d1['DiffE'], color = 'r', label='WT - MUTA')
ax1.legend(loc="lower left")

## Span goes from residue before SNP to residue after SNP -- 130
ax2.axvspan(129, 131, alpha=0.2, color='gray')
ax2.bar(d2['Residue'], d2['DiffE'], color = 'orange', label='WT - MUTB')
ax2.legend(loc="lower left")

## Span goes from residue before SNP to residue after SNP -- 145
## g is green to matplotlib
ax3.axvspan(144, 146, alpha=0.2, color='gray')
ax3.bar(d3['Residue'], d3['DiffE'], color = 'g', label='WT - MUTC')
ax3.legend(loc="lower left")

## Span goes from residue before SNP to residue after SNP -- 275
## b is blue to matplotlib
ax4.axvspan(274, 276, alpha=0.2, color='gray')
ax4.bar(d4['Residue'], d4['DiffE'], color = 'b', label='WT - MUTD')
ax4.legend(loc="lower left")

## Span goes from residue before SNP to residue after SNP -- 375
ax5.axvspan(374, 376, alpha=0.2, color='gray')
ax5.bar(d5['Residue'], d5['DiffE'], color = 'purple', label='WT - MUTE')
ax5.legend(loc="lower left")

## Make subplots close to each other and hide x ticks for all but bottom plot
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

## 0 to total residues
plt.xlim(0, 450)

## X-axis labels (because you need to match the biologists)
labelsx2 = [100, 150, 200, 250, 300, 350, 450, '     BR', '', 1000, 1050]
## X-axis locations
placesx2 = [0, 50, 100, 150, 200, 250, 300, 325, 335, 350, 400]

## Put those xticks on the plot
ax1.set_xticks(placesx2)
ax1.set_xticklabels(labelsx2, fontdict=None, minor=False)

# ## Add a ytick
# extratick=[25]
# plt.yticks(list(plt.yticks()[0]) + extratick)

plt.savefig(save_name, dpi=300, bbox_inches = 'tight', pad_inches = 0)
