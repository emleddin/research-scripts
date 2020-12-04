#! /usr/bin/python3
## This script works with data generated from `rmagic-EDA-avg-diffs.r`

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

## Figure name
save_name = 'WT-MUTA_EDA_diff.png'

## Label on plot
sys_lab = "WT - MUTA"

## Threshold val
thresh=10

## Residue of interest to highlight with gray bar
ROI=150

## Total Residues (set x-axis from 0 to tot_res)
tot_res=450

## Read in data
## header = 0 reads header in first row because Python starts at 0
d1 = pd.read_csv('WT-MUTA_total_interaction_res220_avg.dat', \
 delim_whitespace=True, header=0)

## X-axis labels (because you need to match the biologists)
labelsx2 = [100, 150, 200, 250, 300, 350, 450, '     BR', '', 1000, 1050]
## X-axis locations
placesx2 = [0, 50, 100, 150, 200, 250, 300, 325, 335, 350, 400]

#---------- No need to modify past the curtain ---------------#

neg_thresh=np.negative(thresh)

## Define the column names to use throughout
## This renames the columns in case something was weird
d1.columns = ['Residue', 'DiffE', 'AvgSTDEV']

## Change font size
plt.rcParams.update({'font.size': 24})
plt.rcParams.update({'figure.autolayout': True})

fig = plt.gcf()
#fig.set_size_inches(11,8.5) WxH
fig.set_size_inches(22,17)
ax = plt.subplot(111)

## Add the subplots
## Span goes from residue before ROI to residue after ROI
## g is green to matplotlib
ax.axvspan(ROI-1, ROI+1, alpha=0.2, color='gray')
sys = ax.bar(d1['Residue'], d1['DiffE'], color = 'g', label=sys_lab)

##
ax.axhline(y=thresh, color='gray', linestyle='--', dashes=(5, 5))
ax.axhline(y=neg_thresh, color='gray', linestyle='--', dashes=(5, 5))

## Place the legend
ax.legend(loc="lower left")

## Put those xticks on the plot
ax.set_xticks(placesx2)
ax.set_xticklabels(labelsx2, fontdict=None, minor=False)

## Set x-axis
ax.set_xlim(0, tot_res)
ax.tick_params(axis='both', which='major', pad=10, length=10)

plt.ylabel('Energy (kcal/mol)')
plt.xlabel('Residue Number')
plt.legend(handles=[sys])

#plt.savefig(save_name, dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.savefig(save_name, dpi=300)
plt.close(save_name)
