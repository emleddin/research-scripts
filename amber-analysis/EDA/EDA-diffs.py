#! /usr/bin/python3
## This script works with data generated from multiple comparisons
## of `rmagic-EDA-avg-diffs.r`

## Needed on MacOS for centering error bars correctly
## See https://github.com/matplotlib/matplotlib/issues/3400
#import matplotlib as mpl
#mpl.use( "cairo" )

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

## Color for the barplots -- use Matplotlib color names or a hex code
## r is red, g is green, b is blue
color_1 = 'r'
color_2 = 'orange'
color_3 = 'g'
color_4 = 'b'
color_5 = 'purple'

## Set the residues of interest to highlight with a gray bar for each plot
ROI_1=130
ROI_2=130
ROI_3=145
ROI_4=275
ROI_5=375

## Total Residues (set x-axis from 0 to tot_res)
tot_res=450

## Threshold value -- places a horizontal line on the plot at +/- thresh
thresh=1

## Plot with error bars? True = yes
error_bars=True

## X-axis labels (because you need to match the biologists)
labelsx2 = [100, 150, 200, 250, 300, 350, 450, '     BR', '', 1000, 1050]
## X-axis locations
placesx2 = [0, 50, 100, 150, 200, 250, 300, 325, 335, 350, 400]

#---------- No need to modify past the curtain ---------------#
neg_thresh=np.negative(thresh)

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
ax1.axvspan(ROI_1-1, ROI_1+1, alpha=0.2, color='gray')
ax1.bar(d1['Residue'], d1['DiffE'], color = color_1, label=sys1_lab)
if error_bars == True:
    ax1.errorbar(d1['Residue'], d1['DiffE'], xerr=None, yerr=d1['AvgSTDEV'], \
     ecolor='black', capsize=0.5, capthick=0.25, elinewidth=0.25, \
     label='Avg. St. Dev.', ls='none')
ax1.legend(loc="lower left")
## Add threshold lines
ax1.axhline(y=thresh, color='gray', linestyle='--', dashes=(5, 5))
ax1.axhline(y=neg_thresh, color='gray', linestyle='--', dashes=(5, 5))

ax2.axvspan(ROI_2-1, ROI_2+1, alpha=0.2, color='gray')
ax2.bar(d2['Residue'], d2['DiffE'], color = color_2, label=sys2_lab)
if error_bars == True:
    ax2.errorbar(d2['Residue'], d2['DiffE'], xerr=None, yerr=d2['AvgSTDEV'], \
     ecolor='black', capsize=0.5, capthick=0.25, elinewidth=0.25, \
     label='Avg. St. Dev.', ls='none')
ax2.legend(loc="lower left")
## Add threshold lines
ax2.axhline(y=thresh, color='gray', linestyle='--', dashes=(5, 5))
ax2.axhline(y=neg_thresh, color='gray', linestyle='--', dashes=(5, 5))

ax3.axvspan(ROI_3-1, ROI_3+1, alpha=0.2, color='gray')
ax3.bar(d3['Residue'], d3['DiffE'], color = color_3, label=sys3_lab)
if error_bars == True:
    ax3.errorbar(d3['Residue'], d3['DiffE'], xerr=None, yerr=d3['AvgSTDEV'], \
     ecolor='black', capsize=0.5, capthick=0.25, elinewidth=0.25, \
     label='Avg. St. Dev.', ls='none')
ax3.legend(loc="lower left")
## Add threshold lines
ax3.axhline(y=thresh, color='gray', linestyle='--', dashes=(5, 5))
ax3.axhline(y=neg_thresh, color='gray', linestyle='--', dashes=(5, 5))

ax4.axvspan(ROI_4-1, ROI_4+1, alpha=0.2, color='gray')
ax4.bar(d4['Residue'], d4['DiffE'], color = color_4, label=sys4_lab)
if error_bars == True:
    ax4.errorbar(d4['Residue'], d4['DiffE'], xerr=None, yerr=d4['AvgSTDEV'], \
     ecolor='black', capsize=0.5, capthick=0.25, elinewidth=0.25, \
     label='Avg. St. Dev.', ls='none')
ax4.legend(loc="lower left")
## Add threshold lines
ax4.axhline(y=thresh, color='gray', linestyle='--', dashes=(5, 5))
ax4.axhline(y=neg_thresh, color='gray', linestyle='--', dashes=(5, 5))

ax5.axvspan(ROI_5-1, ROI_5+1, alpha=0.2, color='gray')
ax5.bar(d5['Residue'], d5['DiffE'], color = color_5, label=sys5_lab)
if error_bars == True:
    ax5.errorbar(d5['Residue'], d5['DiffE'], xerr=None, yerr=d5['AvgSTDEV'], \
     ecolor='black', capsize=0.5, capthick=0.25, elinewidth=0.25, \
     label='Avg. St. Dev.', ls='none')
ax5.legend(loc="lower left")
## Add threshold lines
ax5.axhline(y=thresh, color='gray', linestyle='--', dashes=(5, 5))
ax5.axhline(y=neg_thresh, color='gray', linestyle='--', dashes=(5, 5))

## Make subplots close to each other and hide x ticks for all but bottom plot
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

## Set x-axis
plt.xlim(0, tot_res)

## Add labels (use ax3 to center on y-axis)
ax3.set_ylabel('Energy (kcal/mol)')
plt.xlabel('Residue Number')

## Put those xticks on the plot
ax1.set_xticks(placesx2)
ax1.set_xticklabels(labelsx2, fontdict=None, minor=False)

# ## Add a ytick
# extratick=[25]
# plt.yticks(list(plt.yticks()[0]) + extratick)

plt.savefig(save_name, dpi=300, bbox_inches = 'tight', pad_inches = 0)
#plt.savefig(save_name, dpi=300)
plt.close(save_name)
