#! /usr/bin/python3
## This script works with data generated from `rmagic-EDA-avg-diffs.r`

## Needed on MacOS for centering error bars correctly
## See https://github.com/matplotlib/matplotlib/issues/3400
#import matplotlib as mpl
#mpl.use( "cairo" )

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

## Figure file name
save_name = 'WT-MUTA_EDA_diff.png'

## Label on plot
sys_lab = "WT - MUTA"

## Threshold value -- places a horizontal line on the plot at +/- thresh
thresh=1

## Residue of interest to highlight with gray bar
ROI=150

## Total Residues (set x-axis from 0 to tot_res)
tot_res=450

## Plot with error bars? True = yes
error_bars=True

## Color for the barplot -- use Matplotlib color names or a hex code
bar_color='g'

## Read in data
## header = 0 reads header in first row because Python starts at 0
d1 = pd.read_csv('WT-MUTA_total_interaction_res220_avg.dat', \
 delim_whitespace=True, header=0)

## X-axis labels (because you need to match the biologists)
labelsx2 = [100, 150, 200, 250, 300, 350, 400, '     BR', '', 1000, 1050]
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

## No error bars
if error_bars == False:
    sys = ax.bar(d1['Residue'], d1['DiffE'], color = bar_color, label=sys_lab)
## Yes error bars
else:
    e_bar = {'ecolor': 'black',
             'capsize': 2,
             'capthick': 0.5,
             'elinewidth': 0.5,
             'label': 'Avg. St. Dev.'}

    sys = ax.bar(d1['Residue'], d1['DiffE'], yerr=d1['AvgSTDEV'], \
     color = bar_color, label=sys_lab, error_kw=e_bar)

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

plt.ylabel('$\Delta\Delta$E (kcal/mol)')
plt.xlabel('Residue Number')

if error_bars == False:
    plt.legend(handles=[sys])
else:
    ## Have plot label above error bars in legend
    handles,labels = ax.get_legend_handles_labels()
    handles = [handles[1], handles[0]]
    labels = [labels[1], labels[0]]
    plt.legend(handles, labels)

#plt.savefig(save_name, dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.savefig(save_name, dpi=300)
plt.close(save_name)
