#! /usr/bin/python3
## This script works with data generated from running `rmagic-EDA-avg-diffs.r`
## or on 2 different systems

## Needed for centering error bars correctly
## See https://github.com/matplotlib/matplotlib/issues/3400
import matplotlib as mpl
mpl.use( "cairo" )

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

## Figure file name
save_name = 'WT-SNP_EDA_scatter.png'

## Labels on plot
sys1_lab = 'Sys A'
sys2_lab = 'Sys B'
sys3_lab = 'Sys C'
sys4_lab = 'Sys D'
sys5_lab = 'Sys E'

## Color for the barplots -- use Matplotlib color names or a hex code
## These were chosen by https://personal.sron.nl/~pault/
color_1 = '#EE3377'
color_2 = '#EE7733'
color_3 = '#009988'
color_4 = '#33BBEE'
color_5 = '#0077BB'

## Chose markers
mark_1 = 'o'
mark_2 = 's'
mark_3 = 'X'
mark_4 = 'p'
mark_5 = 'v'

## Size of markers
marker_size = 15

## Chose variants
ROI_1 = 130
#ROI_2 = 130
ROI_3 = 145
ROI_4 = 275
ROI_5 = 375

## Threshold value -- places a horizontal line on the plot at +/- thresh
thresh=1

## Total Residues (set x-axis from 0 to n_res)
## In this case, nres + 5 to get a good cushion
tot_res = 455

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

## Set y-limits
ylim2 = [-20,40]
ylim = [-60,-40]

# ## X-axis labels (because you need to match the biologists)
labelsx2 = [100, 150, 200, 250, 300, 350, 400, '     BR', '', 1000, 1050, 2000]
# ## X-axis locations
placesx2 = [0, 50, 100, 150, 200, 250, 300, 325, 335, 350, 400, 450]

#---------- No need to modify past the curtain ---------------#

neg_thresh=np.negative(thresh)

## This renames the columns in case something was weird
d1.columns = ['Residue', 'DiffE', 'AvgStdDev']
d2.columns = ['Residue', 'DiffE', 'AvgStdDev']
d3.columns = ['Residue', 'DiffE', 'AvgStdDev']
d4.columns = ['Residue', 'DiffE', 'AvgStdDev']
d5.columns = ['Residue', 'DiffE', 'AvgStdDev']

# ## If df['DiffE'] value is less than threshold, set it to NaN
d1['DiffE'] = np.where(d1['DiffE'].between(neg_thresh,thresh), float('nan'),\
 d1['DiffE'])
d1['AvgStdDev'] = np.where(d1['DiffE'].between(neg_thresh,thresh), float('nan'),\
 d1['AvgStdDev'])
d2['DiffE'] = np.where(d2['DiffE'].between(neg_thresh,thresh), float('nan'),\
 d2['DiffE'])
d2['AvgStdDev'] = np.where(d2['DiffE'].between(neg_thresh,thresh), float('nan'),\
 d2['AvgStdDev'])
d3['DiffE'] = np.where(d3['DiffE'].between(neg_thresh,thresh), float('nan'),\
 d3['DiffE'])
d3['AvgStdDev'] = np.where(d3['DiffE'].between(neg_thresh,thresh), float('nan'),\
 d3['AvgStdDev'])
d4['DiffE'] = np.where(d4['DiffE'].between(neg_thresh,thresh), float('nan'),\
 d4['DiffE'])
d4['AvgStdDev'] = np.where(d4['DiffE'].between(neg_thresh,thresh), float('nan'),\
 d4['AvgStdDev'])
d5['DiffE'] = np.where(d5['DiffE'].between(neg_thresh,thresh), float('nan'),\
 d5['DiffE'])
d5['AvgStdDev'] = np.where(d5['DiffE'].between(neg_thresh,thresh), float('nan'),\
 d5['AvgStdDev'])

## Change font size
plt.rcParams.update({'font.size': 24})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'lines.markersize': marker_size})

fig = plt.figure()
#fig.set_size_inches(11,8.5) WxH
fig.set_size_inches(22,17)

## https://stackoverflow.com/questions/17976103/matplotlib-broken-axis-example-uneven-subplot-size
## Set up Gridspec y-axes
ylim2_rat = (ylim2[1]-ylim2[0])/(ylim[1]-ylim[0]+ylim2[1]-ylim2[0])
ylim_rat = (ylim[1]-ylim[0])/(ylim[1]-ylim[0]+ylim2[1]-ylim2[0])

# spec = mpl.gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[3,1])
spec = mpl.gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[ylim2_rat,ylim_rat])
ax2 = fig.add_subplot(spec[0])
ax1 = fig.add_subplot(spec[1])

## Plot on bottom y-axis
ax1.errorbar(d1['Residue'], d1['DiffE'], yerr=d1['AvgStdDev'],\
 fmt=mark_1, color=color_1, label=sys1_lab, ecolor='black', capsize=2,\
 capthick=0.5, elinewidth=0.5)
ax1.errorbar(d2['Residue'], d2['DiffE'], yerr=d2['AvgStdDev'],\
 fmt=mark_2, color=color_2, label=sys2_lab, ecolor='black', capsize=2,\
 capthick=0.5, elinewidth=0.5)
ax1.errorbar(d3['Residue'], d3['DiffE'], yerr=d3['AvgStdDev'],\
 fmt=mark_3, color=color_3, label=sys3_lab, ecolor='black', capsize=2,\
 capthick=0.5, elinewidth=0.5)
ax1.errorbar(d4['Residue'], d4['DiffE'], yerr=d4['AvgStdDev'],\
 fmt=mark_4, color=color_4, label=sys4_lab, ecolor='black', capsize=2,\
 capthick=0.5, elinewidth=0.5)
ax1.errorbar(d5['Residue'], d5['DiffE'], yerr=d5['AvgStdDev'],\
 fmt=mark_5, color=color_5, label=sys5_lab, ecolor='black', capsize=2,\
 capthick=0.5, elinewidth=0.5)

## Plot on top y-axis
ax2.errorbar(d1['Residue'], d1['DiffE'], yerr=d1['AvgStdDev'],\
 fmt=mark_1, color=color_1, label=sys1_lab, ecolor='black', capsize=2,\
 capthick=0.5, elinewidth=0.5)
ax2.errorbar(d2['Residue'], d2['DiffE'], yerr=d2['AvgStdDev'],\
 fmt=mark_2, color=color_2, label=sys2_lab, ecolor='black', capsize=2,\
 capthick=0.5, elinewidth=0.5)
ax2.errorbar(d3['Residue'], d3['DiffE'], yerr=d3['AvgStdDev'],\
 fmt=mark_3, color=color_3, label=sys3_lab, ecolor='black', capsize=2,\
 capthick=0.5, elinewidth=0.5)
ax2.errorbar(d4['Residue'], d4['DiffE'], yerr=d4['AvgStdDev'],\
 fmt=mark_4, color=color_4, label=sys4_lab, ecolor='black', capsize=2,\
 capthick=0.5, elinewidth=0.5)
ax2.errorbar(d5['Residue'], d5['DiffE'], yerr=d5['AvgStdDev'],\
 fmt=mark_5, color=color_5, label=sys5_lab, ecolor='black', capsize=2,\
 capthick=0.5, elinewidth=0.5)

## Set y-limits
ax1.set_ylim(ylim)
ax2.set_ylim(ylim2)

## Set y-axes every 10. Use +1 because arange isn't inclusive.
ax1.set_yticks(np.arange(ylim[0], ylim[1]+1, 10))
ax2.set_yticks(np.arange(ylim2[0], ylim2[1]+1, 10))

## Set up edge visibility
ax1.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.xaxis.tick_top()
ax2.tick_params(labeltop=False, top=False)
ax1.xaxis.tick_bottom()

## Put those xticks on the plot
ax1.set_xticks(placesx2)
ax1.set_xticklabels(labelsx2, fontdict=None, minor=False)
## 0 to total residues
ax1.set_xlim(0, tot_res)
ax2.set_xlim(0, tot_res)

## Add threshold lines
ax2.axhline(y=thresh, color='gray', linestyle='--', dashes=(5, 5))
ax2.axhline(y=neg_thresh, color='gray', linestyle='--', dashes=(5, 5))

## Add variant places
ax1.axvspan(ROI_1-1, ROI_1+1, alpha=0.2, color='gray')
#ax1.axvspan(ROI_2-1, ROI_2+1, alpha=0.2, color='gray')
ax1.axvspan(ROI_3-1, ROI_3+1, alpha=0.2, color='gray')
ax1.axvspan(ROI_4-1, ROI_4+1, alpha=0.2, color='gray')
ax1.axvspan(ROI_5-1, ROI_5+1, alpha=0.2, color='gray')

ax2.axvspan(ROI_1-1, ROI_1+1, alpha=0.2, color='gray')
#ax2.axvspan(ROI_2-1, ROI_2+1, alpha=0.2, color='gray')
ax2.axvspan(ROI_3-1, ROI_3+1, alpha=0.2, color='gray')
ax2.axvspan(ROI_4-1, ROI_4+1, alpha=0.2, color='gray')
ax2.axvspan(ROI_5-1, ROI_5+1, alpha=0.2, color='gray')


## Add the diagonal marks at the break point. Do some magic because of gridspec
angle_mod = 0.003
length_mod = 0.001
#
kwargs = dict(color='k', clip_on=False)
xlim = ax1.get_xlim()
dx = angle_mod*(xlim[1]-xlim[0])
dy = length_mod*(ylim2[1]-ylim2[0])/ylim2_rat
ax2.plot((xlim[0]-dx,xlim[0]+dx), (ylim2[0]-dy,ylim2[0]+dy), **kwargs)
ax2.plot((xlim[1]-dx,xlim[1]+dx), (ylim2[0]-dy,ylim2[0]+dy), **kwargs)
dy = length_mod*(ylim[1]-ylim[0])/ylim_rat
ax1.plot((xlim[0]-dx,xlim[0]+dx), (ylim[1]-dy,ylim[1]+dy), **kwargs)
ax1.plot((xlim[1]-dx,xlim[1]+dx), (ylim[1]-dy,ylim[1]+dy), **kwargs)
ax1.set_xlim(xlim)
ax2.set_xlim(xlim)

## Set x- and y-axis labels. Use text for y because of broken axis
plt.xlabel('Residue Number', labelpad=20)
fig.text(0, 0.5, 'Energy (kcal/mol)', va='center', rotation='vertical')

## Add two-column legend
ax1.legend(loc='lower center', ncol=2)

#plt.savefig(save_name, dpi=300, bbox_inches = 'tight', pad_inches = 0)
plt.savefig(save_name, dpi=300)
plt.close(save_name)
