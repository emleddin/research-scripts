'''
This script is designed to take two correlation matrix data files produced
by cpptraj (part of the AmberTools suite) and provide multiple comparisons
between them.  Given the nature of how cpptraj represents correlated and
anti-correlated movement, it is necessary to consider the original values
as well as the resulting difference.
The code below produces four subplots in a single image file.  The top two
subplots show the degree of change in correlated movement (left) or anti-
correlated movement (right).  The bottom plots show cases where a residue
pair has switched between correlated and anti-correlated movement, and the
magnitude of this change.
In its current form, the script requires modification by the user in three
fields: two input files to compare and an output filename to produce.


Original: Mark Hix (GitHub: markahix)
https://github.com/markahix/Basic-Scripts/blob/master/cpptraj_plots/Corr_Decomp_Analysis.py
Update: Emmett (GitHub: emleddin)
'''

##################
# Library Import #
##################
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy.linalg as la
import statsmodels.api as sm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
from mpl_toolkits.axes_grid1 import AxesGrid

####################
# File Definitions #
####################
## The script does B - A
B_Filename="./MUT.dat"
A_Filename="./WT.dat"
outputfile='./MUT-WT_Correlation_Comparisons.png'

## Define your color map
## Bcolor is *more* B in B - A
## Acolor is *more* A in B - A
Bcolor = '#5E3C99'   # purple
midcolor = '#F7F7F7' # white
Acolor = '#E66101'   # orange

## Fix labels if you'd like
set_places = True

## Define `places` and `labels`
##  REQUIRED if `set_places = True`
## Explicitly define axis ticks/labels (e.g., to match biological numbering)
## Explicity choose where to put x and y ticks
places = [0, 100, 200, 300, 333, 347, 400, 430]

## Define those very x and y tick labels
labels = [1130, 1230, 1330, 1430, ' ', 1842, ' ', 'DNA']

## Should the script save an individual graph for each component?
## Note: To modify the save names and whether they have titles, go behind the
## curtain and modify the last ~15 lines of this script!!!
save_singles = True

##--------------------------------------------------------------------------##
##--------------------------- Behind the Curtain ---------------------------##
##--------------------------------------------------------------------------##

################################################################
# The shiftedColorMap function was taken without change from : #
# https://gist.github.com/phobson/7916777                      #
# Thanks to Paul Hobson!                                       #
################################################################

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.
    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = mpl.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

################
# File Loading #
################
Corr_A = np.genfromtxt(A_Filename, delimiter=None)
Corr_B = np.genfromtxt(B_Filename, delimiter=None)

Res_A=len(Corr_A)+1
Res_B=len(Corr_B)+1

#########################
# Matrix Initialization #
#########################
Corr_Diff=np.array(np.zeros(Res_A*Res_A,float)).reshape(Res_A,Res_A)
Anti_Diff=np.array(np.zeros(Res_A*Res_A,float)).reshape(Res_A,Res_A)
Corr_to_Anti=np.array(np.zeros(Res_A*Res_A,float)).reshape(Res_A,Res_A)
Anti_to_Corr=np.array(np.zeros(Res_A*Res_A,float)).reshape(Res_A,Res_A)

#####################
# Matrix Processing #
#####################
for i in range(Res_A-1):
    for j in range(Res_A-1):
        if float(Corr_A[i][j])>=0:
            if float(Corr_B[i][j])>=0:
                Corr_Diff[i][j]=float(Corr_B[i][j])-float(Corr_A[i][j])
            elif float(Corr_B[i][j])<0:
                Corr_to_Anti[i][j]=-float(Corr_B[i][j])+float(Corr_A[i][j])
        elif float(Corr_A[i][j])<0:
            if float(Corr_B[i][j])>=0:
                Anti_to_Corr[i][j]=float(Corr_B[i][j])-float(Corr_A[i][j])
            elif float(Corr_B[i][j])<0:
                Anti_Diff[i][j]=-float(Corr_B[i][j])+float(Corr_A[i][j])

# Setting minimum and maximum ranges to ensure multiple data sets processed via this script are visually comparable.
Corr_Diff[Res_A-1][Res_A-1]=1.0
Anti_Diff[Res_A-1][Res_A-1]=1.0
Corr_to_Anti[Res_A-1][Res_A-1]=2.0
Anti_to_Corr[Res_A-1][Res_A-1]=2.0

Corr_Diff[Res_A-2][Res_A-1]=-1.0
Anti_Diff[Res_A-2][Res_A-1]=-1.0

###################
# Figure Plotting #
###################

plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', titlepad=1)

plotted=plt.figure(1,figsize=(10,8),dpi=100)
corrplot=plt.subplot(2,2,1)
antiplot=plt.subplot(2,2,2)
corranti=plt.subplot(2,2,3)
anticorr=plt.subplot(2,2,4)

## Define the custom color map
test = mpl.colors.LinearSegmentedColormap.from_list('test',[Acolor,midcolor,Bcolor])
cm.register_cmap(cmap=test)

### Note: these will always have 1 extra row/column with the max!!
#corr_shifted_cmap=shiftedColorMap(cm.bwr_r,0.0,1 - Corr_Diff.max() / (Corr_Diff.max() + abs(Corr_Diff.min())),1.0,'shiftedcorr')
corr_shifted_cmap=shiftedColorMap(test,0.0,1 - Corr_Diff.max() / (Corr_Diff.max() + abs(Corr_Diff.min())),1.0,'shiftedcorr')
corrplotax=corrplot.matshow(Corr_Diff,cmap=corr_shifted_cmap, origin="lower")
corrplot.set_title('Correlated Movement')
if set_places == True:
    corrplot.axes.get_xaxis()
    corrplot.set_xticks(places)
    corrplot.set_xticklabels(labels, fontdict=None, minor=False)
corrplot.tick_params(axis='x', which='both', bottom=True, top=False, labeltop=False, labelbottom=True)
if set_places == True:
    corrplot.axes.get_yaxis()
    corrplot.set_yticks(places)
    corrplot.set_yticklabels(labels, fontdict=None, minor=False)
plotted.colorbar(corrplotax, ticks=[Corr_Diff.min(),0, Corr_Diff.max()],ax=corrplot)

#anti_shifted_cmap=shiftedColorMap(cm.bwr_r,0.0,1 - Anti_Diff.max() / (Anti_Diff.max() + abs(Anti_Diff.min())),1.0,'shiftedcorr')
anti_shifted_cmap=shiftedColorMap(test,0.0,1 - Anti_Diff.max() / (Anti_Diff.max() + abs(Anti_Diff.min())),1.0,'shiftedcorr')
antiplotax=antiplot.matshow(Anti_Diff,cmap=anti_shifted_cmap, origin="lower")
antiplot.set_title('Anticorrelated Movement')
if set_places == True:
    antiplot.axes.get_xaxis()
    antiplot.set_xticks(places)
    antiplot.set_xticklabels(labels, fontdict=None, minor=False)
antiplot.tick_params(axis='x', which='both', bottom=True, top=False, labeltop=False, labelbottom=True)
if set_places == True:
    antiplot.axes.get_yaxis()
    antiplot.set_yticks(places)
    antiplot.set_yticklabels(labels, fontdict=None, minor=False)
plotted.colorbar(antiplotax, ticks=[Anti_Diff.min(),0, Anti_Diff.max()],ax=antiplot)

#corrantiax=corranti.matshow(Corr_to_Anti,cmap=cm.Greys_r)
corrantiax=corranti.matshow(Corr_to_Anti,cmap=cm.Greys, origin="lower")
corranti.set_title('Correlated to Anticorrelated')
if set_places == True:
    corranti.axes.get_xaxis()
    corranti.set_xticks(places)
    corranti.set_xticklabels(labels, fontdict=None, minor=False)
corranti.tick_params(axis='x', which='both', bottom=True, top=False, labeltop=False, labelbottom=True)
if set_places == True:
    corranti.axes.get_yaxis()
    corranti.set_yticks(places)
    corranti.set_yticklabels(labels, fontdict=None, minor=False)
plotted.colorbar(corrantiax, ticks=[0, Corr_to_Anti.max()],ax=corranti).ax.set_yticklabels(['-1','+1'])

#anticorrax=anticorr.matshow(Anti_to_Corr,cmap=cm.Greys_r, origin="lower")
anticorrax=anticorr.matshow(Anti_to_Corr,cmap=cm.Greys, origin="lower")
anticorr.set_title('Anticorrelated to Correlated')
if set_places == True:
    anticorr.axes.get_xaxis()
    anticorr.set_xticks(places)
    anticorr.set_xticklabels(labels, fontdict=None, minor=False)
anticorr.tick_params(axis='x', which='both', bottom=True, top=False, labeltop=False, labelbottom=True)
if set_places == True:
    anticorr.axes.get_yaxis()
    anticorr.set_yticks(places)
    anticorr.set_yticklabels(labels, fontdict=None, minor=False)
plotted.colorbar(anticorrax, ticks=[0, Anti_to_Corr.max()],ax=anticorr).ax.set_yticklabels(['-1','+1'])


plotted.tight_layout()
#plotted.savefig(outputfile, dpi=1000,layout="tight")
plotted.savefig(outputfile, dpi=600)


def mc_plot_indiv(data,specific_title,sing_cmap,outfile,set_places):
    """Generate a matrix correlation plot"""
    ## Rotate the data 90 degrees so it plots correctly
    ## Necessary here since we're using a different plotting method
    data = np.rot90(data)
    ## Plot the data
    global test
    sm.graphics.plot_corr(data,normcolor=(-1.0,1.0),cmap=sing_cmap)
    ## Fix the axes
    ax = plt.gca()
    ax.axes.get_xaxis()
    if set_places == True:
        global places, labels
        ax.set_xticks(places)
        ax.set_xticklabels(labels, fontdict=None, minor=False)
    ax.axes.get_yaxis()
    if set_places == True:
        ax.set_yticks(places)
        ax.set_yticklabels(labels, fontdict=None, minor=False)
    ## Do not print generic 'Correlation Matrix' label
    ax.set_title(specific_title)
    ## Default DPI is 100, and high-res for publication is 300 minimum
    plt.savefig(outfile, dpi=600)
    plt.close(outfile)

## Generate the single plots
if save_singles == True:
    ## Set up the individuals
    ## Data, Title on Plot ('' if none wanted), the colormap to use, and the
    ##  output file.
    single_sets = [
    (Corr_Diff, "Correlated Movement", corr_shifted_cmap, "Corr.png"),
    (Anti_Diff, "Anticorrelated Movement", anti_shifted_cmap, "Anti.png"),
    (Corr_to_Anti, "Correlated to Anticorrelated", cm.Greys, "Corr2Anti.png"),
    (Anti_to_Corr, "Anticorrelated to Correlated", cm.Greys, "Anti2Corr.png")
    ]
    for data, specific_title, sing_cmap, outfile in single_sets:
        mc_plot_indiv(data,specific_title,sing_cmap,outfile,set_places)
