import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

#--------------#
# Reading Data #
#--------------#
## Read in the data files generated with cpptraj
dwt1 = np.genfromtxt("WT_protein_system_r1_corr_mat.dat",delimiter=None)
dwt2 = np.genfromtxt("WT_protein_system_r2_corr_mat.dat",delimiter=None)
dwt3 = np.genfromtxt("WT_protein_system_r3_corr_mat.dat",delimiter=None)
dwt4 = np.genfromtxt("WT_protein_system_r4_corr_mat.dat",delimiter=None)

da1 = np.genfromtxt("MUT_A_system_r1_corr_mat.dat",delimiter=None)
da2 = np.genfromtxt("MUT_A_system_r2_corr_mat.dat",delimiter=None)
da3 = np.genfromtxt("MUT_A_system_r3_corr_mat.dat",delimiter=None)

#-------------#
# Saving Data #
#-------------#
## Subtract each pair you care about forwards and backwards
## Because Python starts at zero, subtract 1 from the matrix SNP position
dwt1_da1 = np.subtract(dwt1,da1)
da1_dwt1 = np.subtract(da1,dwt1)
np.savetxt('WT_protein_system_minus_MUT_A_system_123.txt',dwt1_da1[122],fmt='%1.2f')
np.savetxt('WT_protein_system_minus_MUT_A_system_225.txt',dwt1_da1[224],fmt='%1.2f')

#------------------#
# Setting Up Plots #
#------------------#
## Define a list of tuples with (data, outfile) for the self-plots
##  as `self_datasets`. Because it's for each individual file, plot
##  each replicate to check.
self_datasets = [
  (dwt1, "WT_protein_system-1_mc.png"),
  (dwt2, "WT_protein_system-2_mc.png"),
  (dwt3, "WT_protein_system-3_mc.png"),
  (dwt4, "WT_protein_system-4_mc.png"),
  (da1, "MUT_A_system-1_mc.png"),
  (da2, "MUT_A_system-2_mc.png"),
  (da3, "MUT_A_system-3_mc.png"),
]

## Do you want to do cross plots as well? (ex: A-B and B-A)
cross_plots = False

## If you want cross plots, define them here as `cross_datasets`
## Define a list of tuples with (data, outfile)
## These are for the differences between systems -- plot both A-B & B-A
# cross_datasets = [
#   (dwt1_da1, "WT-protein-system_minus_MUT-A-system.png"),
#   (da1_dwt1, "MUT-A-system_minus_WT-protein-system.png"),
# ]

## If you want to use a custom ColorMap -- defined in the function NewCmap below
custom_cmap = True

## Set these RGB values for the Bottom, Middle, and Top colors you want
##  in the diverging colormap.
Bot_RGB = [ 94,  60, 153] # Dark Purple
Mid_RGB = [247, 247, 247] # White
Top_RGB = [230,  97,   1] # Dark Orange

## Do you want to explicilty set the axis labels?
set_places = True

## Explicitly define axis ticks/labels (e.g., to match biological numbering)
## Explicity choose where to put x and y ticks
places = [0, 100, 200, 300, 333, 347, 400, 430]

## Define those very x and y tick labels
labels = [1130, 1230, 1330, 1430, ' ', 1842, ' ', 'DNA']

##---------------------------------------------------------------------------##
##---------------------------- Behind the Curtain ---------------------------##
##---------------------------------------------------------------------------##

#-----------------#
# Set Up Colormap #
#-----------------#
def NewCmap(Bot_RGB, Mid_RGB, Top_RGB):
    """Builds a new diverging colormap.
    Parameters
    ----------
    Bot_RGB : list (3 items)
        The RGB values for the bottom-most color for the colorbar.
        Black would be [0,0,0].
    Mid_RGB : list (3 items)
        The middle color, where things diverge.
    Top_RGB : list (3 items)
        The top-most color for the colorbar.
    Returns
    -------
    newcmap : matplotlib.colors.ListedColormap
        A diverging colormap built using the input colors.
    """
    ## Set up RGB colorspace
    N = 256
    ## Build the Bottom to Mid Gradient
    bvals = np.ones((N, 4))
    bvals[:, 0] = np.linspace(Bot_RGB[0]/N, Mid_RGB[0]/N, N)
    bvals[:, 1] = np.linspace(Bot_RGB[1]/N, Mid_RGB[1]/N, N)
    bvals[:, 2] = np.linspace(Bot_RGB[2]/N, Mid_RGB[2]/N, N)
    bottom_color = ListedColormap(bvals)
    #
    ## Build the Mid to Top Gradient
    tvals = np.ones((N, 4))
    tvals[:, 0] = np.linspace(Mid_RGB[0]/N, Top_RGB[0]/N, N)
    tvals[:, 1] = np.linspace(Mid_RGB[1]/N, Top_RGB[1]/N, N)
    tvals[:, 2] = np.linspace(Mid_RGB[2]/N, Top_RGB[2]/N, N)
    top_color = ListedColormap(tvals)
    #
    ## Save the Bottom - Top gradient
    newcolors = np.vstack((bottom_color(np.linspace(0,1,128)),
                              top_color(np.linspace(0,1,128))))
    newcmap = ListedColormap(newcolors, name="MyCmap")
    plt.register_cmap(cmap=newcmap)
    # How to set-up reverse, if wanted
    # revnewcolors = np.vstack((top_color(np.linspace(1,0,128)),
    #                           bottom_color(np.linspace(1,0,128))))
    # rev_newcmap = ListedColormap(revnewcolors, name="MyRevCmap")
    # plt.register_cmap(cmap=rev_newcmap)
    return newcmap

if custom_cmap == True:
    newcmap = NewCmap(Bot_RGB, Mid_RGB, Top_RGB)
    ## Reverse, if wanted
    # newcmap, rev_newcmap = NewCmap(Bot_RGB, Mid_RGB, Top_RGB)

#--------------------#
# Define How to Plot #
#--------------------#
def mc_plot(data,outfile,set_places,custom_cmap):
    """Generate a matrix correlation plot"""
    ## Rotate the data 90 degrees so it plots correctly
    data = np.rot90(data)
    ## Plot the data
    if custom_cmap == True:
        global newcmap
        sm.graphics.plot_corr(data,normcolor=(-1.0,1.0),cmap='MyCmap')
    else:
        sm.graphics.plot_corr(data,normcolor=(-1.0,1.0),cmap='RdYlBu')
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
    ax.set_title('')
    ## Default DPI is 100, and high-res for publication is 300 minimum
    plt.savefig(outfile, dpi=600)
    plt.close(outfile)

#-------------------#
# Create Self-Plots #
#-------------------#
for data,outfile in self_datasets:
    mc_plot(data,outfile,set_places,custom_cmap)

#--------------------#
# Create Cross-Plots #
#--------------------#
if cross_plots == True:
    for data,outfile in cross_datasets:
        mc_plot(data,outfile,set_places,custom_cmap)
