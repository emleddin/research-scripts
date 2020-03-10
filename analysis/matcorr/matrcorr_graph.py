import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statsmodels.api as sm
from tables import *
from matplotlib.colors import LinearSegmentedColormap

## Read in the data files generated with cpptraj
dwt1 = np.genfromtxt("WT_protein_system_r1_corr_mat.dat",delimiter=None)
dwt2 = np.genfromtxt("WT_protein_system_r2_corr_mat.dat",delimiter=None)
dwt3 = np.genfromtxt("WT_protein_system_r3_corr_mat.dat",delimiter=None)
dwt4 = np.genfromtxt("WT_protein_system_r4_corr_mat.dat",delimiter=None)

da1 = np.genfromtxt("MUT_A_system_r1_corr_mat.dat",delimiter=None)
da2 = np.genfromtxt("MUT_A_system_r2_corr_mat.dat",delimiter=None)
da3 = np.genfromtxt("MUT_A_system_r3_corr_mat.dat",delimiter=None)


#-----------------------------------------------------------------------------#
#                               Saving Data                                   #
#-----------------------------------------------------------------------------#
## Subtract each pair you care about forwards and backwards
## Because Python starts at zero, subtract 1 from the matrix SNP position
dwt1_da1 = np.subtract(dwt1,da1)
da1_dwt1 = np.subtract(da1,dwt1)
np.savetxt('WT_protein_system_minus_MUT_A_system_123.txt',dwt1_da1[122],fmt='%1.2f')
np.savetxt('WT_protein_system_minus_MUT_A_system_225.txt',dwt1_da1[224],fmt='%1.2f')


#-----------------------------------------------------------------------------#
#                               Self-Plots                                    #
#-----------------------------------------------------------------------------#

## Uncomment placesx2, placesy2, labelsx2, labelsy2 to explicitly define
## axis labels (e.g., to match real biological numbering)
# #Explicity choose where to put x and y ticks
# placesx2 = [0, 100, 200, 300, 333, 347, 400, 430]
# placesy2 = [25, 55, 108, 122, 155, 255, 355, 455]
# ## Note: we're not using the inverted y axis
# ## so therefore, this starts at bottom left
#
# #Define those very x and y tick labels
# labelsx2 = [1130, 1230, 1330, 1430, ' ', 1842, ' ', 'DNA']
# labelsy2 = ['DNA', 1895, 1842, 1463, 1430, 1330, 1230, 1130]

def mc_plot(data,outfile):
    """Generate a matrix correlation plot"""
    # global placesx2, placesy2, labelsx2, labelsy2
    sm.graphics.plot_corr(data,normcolor=(-1.0,1.0),cmap='RdYlBu')
    ax = plt.gca()
    ax.axes.get_xaxis()
    # ax.set_xticks(placesx2)
    # ax.set_xticklabels(labelsx2, fontdict=None, minor=False)
    ax.axes.get_yaxis()
    # ax.set_yticks(placesy2)
    # ax.set_yticklabels(labelsy2, fontdict=None, minor=False)
    ax.set_title('')
    plt.savefig(outfile)
    plt.close(outfile)


## Define a list of tuples with (data, outfile)
## This is for each individual file -- plot each replicate to check
self_datasets = [
  (dwt1, "WT_protein_system-1_mc.png"),
  (dwt2, "WT_protein_system-2_mc.png"),
  (dwt3, "WT_protein_system-3_mc.png"),
  (dwt4, "WT_protein_system-4_mc.png"),
  (da1, "MUT_A_system-1_mc.png"),
  (da2, "MUT_A_system-2_mc.png"),
  (da3, "MUT_A_system-3_mc.png"),
]

for data,outfile in self_datasets:
    mc_plot(data,outfile)

#-----------------------------------------------------------------------------#
#                              Cross Plots                                    #
#-----------------------------------------------------------------------------#

## Define a list of tuples with (data, outfile)
## These are for the differences between systems -- plot both A-B & B-A
cross_datasets = [
  (dwt1_da1, "WT-protein-system_minus_MUT-A-system.png"),
  (da1_dwt1, "MUT-A-system_minus_WT-protein-system.png"),
]

for data,outfile in cross_datasets:
    mc_plot(data,outfile)
