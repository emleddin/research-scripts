import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#--------------#
# Reading Data #
#--------------#
## Check your files to verify that it's 8 lines of header!
dwt1 = np.genfromtxt("WT_protein_system_r1_100.nmd", delimiter=None,skip_header=8)
dwt2 = np.genfromtxt("WT_protein_system_r2_100.nmd", delimiter=None,skip_header=8)

#------------------#
# Setting Up Plots #
#------------------#
## Define a list of tuples with (data, outfile)
## This is for each individual file -- plot each replicate to check
ev_datasets = [
  (dwt1, "WT_protein_system_r1_eigenplot.png"),
  (dwt2, "WT_protein_system_r2_eigenplot.png"),
]

##---------------------------------------------------------------------------##
##---------------------------- Behind the Curtain ---------------------------##
##---------------------------------------------------------------------------##
## Set up MPL parameters
# plt.rcParams.update({'font.size': 22})
plt.rcParams["figure.figsize"] = (5,3)
plt.rcParams.update({'figure.autolayout': True})

#--------------------#
# Define How to Plot #
#--------------------#
def ev_plot(data, outfile):
    ## Newer cpptraj does the eigenvalue column differently....
    ## It's saves a scaled value for the square-root of the inverse eigenvalue
    eigenrank = data[:,1]
    sqr_inv_ev = data[:,2]
    ## Undo square
    inv_ev = sqr_inv_ev ** 2
    ## Inverse the inverse eigenvalue
    scaled_eigenvalue = 1/inv_ev
    ## Descale
    eigenvalue = scaled_eigenvalue / np.sum(scaled_eigenvalue)
    ## Create percentage
    eigenvalue_percent = eigenvalue * 100
    ## Debugging
    # print(f"First eigenvalue contribution is {eigenvalue_percent}%")
    #x-axis 0 to 10; y-axis 0 to 100
    plt.axis([0,10,0,100])
    plt.xlabel('Mode Number')
    plt.ylabel('Percentage of\nTotal Motion (%)')
    plt.plot(eigenrank,eigenvalue_percent,marker='o',c='black',linewidth=2.0)
    plt.savefig(outfile, dpi=600)
    plt.close(outfile)
    plt.gcf().clear()

#--------------#
# Create Plots #
#--------------#
for data,outfile in ev_datasets:
    ev_plot(data, outfile)
