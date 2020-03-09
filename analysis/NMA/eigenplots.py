import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statsmodels.api as sm
from matplotlib.mlab import griddata
from tables import *
from matplotlib.colors import LinearSegmentedColormap

data1 = np.genfromtxt("WT_protein_system_r1_100.nmd", delimiter=None,skip_header=9)

plt.rcParams.update({'font.size': 22})
plt.rcParams.update({'figure.autolayout': True})

data1eigenrank = data1[:,1]
data1eigenvalue = data1[:,2]
plt.tight_layout()
#x-axis 0 to 10; y-axis 0 to 50
plt.axis([0,10,0,50])
plt.xlabel('Mode Number')
plt.ylabel('Percentage of\nTotal Motion (%)')
plt.plot(data1eigenrank,data1eigenvalue,marker='o',c='black',linewidth=2.0)
plt.savefig('WT_protein_system_eigenplot.png')
plt.gcf().clear()
