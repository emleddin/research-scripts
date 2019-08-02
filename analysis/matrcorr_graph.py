import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statsmodels.api as sm
from tables import *
from matplotlib.colors import LinearSegmentedColormap

data1 = np.genfromtxt("WT-system/WT_protein_system_matrix_correl.dat",delimiter=None)
data2 = np.genfromtxt("MUT-A-system/MUT-A-system_matrix_correl.dat",delimiter=None)
data3 = np.genfromtxt("MUT-B-system/MUT-B-system_matrix_correl.dat",delimiter=None)
data4 = np.genfromtxt("MUT-C-system/MUT-C-system_matrix_correl.dat",delimiter=None)

## Saving Data
data12 = np.subtract(data1,data2)
data21 = np.subtract(data2,data1)
np.savetxt(’WT_minus_MUT-A_436.txt’,data12[435],fmt=’%1.2f’)

data13 = np.subtract(data1,data3)
data31 = np.subtract(data3,data1)
np.savetxt(’MUT-A_minus_MUT-B_436.txt’,data13[435],fmt=’%1.2f’)

data14 = np.subtract(data1,data4)
data41 = np.subtract(data4,data1)
np.savetxt(’WT_minus_MUT-C_436.txt’,data14[435],fmt=’%1.2f’)

data23 = np.subtract(data2,data3)
data32 = np.subtract(data3,data2)
np.savetxt(’MUT-A_minus_MUT-B_436.txt’,data23[435],fmt=’%1.2f’)

data24 = np.subtract(data2,data4)
data42 = np.subtract(data4,data2)
np.savetxt(’MUT-A_minus_MUT-D_436.txt’,data24[435],fmt=’%1.2f’)

data34 = np.subtract(data3,data4)
data43 = np.subtract(data4,data3)
np.savetxt(’MUT-B_minus_MUT-C_436.txt’,data34[435],fmt=’%1.2f’)

## Self Plots
sm.graphics.plot_corr(data1,normcolor=(-1.0,1.0),cmap=’RdYlBu’)
ax = plt.gca()
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
plt.savefig(’WT_protein_system_mc.png’)
plt.close(’WT_protein_system_mc.png’)

sm.graphics.plot_corr(data2,normcolor=(-1.0,1.0),cmap=’RdYlBu’)
ax = plt.gca()
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
plt.savefig(’MUT-A-system_mc.png’)
plt.close(’MUT-A-system_mc.png’)

sm.graphics.plot_corr(data3,normcolor=(-1.0,1.0),cmap=’RdYlBu’)
ax = plt.gca()
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
plt.savefig(’MUT-B-system_mc.png’)
plt.close(’MUT-B-system_mc.png’)

sm.graphics.plot_corr(data4,normcolor=(-1.0,1.0),cmap=’RdYlBu’)
ax = plt.gca()
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
plt.savefig(’MUT-C-system_mc.png’)
plt.close(’MUT-C-system_mc.png’)

## Actual Cross Plots
sm.graphics.plot_corr(data12,normcolor=(-1.0,1.0),cmap=’RdYlBu’)
ax = plt.gca()
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
plt.savefig(’WT_minus_MUT-A_436.png’)
plt.close(’WT_minus_MUT-A_436.png’)

sm.graphics.plot_corr(data21,normcolor=(-1.0,1.0),cmap=’RdYlBu’)
ax = plt.gca()
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
plt.savefig(’MUT-A_minus_WT_436.png’)
plt.close(’MUT-A_minus_WT_436.png’)

sm.graphics.plot_corr(data13,normcolor=(-1.0,1.0),cmap=’RdYlBu’)
ax = plt.gca()
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
plt.savefig(’WT_minus_MUT-B_436.png’)
plt.close(’WT_minus_MUT-B_436.png’)

sm.graphics.plot_corr(data31,normcolor=(-1.0,1.0),cmap=’RdYlBu’)
ax = plt.gca()
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
plt.savefig(’MUT-B_minus_WT_436.png’)
plt.close(’MUT-B_minus_WT_436.png’)

sm.graphics.plot_corr(data24,normcolor=(-1.0,1.0),cmap=’RdYlBu’)
ax = plt.gca()
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
plt.savefig(’MUT-A_minus_MUT-C_436.png’)
plt.close(’MUT-A_minus_MUT-C_436.png’)
sm.graphics.plot_corr(data42,normcolor=(-1.0,1.0),cmap=’RdYlBu’)
ax = plt.gca()
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
plt.savefig(’MUT-C_minus_MUT-A_436.png’)
plt.close(’MUT-C_minus_MUT-A_436.png’)

sm.graphics.plot_corr(data34,normcolor=(-1.0,1.0),cmap=’RdYlBu’)
ax = plt.gca()
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
plt.savefig(’MUT-B_minus_MUT-C_436.png’)
plt.close(’MUT-B_minus_MUT-C_436.png’)

sm.graphics.plot_corr(data43,normcolor=(-1.0,1.0),cmap=’RdYlBu’)
ax = plt.gca()
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)
plt.savefig(’MUT-C_minus_MUT-B_436.png’)
plt.close(’MUT-C_minus_MUT-B_436.png’)
