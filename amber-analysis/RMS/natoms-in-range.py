import MDAnalysis as mda
import pandas as pd
import matplotlib.pyplot as plt
from cycler import cycler

## Topology/trajectory combo for universe
top = 'trajectory.prmtop'
traj = 'trajectory.nc'

## Residues in group to find ions within X dist around (protein/nucleic/drug)
residues_in = 'resid 1:1500'

## Ion Groups
MG_in = 'resid 1501:1510'
K_in = 'resid 1511:1520'
CL_in = 'resid 1521:1540'

## What's the farthest away from the residues the ions can be? (Angstrom)
ion_dist = 7

## Figure name
save_name = "ion_distances.png"

## Multiply timesteps by this fraction to convert frames to nanoseconds
frame2ns = 1/1000

#---------- No need to modify past the curtain ---------------#
## Read in files
## NetCDF may need to specify "format=NCDF" and load netCDF4
system = mda.Universe(top, traj)

## Set up Ion Groups
residues = system.select_atoms(residues_in)
MG = system.select_atoms(MG_in)
K = system.select_atoms(K_in)
CL = system.select_atoms(CL_in)

## Must have space for adding string!
# close_ions = system.select_atoms("around 7 resid 1:24")
# close_ions = system.select_atoms("around 7 " + residues_in)
## Get how many are within 7 A of system over simulation
#traj_ions = system.select_atoms("around " str(ion_dist) + " " +
# residues_in, updating=True)

## Get the groups of atoms matching critera per frame of trajectory
## Ex: "(around 7 resid 1:50) and (name MG)"
## The parentheses are very important!!!!!
traj_MG = system.select_atoms("(around " + str(ion_dist) + " " +
 residues_in + ") and (name MG)", updating=True)
traj_K = system.select_atoms("(around " + str(ion_dist) + " " +
 residues_in + ") and (name K)", updating=True)
traj_CL = system.select_atoms("(around " + str(ion_dist) + " " +
 residues_in + ") and (name CL)", updating=True)

## Initialize list
d = []
for ts in system.trajectory:
    d.append((ts.time, len(traj_MG), len(traj_K), len(traj_CL)))

## Create a DataFrame
ion_df = pd.DataFrame(d, columns=('Timestep', 'MG', 'K', 'CL'))

## Set Up Graphs
## Set up different color axes based on https://personal.sron.nl/~pault/
plt.rcParams.update({'axes.prop_cycle': cycler('color', ['#332288', '#88CCEE',
'#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499'])})

## Plot the Data
plt.plot(ion_df["Timestep"]*frame2ns, ion_df["MG"], label="$Mg^{2+}$")
plt.plot(ion_df["Timestep"]*frame2ns, ion_df["K"], label="$K^+$")
plt.plot(ion_df["Timestep"]*frame2ns, ion_df["CL"], label="$Cl^-$")

plt.xlabel('Time (ns)', labelpad=5)
plt.ylabel('Number of Ions', labelpad=10)
plt.legend(loc='best')

plt.savefig(save_name, dpi=300)
plt.close(save_name)
