import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from MDAnalysis.analysis import align, rms, hbonds

#------------------------#
#       Input Files      #
#------------------------#
## Load in a PDB for bonding info to use a parmtop
pdb = 'WT_protein_system_cryst.pdb'

## Load in the trajectory files
file1 = 'WT_protein_system_NVT.arc_1'
file2 = 'WT_protein_system_NVT.arc_2'
file3 = 'WT_protein_system_NVT.arc_3'
file4 = 'WT_protein_system_NVT.arc_4'
file5 = 'WT_protein_system_NVT.arc_5'
file6 = 'WT_protein_system_NVT.arc_6'
file7 = 'WT_protein_system_NVT.arc_7'
file8 = 'WT_protein_system_NVT.arc_8'
file9 = 'WT_protein_system_NVT.arc_9'
file10 = 'WT_protein_system_NVT.arc_10'
file11 = 'WT_protein_system_NVT.arc_11'
file12 = 'WT_protein_system_NVT.arc_12'
file13 = 'WT_protein_system_NVT.arc_13'
file14 = 'WT_protein_system_NVT.arc_14'
file15 = 'WT_protein_system_NVT.arc_15'
file16 = 'WT_protein_system_NVT.arc_16'
file17 = 'WT_protein_system_NVT.arc_17'
file18 = 'WT_protein_system_NVT.arc_18'
file19 = 'WT_protein_system_NVT.arc_19'
file20 = 'WT_protein_system_NVT.arc_20'

#--------------------------------#
#     Outfiles and Run Info      #
#--------------------------------#
pdb_out = 'WT_protein_system_final_frame.pdb'
rmsd_out = 'WT_protein_system_rmsd.dat'
hbond_out = 'WT_protein_system_hbond.dat'
rmsf_out = 'WT_protein_system_rmsf.dat'

rmsd_plot = 'WT_protein_system_RMSD.png'
rmsf_plot = 'WT_protein_system_RMSF.png'

## Select residues/atoms of interest to calculate RMSD, hbonds, & RMSF
ROI = 'resid 1:300'

## Load in the trajectory
system = mda.Universe(pdb, file1, file2, file3, file4, file5, file6, file7, \
file8, file9, file10, file11, file12, file13, file14, file15, file16, file17, \
file18, file19, file20, format='ARC', dt=0.5)

## Reference (Crystal PDB)
ref = mda.Universe(pdb)

#----------------------#
#          PDB         #
#----------------------#
# ## Remove the disgusting segment IDs
# for atom in system.atoms:
#     atom.segment.segid = ''
#
# ## Write out the final frame of a trajectory as PDB
# with mda.Writer(outfile, multiframe=False) as pdb:
#     system.trajectory[-1]
#     pdb.write(system)

#----------------------#
#         RMSD         #
#----------------------#
## Select atoms for RMSD
resnames = system.select_atoms(ROI)
resnames_ref = ref.select_atoms(ROI)

## Get the RMSDs
rmsd = mda.analysis.rms.RMSD(resnames, resnames_ref)
rmsd.run()

## Create a Pandas dataframe with info
rmsd_df =  pd.DataFrame(columns = ['#Frame', 'RMSD'])
for i in range(len(rmsd.rmsd)):
    rmsd_df = rmsd_df.append({
    ## First indice is frame, then info is Frame, Time, RMSD aka [0, 1, 2]
    "#Frame" : rmsd.rmsd[i][0],
    "RMSD" : rmsd.rmsd[i][2]
    }, ignore_index=True)

## Round to 4 digits
rmsd_df = rmsd_df.round(4)

rmsd_df.to_csv(rmsd_out, sep='\t', index=False, encoding='utf8', header=True)

#----------------------#
#         RMSF         #
#----------------------#
## Select group of residues
resnames = system.select_atoms(ROI)
## Get residue weights
res_weights = mda.lib.util.get_weights(resnames, weights='mass')

## Run for individual residues
rmsf = mda.analysis.rms.RMSF(resnames)
rmsf.run()

## Apply weights to the by-atom fluctuations
weight_rmsf = []
for atom in range(len(rmsf.rmsf)):
    ## weight_rmsf = fluctuation * atom_mass / residue_mass
    weight_rmsf.append( (rmsf.rmsf[atom] *
     ( res_weights[atom] / rmsf.atomgroup[atom].residue.mass)) )

## Combine into by-residue value
byres_rmsf = resnames.accumulate(weight_rmsf, compound='residues')

## Create a Pandas dataframe with info
rmsf_df =  pd.DataFrame(columns = ['#Res', 'AtomicFlx'])
for i in range(len(byres_rmsf)):
    rmsf_df = rmsf_df.append({
    "#Res" : resnames.residues[i].resnum,
    "AtomicFlx" : byres_rmsf[i]
    }, ignore_index=True)

## Round to 4 digits
rmsf_df = rmsf_df.round(4)

rmsf_df.to_csv(rmsf_out, sep='\t', index=False, encoding='utf8', header=True)

#----------------------#
#     Set-up HBA       #
#----------------------#
## Test for AMBER
# class HydrogenBondAnalysis_AMBER(HydrogenBondAnalysis):
class HydrogenBondAnalysis_AMBER(mda.analysis.hbonds.HydrogenBondAnalysis):
    # use tuple(set()) here so that one can just copy&paste names from the
    # table; set() takes care for removing duplicates. At the end the
    # DEFAULT_DONORS and DEFAULT_ACCEPTORS should simply be tuples.
    #
    #: default heavy atom names whose hydrogens are treated as *donors*
    #: (see :ref:`Default atom names for hydrogen bonding analysis`);
    #: use the keyword `donors` to add a list of additional donor names.
    DEFAULT_DONORS = {
        'AMBER': tuple(set([
            'N','OG','OG1','OE2','OH','NH1','NH2','NZ','OH','N2','ND1','NE', \
            'NE1','NE2','N2','N3','N4','N6','O5\'', 'OD2', ])),
        'CHARMM27': tuple(set([
            'N', 'OH2', 'OW', 'NE', 'NH1', 'NH2', 'ND2', 'SG', 'NE2', 'ND1', \
            'NZ', 'OG', 'OG1', 'NE1', 'OH'])),
        'GLYCAM06': tuple(set(['N', 'NT', 'N3', 'OH', 'OW'])),
        'other': tuple(set([]))}
    #
    #: default atom names that are treated as hydrogen *acceptors*
    #: (see :ref:`Default atom names for hydrogen bonding analysis`);
    #: use the keyword `acceptors` to add a list of additional acceptor names.
    DEFAULT_ACCEPTORS = {
        'AMBER': tuple(set([
            'O','OD1','OD2','OE1','OE2','N','ND1','NE2','NZ','OG','OG1','O1',\
            'O2','O4','O6','N1','N3','N6','N7','O1P','O2P','O3\'','O4\'','OH'])),
        'CHARMM27': tuple(set([
            'O', 'OC1', 'OC2', 'OH2', 'OW', 'OD1', 'OD2', 'SG', 'OE1', 'OE1', \
            'OE2', 'ND1', 'NE2', 'SD', 'OG', 'OG1', 'OH'])),
        'GLYCAM06': tuple(set(['N', 'NT', 'O', 'O2', 'OH', 'OS', 'OW', 'OY', 'SM'])),
        'other': tuple(set([]))}

#----------------------#
#         HBA          #
#----------------------#
# hbond_df =  pd.DataFrame(columns = ['#Frame', 'RMSD'])
# hbond = mda.analysis.hbonds.HydrogenBondAnalysis(system, distance=3.0, \
#  angle=120.0)

hbond = HydrogenBondAnalysis_AMBER(system, ROI, ROI, distance=3.0, angle=135.0, \
 forcefield='AMBER')
hbond.run()

## Save by frequency
new_hbond_df = pd.DataFrame(hbond.count_by_type())
new_hbond_df = new_hbond_df.sort_values(by=['frequency'], ascending=False)

new_hbond_df.to_csv(hbond_out, sep='\t', index=False, encoding='utf8', header=True)

#----------------------#
#        Plots         #
#----------------------#
plt.plot(rmsd_df['#Frame'], rmsd_df['RMSD'])
plt.xlabel('Time')
plt.ylabel('RMSD ($\AA$)')
plt.savefig(rmsd_plot, dpi=300)

plt.plot(rmsf_df['#Res'], rmsf_df['AtomicFlx'])
plt.xlabel('Residue')
plt.ylabel('RMSF ($\AA$)')
plt.savefig(rmsf_plot, dpi=300)


## Helpful things
## How many frames total?
# u.trajectory.n_frames
# 2999
# u.trajectory.totaltime
# 1499.0                        ## Based on file1, file2...

## Read back in data frame
# rmsd_file = pd.read_csv(rmsd_out, sep='\t')

## Plot Opt 1
# plt.plot(rmsd_file['#Frame'], rmsd_file['RMSD'])
# plt.xlabel('Time')
# plt.ylabel('RMSD ($\AA$)')

## Plot Opt 2
# plt.plotfile(rmsd_out, cols=(0,1), delimiter='\t', \
#  names=('Time','RMSD ($\AA$)'))
# plt.savefig(rmsd_plot, dpi=300)
