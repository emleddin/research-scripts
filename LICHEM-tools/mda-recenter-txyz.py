import MDAnalysis as mda
import MDAnalysis.transformations
import pandas as pd
import numpy as np

## in_PDB is the original PDB, in_txyz is the converted TINKER XYZ
in_pdb = "WT_protein_system_frame_23456.pdb"
in_txyz = "WT_protein_system_frame_23456_convert.xyz"

## Out_xyz is in XYZ format, out_txyz is in TINKER XYZ format
out_xyz = "test_xyzfile.xyz"
out_txyz = "WT_protein_system_frame_23456_convert_recenter.xyz"

## Read in the PDB as topology and TXYZ as "trajectory"
system = mda.Universe(in_pdb, in_txyz, format='TXYZ', in_memory=True)

## Remove the disgusting segment IDs
for atom in system.atoms:
    atom.segment.segid = ''

## Load in the PDB alone to get the dimensions for the system
pdb = mda.Universe(in_pdb)
system.dimensions = pdb.dimensions

## Translate all the atoms to the origin
system.atoms.translate(-system.select_atoms('all').center_of_mass())

## Writes out a regular XYZ, not a TXYZ
system.atoms.write(out_xyz)

## Read in the original Tinker XYZ
df_txyz = pd.read_csv(in_txyz, skiprows=1, header=None, sep=r'\s{2,}', \
 engine='python')
df_txyz.columns = ['Index', 'A_Name', 'X_Coord', 'Y_Coord', 'Z_Coord', 'Type', \
 'Connect1', 'Connect2', 'Connect3', 'Connect4']

## Read in the generated XYZ with COM at origin
xyz_cols = ['A_Name', 'X_Coord', 'Y_Coord', 'Z_Coord']
df_xyz = pd.read_csv(out_xyz, skiprows=2, header=None, sep=r'\s{2,}', \
 names=xyz_cols, engine='python')

## Combine these specific columns to make a Tinker XYZ with new coordinates
combo = pd.concat([df_txyz['Index'], df_txyz['A_Name'], df_xyz['X_Coord'], \
 df_xyz['Y_Coord'], df_xyz['Z_Coord'], df_txyz['Type'], df_txyz['Connect1'], \
 df_txyz['Connect2'], df_txyz['Connect3'], df_txyz['Connect4']], axis=1)

## Make NaNs blank
combo = combo.replace(np.nan, '', regex=True)

## Left-align atom names
combo['A_Name'] = (lambda s: s.str.ljust(s.str.len().max()))(combo['A_Name'])

## Get total atoms for TXYZ
natom = len(combo['Index'])

## Save the new TXYZ dataframe as a string for printing
print_string = combo.to_string(index=False, index_names=False, header=False)

## Write the new Tinker XYZ file with transformed coords
with open(out_txyz, 'w+') as file:
    file.write(str(natom)+"\n")
    file.write(print_string+"\n")
    file.close()
