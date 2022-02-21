#!/usr/env/python
import MDAnalysis as mda

## Define a list of tuples with (Aligned_PDB, Temporary_Coordinate_File,
##  Input_NMD, Fixed_NMD) as `datasets`.
datasets = [
  ("WT_aligned_structure.pdb", "WT_coords.txt", "WT_protein_system_100.nmd", "fixed_WT_protein_system_100.nmd"),
  ("MUTA_aligned_structure.pdb", "MUTA_coords.txt", "MUTA_protein_system_100.nmd", "fixed_MUTA_protein_system_100.nmd"),
  ("MUTB_aligned_structure.pdb", "MUTB_coords.txt", "MUTB_protein_system_100.nmd", "fixed_MUTB_protein_system_100.nmd"),
  ("MUTD_aligned_structure.pdb", "MUTD_coords.txt", "MUTD_protein_system_100.nmd", "fixed_MUTD_protein_system_100.nmd"),
]

## This assumes all systems use the name selections!!!
## If you need them to be different, use the non-multi version of this script.
## Use MDAnalysis selection language to match what residues are currently in
##  your NMD file. Check the `resid` list, since cpptraj may have skipped some!
nmd_atom_selection = "name CA P C4\' C2 and resid 1:301"

## --------------------------- Behind the Curtain ---------------------------#

def fix_NMD(frame_pdb, coord_line_file, in_nmd, out_nmd):
    ## Read in the PDB
    fin_frame = mda.Universe(frame_pdb)
    ## Pull the selected atoms from the structure
    nmd_atoms = fin_frame.select_atoms(nmd_atom_selection)
    ## Create a file that just contains the new coordinates
    with open(coord_line_file, "w+") as f:
        f.write("coordinates")
        for atom in nmd_atoms.atoms.positions:
            f.write(f" {atom[0]:>8.3f} {atom[1]:>8.3f} {atom[2]:>8.3f}")
        f.write(" \n")

    # Open that coordinates file
    fcoord = open(coord_line_file, "r")
    coord_line = fcoord.readlines()
    fcoord.close()

    ## Open the original NMD file and replace the line
    orig_nmd = open(in_nmd, "r")
    nmd_lines = orig_nmd.readlines()
    for i,line in enumerate(nmd_lines):
        if line.lower().startswith('coordinates'):
            nmd_lines[i] =  coord_line[0]
    orig_nmd.close()

    ## Write out the fixed NMD file
    new_nmd = open(out_nmd, "w+")
    new_nmd.writelines(nmd_lines)
    new_nmd.close()

## Run the script
for frame_pdb, coord_line_file, in_nmd, out_nmd in datasets:
    fix_NMD(frame_pdb, coord_line_file, in_nmd, out_nmd)
