import MDAnalysis as mda

frame_pdb = "aligned_structure.pdb"
coord_line_file = "coords.txt"

in_nmd = "WT_protein_system_100.nmd"
out_nmd = "fixed_WT_protein_system_100.nmd"

## Use MDAnalysis selection language to match what residues are currently in
##  your NMD file. Check the `resid` list, since cpptraj may have skipped some!
nmd_atom_selection = "name CA P C4\' C2 and resid 1:301"

## --------------------------- Behind the Curtain ---------------------------#

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
