import MDAnalysis as mda

pdb = "WT_protein_system_frame_23456.pdb"
out_pdb = "WT_protein_system_frame_23456_centered.pdb"

## Read in the PDB file
system = mda.Universe(pdb)

## Translate all the atoms to the origin
new = system.atoms.translate(-system.select_atoms('all').center_of_mass())

## Writes out the PDB with COM at origin
new.atoms.write(out_pdb)
