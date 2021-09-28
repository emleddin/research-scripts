import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_types

## in_pqr: Use the PQR from PDB2PQR (PROPKA)
## in_pdb: Use the PDB you submitted to PDB2PQR (PROPKA)
## out_pdb: The output PDB updated with the PROPKA information
in_pqr="propka-job-output.pqr"
in_pdb="wildtype-propka-input.pdb"
out_pdb="wildtype-propka-integrated.pdb"

## Read in the input structures
ph_list = mda.Universe(in_pqr)
pdb = mda.Universe(in_pdb)

## Guess the elements for ph_list (PQR)
guessed_elements = guess_types(ph_list.atoms.names)
ph_list.add_TopologyAttr('elements', guessed_elements)

## Select any non-hydrogen from the PQR
## Molprobity removes Hydrogens, LEaP will add hydrogens, just...
## don't use the ones from here
ph_noh = ph_list.select_atoms("not element H")

## Create an residue group with the updated residues
## Use `-1` because pdb starts at 0 and ph_noh starts at 1
ph_ag = mda.ResidueGroup(ph_noh.residues.resids-1, pdb)

## Update the positions and resnames
## Positions change to optmize hydrogen bonding network
## Resnames could become ASH/CYM/CYX/GLH/HID/HIP/LYN in AMBER
ph_ag.atoms.positions = ph_noh.atoms.positions
ph_ag.residues.resnames = ph_noh.residues.resnames

## Remove SegIDs (makes it easier to check "ATOM" lines with grep)
for atom in pdb.atoms:
    atom.segment.segid = ''

## Write out the PDB
pdb.atoms.write(out_pdb)
