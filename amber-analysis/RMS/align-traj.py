import MDAnalysis as mda
from MDAnalysis.analysis import align
​
oriented_pdb = "WT-protein-system-image-orientation.pdb"
​
## Option: read in PDB or prmtop/rst combo
## PDB cannot have asterisk in output
# frame_pdb = "WT-protein-system-final-frame.pdb"
frame_prm = "WT-protein-system.prmtop"
## RST should be rewritten using cpptraj as ASCII RST
## Ex: trajout final-frame.rst restart
frame_rst = "final-frame.rst"
​
## Set up range of waters to keep (5 angstrom around residue)
wat_range = "(resname WAT) and (around 5 resid 123)"
​
## Select resid range you want in final PDB
res_range = "resid 1:500"
​
## Name of output PDB
out_pdb = "test.pdb"
​
## Read in the reference and final frame structures
ref = mda.Universe(oriented_pdb)
# fin_frame = mda.Universe(frame_pdb)
fin_frame = mda.Universe(frame_prm, frame_rst, format="RESTRT")
​
## Align the final structure to the reference
align.alignto(fin_frame, ref, select="name CA", match_atoms=True)
​
## Select any water residues within given range + the residue range
struct_wat = fin_frame.select_atoms(wat_range + " or " + res_range)
​
## Write the output PDB
struct_wat.write(out_pdb)
