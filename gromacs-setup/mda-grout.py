import MDAnalysis as mda

## Load in the Universe
u = mda.Universe("system.gro")

## Remove the segment IDs (aka `SYSTEM`) for prettier AtomGroup printing
for atom in u.atoms:
    atom.segment.segid = ''

## Select specific ions you want
AG = u.select_atoms("not resname SOL")

## Write out a PDB for easy checking
AG.write("system_no_SOL.pdb")

## Can't rewrite topology, but that's ok
AG.write("system_no_SOL.gro", reindex=False)
