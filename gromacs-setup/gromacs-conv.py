import parmed as pmd

## Load the ligand
sys = pmd.load_file("system.prmtop", "system.inpcrd")

## Save to GROMACS
sys.save('system.top', format='gromacs')
sys.save('system.gro')
