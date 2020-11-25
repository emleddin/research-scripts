# LICHEM Tools

This directory contains a number of scripts for working with
[LICHEM](https://github.com/CisnerosResearch/LICHEM).

## `create-reg.py`
VMD and LICHEM use one numbering system, whereas TINKER uses another.
This script provides a skeleton for building a `regions.inp` to be used
with LICHEM, as well as a file mapping the `BASIS` file to the numbering in
VMD and TINKER.
It needs heavy modification to make it apply to your QM/MM system, because you
need to define the unique QM region.
The benefit to scripting this process is being able to use the MDAnalysis
selection language.
With atom selection commands, you can explicitly list the atoms in your QM
region using whatever means is most logical to you.

## Ways to Center a Structure on the Origin
Using `cpptraj` to write out a file without the `center :1-XXX origin mass`
command means that the snapshot will likely not have the protein's center of
mass around the origin as AMBER centers it within the periodic box.
The issue with this is that TINKER expects the structure to be centered around
the origin, and an uncentered system will lead to massive issues upon visual
inspection of the optimized structure, as well as in the overall QM/MM energy.
The options for the regions file are set in the `make_regions` function.

### `mda-recenter-pdb.py`
Reads in a PDB and writes out a recentered PDB.

### `mda-recenter-txyz.py`
Uses a PDB and a TINKER XYZ to write out a recentered TINKER XYZ file and a
traditional XYZ file.

## Saving a PDB of QM atoms and Merging QM Position Changes into an XYZ
GaussView's angle and bond tools are great for modifying QM regions, but
incorporating the changes into a QM/MM system can be tedious.
`mda-qm-part1.py`, `mda-qm-part2.py`, and `vmd-qm-part1.py` can be used to
generate PDBs of the QM atoms from a LICHEM QM/MM optimization.
These can then be edited using your favorite PDB manipulation program, and
finally reincorporated into the original XYZ.

### `mda-qm-part1.py`
This script saves a PDB file from a LICHEM XYZ that consists of only the QM
atoms listed in the LICHEM regions file.

### `mda-qm-part2.py`
After running `mda-qm-part1.py` and modifying the output QM PDB file to reflect
the changes you want to incorporate (i.e., generating a product QM/MM structure
from an optimized reactant), this script will replace the positions of the QM
atoms in the original XYZ with the modified QM atom positions.

> **WARNING**:  MDAnalysis will include frame number by default in
  the XYZ molecule_name line. LICHEM cannot handle any information in this line.
  (There is currently a pull request which will fix this on the MDAnalysis end,
  using remarks, and this script will work properly once it is incorporated.)

### `vmd-qm-part1.py`
An alternative to running `mda-qm-part1.py`, this creates a VMD command file
to save the QM PDB, as well as view the QM region of the output XYZ.
The generated command file can be opened from the command line with:
```
$ vmd -e view-qm.vmd
```
