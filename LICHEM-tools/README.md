# LICHEM Tools

GaussView's angle and bond tools are great for modifying QM regions, but
incorporating the changes into a QM/MM system can be tedious.
`mda-qm-part1.py`, `mda-qm-part2.py`, and `vmd-qm-part1.py` can be used to
generate PDBs of the QM atoms from a LICHEM QM/MM optimization.
These can then be edited using your favorite PDB manipulation program, and
finally reincorporated into the original XYZ.

## `mda-qm-part1.py`
This script saves a PDB file from a LICHEM XYZ that consists of only the QM
atoms listed in the LICHEM regions file.

## `mda-qm-part2.py`
After running `mda-qm-part1.py` and modifying the output QM PDB file to reflect
the changes you want to incorporate (i.e., generating a product QM/MM structure
from an optimized reactant), this script will replace the positions of the QM
atoms in the original XYZ with the modified QM atom positions.

<div class="alert alert-danger">
  <strong>Warning!</strong> MDAnalysis will include frame number by default in
  the XYZ molecule_name line. LICHEM cannot handle any information in this line.
  (There is currently a pull request which will fix this on the MDAnalysis end,
  using remarks, and this script will work properly once it is incorporated.)
</div>

## `vmd-qm-part1.py`
An alternative to running `mda-qm-part1.py`, this creates a VMD command file
to save the QM PDB, as well as view the QM region of the output XYZ.
The generated command file can be opened from the command line with:
```
$ vmd -e view-qm.vmd
```
