# LICHEM Tools

This directory contains a number of scripts for working with
[LICHEM](https://github.com/CisnerosResearch/LICHEM).

General Tools:
- [`get-energies.py`](https://github.com/emleddin/research-scripts/tree/main/LICHEM-tools#get-energiespy)
- [`qm-pdb-map.py`](https://github.com/emleddin/research-scripts/tree/main/LICHEM-tools#qm-pdb-mappy)

File conversion and visualization:
- [`create-reg.py`](https://github.com/emleddin/research-scripts/tree/main/LICHEM-tools#create-regpy)
- [`vmd-regions.py`](https://github.com/emleddin/research-scripts/tree/main/LICHEM-tools#vmd-regionspy)
- [`lichem-setup.sh`](https://github.com/emleddin/research-scripts/tree/main/LICHEM-tools#lichem-setupsh)

Modifying Structures:
- [`mda-qm-part1.py`](https://github.com/emleddin/research-scripts/tree/main/LICHEM-tools#mda-qm-part1py)
- [`mda-qm-part2.py`](https://github.com/emleddin/research-scripts/tree/main/LICHEM-tools#mda-qm-part2py)
- [`mda-recenter-pdb.py`](https://github.com/emleddin/research-scripts/tree/main/LICHEM-tools#mda-recenter-pdbpy)
- [`mda-recenter-txyz.py`](https://github.com/emleddin/research-scripts/tree/main/LICHEM-tools#mda-recenter-txyzpy)
- [`swapsies.py`](https://github.com/emleddin/research-scripts/tree/main/LICHEM-tools#swapsiespy)
- [`vmd-qm-part1.py`](https://github.com/emleddin/research-scripts/tree/main/LICHEM-tools#vmd-qm-part1py)

QSM:
- [`dissected.py`](https://github.com/emleddin/research-scripts/tree/main/LICHEM-tools#dissectedpy)
- [`stitching.py`](https://github.com/emleddin/research-scripts/tree/main/LICHEM-tools#stitchingpy)
- [`qsm-energy-diagram.py`](https://github.com/emleddin/research-scripts/tree/main/LICHEM-tools#qsm-energy-diagrampy)

## `create-reg.py`
VMD and LICHEM use one numbering system, whereas TINKER uses another.
This script provides a skeleton for building a `regions.inp` to be used
with LICHEM, as well as the `BASIS` file.
A file mapping the `BASIS` file to the numbering in VMD and TINKER, and the
respective assignments will be created (`BASIS_verification.txt`).

The `select_QM` function needs heavy modification to make it apply to your
QM/MM system, because you need to define the unique QM region.
The `select_higher_basis` function also needs to me modified, so that the atoms
needing a higher level of theory can be specified.
The benefit to scripting this process is being able to use the MDAnalysis
selection language.
With atom selection commands, you can explicitly list the atoms in your QM
region using whatever means is most logical to you.

## `dissected.py`
This script does two things.
First, it can create a `BurstStruct.xyz` that uses the reactant MM for some
beads and the product MM for others (aka swapsies), as opposed to the default
reactant MM for all.
Second, it can create a `BeadStartStruct.xyz` from a `BurstStruct.xyz` file.
The output of each of these created files are marked by `new_`.

> Note:
> The formatting for the `new_BurstStruct.xyz` and the `new_BeadStartStruct.xyz`
> May not match. I'm working on a portable MDAnalysis writer modification for
> the output XYZ file, since it does not match LICHEM's and has fewer digits.

## `get-energies.py`
This script will parse the `LICHM_GaussEnergy` and `LICHM_TINKEREnergy` file
to separate out the QM and MM energy given in the LICHEM log.
These energies can be converted to different units, such as eV and kcal/mol.

## `lichem-setup.sh`
This bash script is intended for multi-system preparation, starting from a PDB
structure.
It is written for 2 separate systems within a base folder, and it will:
- Convert the PDB to a Tinker XYZ
- Create the LICHEM `regions.inp` and Gaussian `BASIS` files
- Convert the Tinker XYZ to LICHEM, creating `connect.inp` and `xyzfile.xyz`
- Write a `.vmd` file with information from the `regions.inp` file
- Copy all the necessary files for running a single point calculation to one
  subfolder

The skeleton for `create-reg.py` needs to be modified for your system!
Additionally, both `create-reg.py` and the chosen
[PDBXYZ](https://github.com/emleddin/pdbxyz-xyzpdb)
script need to be modified to allow system arguments.
An explanation of this is in the script itself.

It is *highly* recommended that you test this process with 1 individual frame
of each to debug it before preparing multiple frames!

## `qsm-energy-diagram.py`
This script parses the QSM output file and creates figures of plots for the
initial and final reaction coordinates.

## `stitching.py`
This script creates a complete `BeadStartStruct.xyz` from two individual
`BeadStartStruct.xyz` files.
This way you can build one path between the reactant and intermediate, and
a second path between the intermediate and the product.

> **Example Case**
> ```bash
> $ lichem -path -b 9 -r reactant.xyz -p intermediate.xyz
> $ mv BeadStartStruct.xyz react-int-BSS.xyz
> $ lichem -path -b 9 -r intermediate.xyz -p product.xyz
> $ mv BeadStartStruct.xyz int-prod-BSS.xyz
> $ python3 stitching.py
> ```

## `swapsies.py`
This script is a combination of `mda-qm-part1.py` and `mda-qm-part2.py`.
"Swapsies"<sup>*</sup> is putting the product QM region into the reactant MM prior to
reoptimizing the product, likely because the MM environment migrated a lot
between independent reactant and product optimizations.
The script takes in the regions file and the reactant and product XYZ files.
From there, it creates a PDB of the replacing QM atoms, and puts those in the
other XYZ. Snip snip! :scissors: :smiley:

The commands for "Reverse Swapsies" (reactant QM in product MM) are also
provided as comments.

> *: "Swapsies" was coined by the esteemed [Madison](MadisonB14).

## `vmd-regions.py`
This script parses the `regions.inp` file and creates VMD macros for:
- QM atoms: `qm`, `QM`, `quantum`
- Pseudobond atoms: `pseudobond`, `pb`, `pseudo`,
- Boundary atoms: `boundary`, `bound`
- Frozen: `frozen`, `f`
- Unfrozen atoms: `unfrozen`, `uf`

You can then read in the resulting file with your XYZ:
```
vmd -xyz LICHEM_output.xyz -e vmd-selections.vmd
```
and use those keywords under `Selected Atoms` in the representations menu.

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

### `qm-pdb-map.py`
If you ever need to Frankenstein a QM region together, it will help to know
what atom is which.
Using this script will make a `PDB_verification.txt` file that will contain
the PDB file's index for QM structures obtained through `mda-qm-part1.py`.

### `vmd-qm-part1.py`
An alternative to running `mda-qm-part1.py`, this creates a VMD command file
to save the QM PDB, as well as view the QM region of the output XYZ.
The generated command file can be opened from the command line with:
```
$ vmd -e view-qm.vmd
```
