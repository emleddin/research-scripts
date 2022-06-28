# LICHEM Tools

This directory contains a number of scripts for preparing job directories
for [LICHEM](https://github.com/CisnerosResearch/LICHEM).

## General

### `lichem-setup.sh`
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

### `cont-dfp.sh`
This script will copy files from a partially optimized structure to
continue the optimization.
An example use is the job not finishing at the wallclock time.

## For QSM

### `prep-qsm.sh`
This script will create the path for QSM from an optimized reactant and
product.
It will *not* modify the `regions.inp` file for QSM.
If you need to rerun this script, be sure to remove all the existing XYZ files
in the QSM directory.

### `cont-qsm.sh`
This script will copy files from an optimized QSM path to continue the
optimization.
An example use is running unrestrained QSM after restrained QSM.
It will *not* modify the `regions.inp` file for QSM.
You can rerun this script without needing to delete files.
