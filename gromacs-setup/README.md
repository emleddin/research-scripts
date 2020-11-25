# GROMACS Set-Up

This directory contains scripts to help prepare and run GROMACS simulations
using AMBER force fields.

## `add-posres.sh`
Example bash and sed commands for adding the position restraints section to a
`.top` file based on the line it needs to follow.
(This is probably only helpful for batch preparation.)

## `gromacs-conv.py`
Convert from an AMBER `.prmtop` and `.inpcrd` to a GROMACS `.top` and `.gro`.

## `mda-grout.py`
Read in a `.gro` file and use `AtomGroup` language to select what to write out
as both a `.pdb` and `.gro`.

##  `solv-sedit.sh`
Example bash and sed commands for adding specific molecule types to a `.top`
file, or using global string replacements for `WAT` or `SOL`.
