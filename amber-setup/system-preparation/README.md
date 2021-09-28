# system-preparation
This directory contains scripts for preparing the biomolecular structure for MD.

## `clean-pdb.sh`
This script is meant to clean up a PDB downloaded from [RCSB](https://www.rcsb.org/)
to use in MD.
By default, it removes `REMARK` and `ANISOU` lines, saves only A-chain residues,
and renames `HOH` to `WAT`.
There are also comments for lines that would remove common inhibitors used in
crystallization.

## `propka-integration.py`
This script take the place of the `grep` command for any AMBER-specific naming
after running [PROPKA/PDB2PQR](https://server.poissonboltzmann.org/pdb2pqr) for
use with an AMBER forcefield.

:warning: **Important!** :warning:

Make sure you use the `AMBER` forcefield and `AMBER` output naming scheme when
you submit your job on the PDB2PQR webserver, or use the command-line options
to specify it on a local installation.
```
--ff=AMBER --ffout=AMBER
```

Because AMBER has specific names for different protonation states, those states
need to be incorporated into the PDB structure given to MolProbity/LEaP.
PDB2PQR also attempts to optimize the hydrogen bonding network which affects
the position of some atoms.
This script integrates any non-hydrogen changes made/suggested by PROPKA into
the PDB you had been working from.

> That `grep` command, for posterity:
> ```
> grep - e 'ASH\|CYM\|CYX\|GLH\|HID\|HIP\|LYN' input.pqr > output.pqr
> ```
