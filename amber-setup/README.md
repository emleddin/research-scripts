# AMBER Set-Up
This directory contains a collection of scripts for preparing systems and
running simulations with AMBER MD.

## Berendsen
This directory contains scripts that use the Berendsen thermostat.

`mineq.sh` is a script for looping through minimization to equilibration using
CPUs. It calls `mdin.1` through `mdin.10`.
`dyanmicsgpu.sh` uses GPUs for production. It calls `mdin.11`.

You can read more about what's happening in `mdin.1` to `mdin.11` at:
[https://emleddin.github.io/comp-chem-website/AMBERguide-example-files.html#NVT](https://emleddin.github.io/comp-chem-website/AMBERguide-example-files.html#NVT)

## Langevin
This directory contains scripts that use the Langevin thermostat.

### hedi-NPT-NVT
`mineq.sh` is a script for looping through minimization to equilibration using
CPUs. It calls `mdin.1` through `mdin.3`.
`dyanmicsgpu.sh` uses GPUs for production. It calls `mdin.4`.

You can read more about the `mdin` files in the `hedi-NPT-NVT` subfolder at:
[https://emleddin.github.io/comp-chem-website/AMBERguide-example-files.html#NVT-NPT](https://emleddin.github.io/comp-chem-website/AMBERguide-example-files.html#NVT-NPT)

### miller-NPT
`SLURM-jobfile` is a jobfile built for running AMBER on an XSEDE resource.
Each of the `mdin` files can be renamed numerically, if you wish to use a
combination of `../hedi-NPT-NVT/mineq.sh` and `../hedi-NPT-NVT/dynamicsgpu.sh`.

You can read more about the `mdin` files in the `miller-NPT` subfolder at:
[https://emleddin.github.io/comp-chem-website/AMBERguide-example-files.html#NPT](https://emleddin.github.io/comp-chem-website/AMBERguide-example-files.html#NPT)

## submission-scripts
This directory contains scripts for running SLURM jobs on SDSC's Comet or
UNT's Talon.

## `clean-pdb.sh`
This script is meant to clean up a PDB downloaded from [RCSB](https://www.rcsb.org/)
to use in MD.
By default, it removes `REMARK` and `ANISOU` lines, saves only A-chain residues,
and renames `HOH` to `WAT`.
There are also comments for lines that would remove common inhibitors used in
crystallization.
