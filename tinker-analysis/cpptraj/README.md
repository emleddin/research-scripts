# cpptraj
These are example scripts for running `cpptraj` on TINKER trajectories.

## `cpptraj_1.sh`
A PBS script for running `cpptraj_strip.in`.

## `cpptraj_strip.in`
Reads in a PDB to understand residue information, all of the trajectory files,
autoimages, and writes out a stripped trajectory in NetCDF format.

## `cpptraj_2.sh`
A PBS script for running `cpptraj_analysis.in`.

## `cpptraj_analysis.in`
Relies on the successful completion of `cpptraj_1.sh`. This reads in the
autoimaged, combined trajectory and analyzes it based on the contained arguments.
