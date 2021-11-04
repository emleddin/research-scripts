# cpptraj
These are example scripts for running `cpptraj`.

## `cpptraj_1.sh`
A PBS script for running `cpptraj_strip.in`.

## `cpptraj_strip.in`
Reads in all of the trajectory files, autoimages, strips the water and ions,
and writes out a stripped trajectory in NetCDF format.

## `cpptraj_2.sh`
A PBS script for running `cpptraj_analysis.in`.

## `cpptraj_analysis.in`
Relies on the successful completion of `cpptraj_1.sh`. This reads in the
stripped trajectory and analyzes it based on the contained arguments.

## `cpptraj-EDA.sh`
A PBS script for running `cpptraj-EDA-nas.in`.

## `cpptraj-EDA-nas.in`
Reads in all of the trajectory files and writes out a complete trajectory in
ASCII format. This can be submitted at the same time as `cpptraj_1.sh`.

## `write-cpptraj.py`
Creates the `cpptraj_` bash and input files for stripping, analysis, and
EDA.
Modify the `write_analy_traj` function to add more analyses.
This script is particularly helpful for batch creation.
