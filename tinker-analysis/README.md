# TINKER Analysis

## cpptraj
Contains information on using the `cpptraj` analysis package for analyzing
trajectory information saved in TINKER ARC or TINKER XYZ files.

## `log_total_energy.py`
This script reads through the logfile from TINKER to graph potential and kinetic
energy by frame.

## `mda-rmsd-rmsf-hbond.py`
This script uses the MDAnalysis package to analyze trajectory information saved
in TINKER ARC or TINKER XYZ files.
It is substantially slower than using cpptraj, but the RMSD/RMSF information
will match it.

## `tinker-cuda.sh`
This is a bash script for doing MD with TINKER.
Instead of one massive ARC file, and individual XYZ cycle files for each frame,
it saves ARC files for each specified time chunk.

If you restart this script, it is CRITICAL that you comment out the first XYZ
copy outside of the loop.
Otherwise, everything it will restart.
