# TINKER Analysis

## `log_total_energy.py`
This script reads through the logfile from TINKER to graph potential and kinetic
energy by frame.

## `tinker-cuda.sh`
This is a bash script for doing MD with TINKER.
Instead of one massive ARC file, and individual XYZ cycle files for each frame,
it saves ARC files for each specified time chunk.

If you restart this script, it is CRITICAL that you comment out the first XYZ
copy outside of the loop.
Otherwise, everything it will restart.
