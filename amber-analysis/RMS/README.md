# RMS

These scripts pertain to root mean square analyses (such as deviation and
fluctuation) of protein complexes.
Total number of hydrogen bonds are also included.

## `align-traj.py`
A `Python3` script will align a structure to a reference and can be used to
select out certain atoms.
As written, it will select all water residues within a specific distance
sphere of a specific residue.

## `avg-rmsd-etc.gnu`
A `gnuplot` script for graphing the RMSD, RMSF, and average hydrogen bond data
and standard deviations generated using `rmagic-rmsd-rmsf-hbond-5.r`.

## `byres-processing.R`
Using `cpptraj` for by-residue data results in a column for each residue with
all the frames as rows.
This script will process that data file (regardless of the number of columns)
and write out a 2-column data set for `Residue` and `Average`.
Residues will use the numbering in the PDB.

## `r-gg-rmsd-rmsf.gnu`
An `R` script for plotting RMSD and RMSF for multiple systems in individual and
composed figures.

## `rmsd-etc.gnu`
A `gnuplot` script for plotting RMSD, RMSF, and average hydrogen bonds present
on separate graphs.
The `u ($1/500):($2)` divides the Frame number column by an appropriate number
(in this case, 500) to get nanoseconds.
It likely will not be `500`, depending on what you used for you mdin files.

## `rmagic-rmsd-rmsf-hbond-5.r`
This script will average RMSDs, RMSFs, and average hydrogen bonds for
replicate trajectories over the same timescale and give standard deviations.
The output files can be plotted using `avg-rmsd-etc.gnu`.

## `moving-averages.py`
A `Python3` script for plotting RMSD, RMSF, and average hydrogen bonds present
on separate graphs.
RMSD and average hydrogen bonds present are plotted as moving averages based
on the previous 1 ns; RMSD is plotted as given.
The `frame2ns` variable converts Frame number column to time in nanoseconds
using a multiplier (in this case, 1/500).
It likely will not be `500`, depending on what you used for you mdin files.

## `natoms-in-range.py`
A `Python3` script for plotting the number of a specific type of ion within
a range of a residue or protein throughout the course of a simulation.
This is not intended to determine which ion, specifically, is that location,
just total the type.

## Rendered Examples
This directory has examples of the Python, Gnuplot, and R figures made from the
above plotting scripts.
