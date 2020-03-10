# Backbone Analysis
These scripts pertain to analyzing a nucleic acid backbone.

## `dihedral-avg.sh`
This script averages the per-residue alpha, beta, gamma, delta, epsilon, zeta,
and chi values of a backbone, then comes up with an overall average.
`f` sets the residue number.
The DX5/RX5 and DX3/RX3 residues will be missing specific angles (alpha and
beta for DX5/RX5, epsilon and zeta for DX3/RX5), which are why they are
separated out.

## `grep-nastruct.sh`
Using `nastruct` with `cpptraj` includes a lot of blank lines.
This script removes those blank lines, and renames the data file.

## `BPstep-awk.sh`
Data from `nastruct` with `cpptraj` comes out as base pair matches.
This script provides the commands to resave specific base pairs.
`NR == 1` selects the header information, and `NR % 4 == 0` chooses every 4th
line thereafter.
Likewise, `NR % 4 == 2` selects every 4th line starting from the second.

## `delta-chi.gnu`
A gnuplot script to plot delta and chi backbone angles.

## `shift-slide.gnu`
A gnuplot script to plot shift and slide data.
