# AMBER analysis

This directory contains different scripts relevant to analyzing AMBER MD
simulations.

## EDA
These are a number of scripts for performing an energy decomposition analysis
for AMBER MD trajectories.
They use the `Residue_E_Decomp_07_15.f90` program.
An OpenMP version of this program can be found on the
[CisnerosRes GitHub](https://github.com/CisnerosResearch/AMBER-EDA).

## HBA
These are a collection of scripts for processing hydrogen bond data.

## NMA
These scripts pertain to normal mode analysis (related to principle component
analysis) of protein complexes.

## RMS
These scripts pertain to root mean square analyses (such as deviation and
fluctuation) of protein complexes.

## backbone
These scripts pertain to analyzing a nucleic acid backbone.

## cpptraj
Example scripts for running `cpptraj`.

## matcorr
These scripts pertain to matrix correlation analysis (instances of correlation,
anti-correlation, and no correlation between residues).

## network-analysis
Contains a script for processing cluster-based data from network analysis.

## `residuenumbers.sh`
A `bash` script for writing out the residue numbers used by biologists so
that they can be "keyed" to the ones in the AMBER PDB file that starts at 1.
The example here had a gap in the biologist numbering (due to the actual
residues not being crystallized), which is why it is broken up.
