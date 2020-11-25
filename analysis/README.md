# AMBER analysis
This directory contains different scripts relevant to analyzing AMBER MD
simulations.

## `residuenumbers.sh`

A `bash` script for writing out the residue numbers used by biologists so
that they can be "keyed" to the ones in the AMBER PDB file that starts at 1.
The example here had a gap in the biologist numbering (due to the actual
residues not being crystallized), which is why it is broken up.
