#!/bin/bash

## Define variables
## f will name your files, RES is the first residue name, RESB is the second residue name
## g is the file generated through rmagic-EDA-avg-diffs.r
f="SYSTEM-WT--MUT-A"
g="/absolute/path/to/the/difference/output/WT-MUTA_total_interaction_resX_avg.dat"
RESA="WT"
RESB="MUT-A"

## Coul
## Set up the information for Chimera
## Note: attribute cannot start with a capital letter
echo "#TET${RESA}TET${RESB}" > ${f}_EDA_tot_chimera.txt
echo "attribute: tet${RESA}${RESB}tot" >> ${f}_EDA_tot_chimera.txt
echo "match mode: 1-to-1" >> ${f}_EDA_tot_chimera.txt
echo "recipient: residues" >> ${f}_EDA_tot_chimera.txt

## This will skip the header and print the residue number ($1) and the total difference ($2)
awk 'NR > 1 {printf "\t:%-3s\t%-22s\n", $1, $2}' $g >> ${f}_EDA_tot_chimera.txt
