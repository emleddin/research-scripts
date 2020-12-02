#!/bin/bash

## Define variables
## f will name your files, RESA is the first residue name, RESB is the
## second residue name, and INFILE is the file created when making the
## correlation plots when using the Python script
f="SYSTEM"
RESA="WT"
RESB="MUTA"
INFILE='/path/to/file/created/by/matr_corr.py'

## Find number of lines in generated matrix value (i.e. number of residues)
## And paste those numbers and data together into one file
a=$(wc $INFILE)
LC=$(echo $a|cut -d' ' -f1)
seq 1 $LC > ${f}-numbers.txt
paste ${f}-numbers.txt $INFILE > ${f}-pastenumbers.txt

## Set up information for Chimera attribues
## Note: attribute cannot start with a capital letter
echo "#SYSTEM${RESA}SYSTEM${RESB}" > ${f}_mc_chimera.txt
echo "attribute: system${RESA}${RESB}mc" >> ${f}_mc_chimera.txt
echo "match mode: 1-to-1" >> ${f}_mc_chimera.txt
echo "recipient: residues" >> ${f}_mc_chimera.txt

## This will print the formatted residue number and value
awk '{printf "\t:%-3s\t%-5s\n", $1, $2}' ${f}-pastenumbers.txt >> ${f}_mc_chimera.txt

## Remove the lists of numbers generated with paste
rm ${f}-numbers.txt
rm ${f}-pastenumbers.txt
