#!/bin/bash

## Residue 1 is ACE cap for linker
echo "" > residuenumbers.txt

## PDB 2 to 444
seq 1132 1463 >> residuenumbers.txt

## PDB 334 to 346
## Linker from 1463 to 1842 in protein

for i in {334..346}; do
     echo "" >> residuenumbers.txt
done

## PDB 400 to 429
seq 1842 1924 >> residuenumbers.txt

##Residue 430 is NME cap for linker
## DNA residues
for i in {430..455}; do
     echo "" >> residuenumbers.txt
done
