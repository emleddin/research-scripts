#!/bin/bash

## Input file titles like: WT_protein_backbone-445-dihedral.dat

## 440,450  has no epsilon or zeta
## 431, 441 has no alpha or beta

## Print the header
printf "%15s %15s %15s %15s %15s %15s %15s %15s\n" "Residue" "Alpha" "Beta" "Gamma" "Delta" "Epsilon" "Zeta" "Chi" > WT_protein_dihedral_backbone.dat

f=431
## Do residue 431 without alpha and beta
GAMMA=$(awk -v f=$f '{ sum += $2 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
DELTA=$(awk -v f=$f '{ sum += $3 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
EPSILON=$(awk -v f=$f '{ sum += $4 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
ZETA=$(awk -v f=$f '{ sum += $5 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
CHI=$(awk -v f=$f '{ sum += $6 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)

printf "%15s %15s %15s %15f %15f %15f %15f %15f\n" "$f" "NA" "NA" "$GAMMA" "$DELTA" "$EPSILON" "$ZETA" "$CHI" >> WT_protein_dihedral_backbone.dat

f=432

## Do residues 432 to 439
while [ $f -lt 440 ]; do

## Take the column average from the cpptraj output
ALPHA=$(awk -v f=$f '{ sum += $2 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
BETA=$(awk -v f=$f '{ sum += $3 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
GAMMA=$(awk -v f=$f '{ sum += $4 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
DELTA=$(awk -v f=$f '{ sum += $5 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
EPSILON=$(awk -v f=$f '{ sum += $6 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
ZETA=$(awk -v f=$f '{ sum += $7 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
CHI=$(awk -v f=$f '{ sum += $8 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)

## Print those averages
printf "%15s %15f %15f %15f %15f %15f %15f %15f\n" "$f" "$ALPHA" "$BETA" "$GAMMA" "$DELTA" "$EPSILON" "$ZETA" "$CHI" >> WT_protein_dihedral_backbone.dat

f=$((f+1))
done

f=440
## Do residue 440 without epsilon and zeta
ALPHA=$(awk -v f=$f '{ sum += $2 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
BETA=$(awk -v f=$f '{ sum += $3 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
GAMMA=$(awk -v f=$f '{ sum += $4 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
DELTA=$(awk -v f=$f '{ sum += $5 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
CHI=$(awk -v f=$f '{ sum += $6 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)

## Print those averages
printf "%15s %15f %15f %15f %15f %15s %15s %15f\n" "$f" "$ALPHA" "$BETA" "$GAMMA" "$DELTA" "NA" "NA" "$CHI" >> WT_protein_dihedral_backbone.dat

f=441
## Do residue 441 without alpha and beta
GAMMA=$(awk -v f=$f '{ sum += $2 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
DELTA=$(awk -v f=$f '{ sum += $3 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
EPSILON=$(awk -v f=$f '{ sum += $4 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
ZETA=$(awk -v f=$f '{ sum += $5 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
CHI=$(awk -v f=$f '{ sum += $6 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)

printf "%15s %15s %15s %15f %15f %15f %15f %15f\n" "$f" "NA" "NA" "$GAMMA" "$DELTA" "$EPSILON" "$ZETA" "$CHI" >> WT_protein_dihedral_backbone.dat

f=442

## Do residues 442 to 449
while [ $f -lt 450 ]; do

## Take the column average from the cpptraj output
ALPHA=$(awk -v f=$f '{ sum += $2 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
BETA=$(awk -v f=$f '{ sum += $3 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
GAMMA=$(awk -v f=$f '{ sum += $4 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
DELTA=$(awk -v f=$f '{ sum += $5 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
EPSILON=$(awk -v f=$f '{ sum += $6 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
ZETA=$(awk -v f=$f '{ sum += $7 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
CHI=$(awk -v f=$f '{ sum += $8 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)

## Print those averages
printf "%15s %15f %15f %15f %15f %15f %15f %15f\n" "$f" "$ALPHA" "$BETA" "$GAMMA" "$DELTA" "$EPSILON" "$ZETA" "$CHI" >> WT_protein_dihedral_backbone.dat

f=$((f+1))
done

f=450
## Do residue 450 without epsilon and zeta
ALPHA=$(awk -v f=$f '{ sum += $2 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
BETA=$(awk -v f=$f '{ sum += $3 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
GAMMA=$(awk -v f=$f '{ sum += $4 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
DELTA=$(awk -v f=$f '{ sum += $5 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)
CHI=$(awk -v f=$f '{ sum += $6 } END {if (NR > 0) print sum / NR }' WT_protein_backbone-${f}-dihedral.dat)

## Print those averages
printf "%15s %15f %15f %15f %15f %15s %15s %15f\n" "$f" "$ALPHA" "$BETA" "$GAMMA" "$DELTA" "NA" "NA" "$CHI" >> WT_protein_dihedral_backbone.dat

## Final average
ALPHA=$(awk '{ sum += $2 } END {if (NR > 0) print sum / NR }' WT_protein_dihedral_backbone.dat)
BETA=$(awk '{ sum += $3 } END {if (NR > 0) print sum / NR }' WT_protein_dihedral_backbone.dat)
GAMMA=$(awk '{ sum += $4 } END {if (NR > 0) print sum / NR }' WT_protein_dihedral_backbone.dat)
DELTA=$(awk '{ sum += $5 } END {if (NR > 0) print sum / NR }' WT_protein_dihedral_backbone.dat)
EPSILON=$(awk '{ sum += $6 } END {if (NR > 0) print sum / NR }' WT_protein_dihedral_backbone.dat)
ZETA=$(awk '{ sum += $7 } END {if (NR > 0) print sum / NR }' WT_protein_dihedral_backbone.dat)
CHI=$(awk '{ sum += $8 } END {if (NR > 0) print sum / NR }' WT_protein_dihedral_backbone.dat)

printf "\n%15s %15f %15f %15f %15f %15s %15s %15f\n" "Avg" "$ALPHA" "$BETA" "$GAMMA" "$DELTA" "$EPSILON" "$ZETA" "$CHI" >> WT_protein_dihedral_backbone.dat
