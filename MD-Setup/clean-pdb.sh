#!/bin/bash

## Define your files
thing=RCSB.pdb
thing_clean=RCSB_clean.pdb

## Clean the RCSB PDB (Remove crystal junk)
grep -e '^ATOM\|^HETATM\|^TER\|^END' $thing > $thing_clean

## The remainder will only act on lines containing the specified info

## Remove any B chain lines
## OSX needs to be sed -i "" BLAH to work
sed -i '/BALA/d' $thing_clean
sed -i '/BARG/d' $thing_clean
sed -i '/BASN/d' $thing_clean
sed -i '/BASP/d' $thing_clean
sed -i '/BCYS/d' $thing_clean
sed -i '/BGLU/d' $thing_clean
sed -i '/BGLN/d' $thing_clean
sed -i '/BGLY/d' $thing_clean
sed -i '/BHIS/d' $thing_clean
sed -i '/BILE/d' $thing_clean
sed -i '/BLEU/d' $thing_clean
sed -i '/BLYS/d' $thing_clean
sed -i '/BMET/d' $thing_clean
sed -i '/BPHE/d' $thing_clean
sed -i '/BPRO/d' $thing_clean
sed -i '/BSER/d' $thing_clean
sed -i '/BTHR/d' $thing_clean
sed -i '/BTRP/d' $thing_clean
sed -i '/BTYR/d' $thing_clean
sed -i '/BVAL/d' $thing_clean

## Rename A chain as only chain
sed -i 's/AALA/ ALA/g' $thing_clean
sed -i 's/AARG/ ARG/g' $thing_clean
sed -i 's/AASN/ ASN/g' $thing_clean
sed -i 's/AASP/ ASP/g' $thing_clean
sed -i 's/ACYS/ CYS/g' $thing_clean
sed -i 's/AGLU/ GLU/g' $thing_clean
sed -i 's/AGLN/ GLN/g' $thing_clean
sed -i 's/AGLY/ GLY/g' $thing_clean
sed -i 's/AHIS/ HIS/g' $thing_clean
sed -i 's/AILE/ ILE/g' $thing_clean
sed -i 's/ALEU/ LEU/g' $thing_clean
sed -i 's/ALYS/ LYS/g' $thing_clean
sed -i 's/AMET/ MET/g' $thing_clean
sed -i 's/APHE/ PHE/g' $thing_clean
sed -i 's/APRO/ PRO/g' $thing_clean
sed -i 's/ASER/ SER/g' $thing_clean
sed -i 's/ATHR/ THR/g' $thing_clean
sed -i 's/ATRP/ TRP/g' $thing_clean
sed -i 's/ATYR/ TYR/g' $thing_clean
sed -i 's/AVAL/ VAL/g' $thing_clean

## Rename HOH as WAT
sed -i 's/HOH/WAT/g' $thing_clean

## Optional: Delete inhibitor/artifact lines
## Check the "Small Molecules" section of RCSB
## And compare against paper/common inhibitors for enzyme class
#sed -i '/ACT /d' $thing_clean
#sed -i '/OGA /d' $thing_clean
#sed -i '/SO4 /d' $thing_clean
#sed -i '/GOL /d' $thing_clean
#sed -i '/FSU /d' $thing_clean
#sed -i '/EDO /d' $thing_clean
