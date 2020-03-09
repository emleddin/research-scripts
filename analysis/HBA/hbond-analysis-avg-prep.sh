#!/bin/bash

##########################################
##        Define your 2 datasets        ##
##      Your h-bond optional cutoff     ##
##         And your outfile name        ##
##########################################

## What two sets are you comparing?
infileA=WT_protein_system_hbond_avg.dat
infileB=WT_protein_system_hbond_avg.dat
infileC=WT_protein_system_hbond_avg.dat

## Tag of differences for outfiles
filename=WT-protein-123-hbond-avg
setA=WT-1
setB=WT-2
setC=WT-3

##########################################
##         Predefined variables         ##
##########################################

## You can change the file names, but
## it'll be annoying to change the variables

outfile1A=hbond-clean-1A.tmp
outfile2A=hbond-clean-2A.tmp
outfile3A=hbond-clean-3A.tmp

outfile1B=hbond-clean-1B.tmp
outfile2B=hbond-clean-2B.tmp
outfile3B=hbond-clean-3B.tmp

outfile1C=hbond-clean-1C.tmp
outfile2C=hbond-clean-2C.tmp
outfile3C=hbond-clean-3C.tmp

outfile4AB=hbond-clean-4AB.tmp
outfile5ABC=hbond-clean-5ABC.tmp
outfile6A=hbond-clean-6A.tmp
outfile6B=hbond-clean-6B.tmp
outfile6C=hbond-clean-6C.tmp
outfile7ABC=hbond-clean-7ABC.tmp
outfile7ABC2=hbond-clean-7ABC2.tmp
outfile8ABC=$filename.dat


outfile8A=hbond-clean-7A.tmp
outfile8B=hbond-clean-7B.tmp
outfile8C=hbond-clean-7C.tmp
outfile9A=$setA-nomatch.dat
outfile9B=$setB-nomatch.dat
outfile9C=$setC-nomatch.dat

##########################################
##              Fileset A               ##
##  Make some files; do some analysis   ##
##########################################

## Clean the data. Remove lines less than
## 1% and print file with the
## 3 columns you want; keep header
## 1=Acceptor 3=Donor 5=Frac

awk 'NR == 1 {print $1,$3,$5}; NR > 1 { if ($5>0.0099) print $1, $3, $5 }' $infileA > $outfile1A

## Sum duplicate acceptor/donor columns

awk 'NR == 1; {s1[$1,$2] = $1; s2[$1,$2] = $2; s3[$1,$2] += $3} END { for (i in s3) print s1[i], s2[i], s3[i]}' $outfile1A > $outfile2A

## Clean up the output. Make alphabetical order
## by acceptor then by donor and print that
## in clean columns (with left-aligned AAs)

sort $outfile2A | awk '{ printf "%-15s %-15s %8s\n", $1, $2, $3 }' > $outfile3A

##########################################
##              Fileset B               ##
##  Make some files; do some analysis   ##
##########################################

## Clean the data. Remove lines less than
## 1% and print file with the
## 3 columns you want; keep header
## 1=Acceptor 3=Donor 5=Frac

awk 'NR == 1 {print $1,$3,$5}; NR > 1 { if ($5>0.0099) print $1, $3, $5 }' $infileB > $outfile1B

## Sum duplicate acceptor/donor columns

awk 'NR == 1; {s1[$1,$2] = $1; s2[$1,$2] = $2; s3[$1,$2] += $3} END { for (i in s3) print s1[i], s2[i], s3[i]}' $outfile1B > $outfile2B

## Clean up the output. Make alphabetical order
## by acceptor then by donor and print that
## in clean columns (with left-aligned AAs)

sort $outfile2B | awk '{ printf "%-15s %-15s %8s\n", $1, $2, $3 }' > $outfile3B


##########################################
##              Fileset C               ##
##  Make some files; do some analysis   ##
##########################################

## Clean the data. Remove lines less than
## 1% and print file with the
## 3 columns you want; keep header
## 1=Acceptor 3=Donor 5=Frac

awk 'NR == 1 {print $1,$3,$5}; NR > 1 { if ($5>0.0099) print $1, $3, $5 }' $infileC > $outfile1C

## Sum duplicate acceptor/donor columns

awk 'NR == 1; {s1[$1,$2] = $1; s2[$1,$2] = $2; s3[$1,$2] += $3} END { for (i in s3) print s1[i], s2[i], s3[i]}' $outfile1C > $outfile2C

## Clean up the output. Make alphabetical order
## by acceptor then by donor and print that
## in clean columns (with left-aligned AAs)

sort $outfile2C | awk '{ printf "%-15s %-15s %8s\n", $1, $2, $3 }' > $outfile3C

##########################################
##             Average Time             ##
## Let's make averaging A, B, & C easy  ##
##########################################

## Get a list of things with matching Acceptor
## and Donor columns between the files

awk 'FNR==NR{a[$1,$2];next} (($1,$2) in a)' $outfile3A $outfile3B | awk '{printf "%-15s %-15s\n", $1, $2}' > $outfile4AB

awk 'FNR==NR{a[$1,$2];next} (($1,$2) in a)' $outfile4AB $outfile3C | awk '{printf "%-15s %-15s\n", $1, $2}' > $outfile5ABC

## Get rows with the matches from A & B

awk 'FNR==NR{a[$1,$2];next} (($1,$2) in a)' $outfile5ABC $outfile3A > $outfile6A
awk 'FNR==NR{a[$1,$2];next} (($1,$2) in a)' $outfile5ABC $outfile3B > $outfile6B
awk 'FNR==NR{a[$1,$2];next} (($1,$2) in a)' $outfile5ABC $outfile3C > $outfile6C

## Print a single file for comparison
## 1=Acceptor 2=Donor 3=A_Frac 6=B_Frac

paste $outfile6A $outfile6B $outfile6C > $outfile7ABC

## Fix the headers

grep -v '^#Acceptor' $outfile7ABC > $outfile7ABC2
printf "%-15s %-15s %8s %-15s %-15s %8s %-15s %-15s %8s\n" "#Acceptor" "Donor" "Frac" "#Acceptor" "Donor" "Frac" "#Acceptor" "Donor" "Frac" > $outfile8ABC
awk '{printf "%-15s %-15s %8s %-15s %-15s %8s %-15s %-15s %8s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9}' $outfile7ABC2 >> $outfile8ABC

##########################################
##             No Match List            ##
##  No match between AB&C? No problem.  ##
##########################################

## Get full list of no matches
awk 'NR==FNR{a[$1,$2];next} !(($1,$2) in a)' $outfile7ABC $outfile3A > $outfile8A
awk 'NR==FNR{a[$1,$2];next} !(($1,$2) in a)' $outfile7ABC $outfile3B > $outfile8B
awk 'NR==FNR{a[$1,$2];next} !(($1,$2) in a)' $outfile7ABC $outfile3C > $outfile8C


## Print the no matches larger than cutoff
## giving it a header

awk -v f=$cutoff 'NR == 1; NR > 1 {if ( (100*$3) > f ) print}' $outfile8A > $outfile9A
awk -v f=$cutoff 'NR == 1; NR > 1 {if ( (100*$3) > f ) print}' $outfile8B > $outfile9B
awk -v f=$cutoff 'NR == 1; NR > 1 {if ( (100*$3) > f ) print}' $outfile8C > $outfile9C


## Remove the temporary data files
rm hbond-clean*.tmp
