#!/bin/bash

##########################################
##        Define your 2 datasets        ##
##      Your h-bond optional cutoff     ##
##         And your outfile name        ##
##########################################

## What two sets are you comparing?
infileA=WT_protein_sys_hbond_avg_WTagain.dat
infileB=MUTA_protein_sys_hbond_avg.dat

## List cutoff as a percentage (20 is 20%)
cutoff=20

## Tag of differences for outfiles
filename=WT-MUTA-hbond-avg
setA=WT_protein_sys
setB=MUTA_protein_sys

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

outfile4AB=hbond-clean-4AB.tmp
outfile5A=hbond-clean-5A.tmp
outfile5B=hbond-clean-5B.tmp
outfile6AB=hbond-clean-6AB.tmp
outfile7AB=hbond-clean-7AB.tmp
outfile8AB=hbond-clean-8AB.tmp
outfile9AB=hbond-clean-9AB.tmp
outfile10AB=$filename-abs-diff.dat

outfile11A=hbond-clean-11A.tmp
outfile11B=hbond-clean-11B.tmp
outfile12A=hbond-clean-12A.tmp
outfile12B=hbond-clean-12B.tmp

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
##           Comparison Time            ##
##      How different are A and B?      ##
##########################################

## Get a list of things with matching Acceptor
## and Donor columns between the files

awk 'FNR==NR{a[$1,$2];next} (($1,$2) in a)' $outfile3A $outfile3B | awk '{printf "%-15s %-15s\n", $1, $2}' > $outfile4AB

## Get rows with the matches from A & B

awk 'FNR==NR{a[$1,$2];next} (($1,$2) in a)' $outfile4AB $outfile3A > $outfile5A
awk 'FNR==NR{a[$1,$2];next} (($1,$2) in a)' $outfile4AB $outfile3B > $outfile5B

## Print a single file for comparison
## 1=Acceptor 2=Donor 3=A_Frac 6=B_Frac

paste $outfile5A $outfile5B > $outfile6AB

## Get the difference and print it

awk 'NR == 1 { print $1 "\t" $2 "\t\t" "Diff_A-B" }; NR>1, NF > 0 { print $1 "\t" $2 "\t" ($6 - $3) }' $outfile6AB > $outfile7AB
(head -n 1 $outfile7AB && tail -n +3 $outfile7AB | cat ) > $filename-diff.dat

## Give it as absolute difference

awk '{ printf "%-15s %-15s %8s\n", $1, $2, ($3 >=0) ? $3 : 0 - $3}' $outfile7AB > $outfile8AB

## Sort the file,
## keeping the header

(head -n 1 $outfile8AB && tail -n +3 $outfile8AB | sort -Rk 3nr ) > $outfile9AB

## Show the difference as a percentage,
## keeping the header

awk 'NR == 1 { print $1 "\t" $2 "\t\t" "AbsDiff%" }; NR>1, NF > 0 { printf "%-15s %-15s %8s\n", $1, $2, $3*100 }' $outfile9AB > $outfile10AB

## Get the cutoff file

awk -v f=$cutoff 'NR == 1; NR > 1 {if ($3> f ) print}' $outfile10AB > $filename-$cutoff-percent-cut.dat

##########################################
##             No Match List            ##
##  No match between A&B?  No problem.  ##
##########################################

## Get full list of no matches
awk 'NR==FNR{a[$1,$2];next} !(($1,$2) in a)' $outfile4AB $outfile3A > $outfile11A
awk 'NR==FNR{a[$1,$2];next} !(($1,$2) in a)' $outfile4AB $outfile3B > $outfile11B


## Print the no matches larger than cutoff
## giving it a header

awk -v f=$cutoff 'NR == 1; NR > 1 {if ( (100*$3) > f ) print}' $outfile11A > $outfile12A
awk -v f=$cutoff 'NR == 1; NR > 1 {if ( (100*$3) > f ) print}' $outfile11B > $outfile12B

## Clean the output and give the files

awk 'NR == 1 { print "#Acceptor" "\t" "Donor" "\t\t" "AbsDiff%" }; NR > 1, NF > 0 { printf "%-15s %-15s %8s\n", $1, $2, $3*100 }' $outfile12A > $setA-$cutoff-percent-cut-nomatch.dat
awk 'NR == 1 { print "#Acceptor" "\t" "Donor" "\t\t" "AbsDiff%" }; NR > 1, NF > 0 { printf "%-15s %-15s %8s\n", $1, $2, $3*100 }' $outfile12B > $setB-$cutoff-percent-cut-nomatch.dat


## Remove the temporary data files
rm hbond-clean*.tmp
