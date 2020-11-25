#!/bin/bash
## Goal: Add positional restraint sections to correct location
## But via a callable script

sys=topfile-before-extension

m1_last=' the last dihedral for this particular molecule '
ESC_m1_last=$(printf '%s\n' "$m1_last" | sed -e 's/[]\/$*.^[]/\\&/g');

m1_res='; m1 position restraints
#ifdef POSRES_m1
#include "posre_m1.itp"
#endif'

ESC_m1_res=$(printf '%s\n' "$m1_res" |
  sed 's/\\/&&/g;s/^[[:blank:]]/\\&/;s/$/\\/')

## Add m1
sed -i -e "/$ESC_m1_last/a\\
\n\n${ESC_m1_res%?}\n" ${sys}.top

m2_last=' the last dihedral for this particular molecule '
ESC_m2_last=$(printf '%s\n' "$m2_last" | sed -e 's/[]\/$*.^[]/\\&/g');

m2_res='; m2 position restraints
#ifdef POSRES_m2
#include "posre_m2.itp"
#endif'

ESC_m2_res=$(printf '%s\n' "$m2_res" |
  sed 's/\\/&&/g;s/^[[:blank:]]/\\&/;s/$/\\/')

## Add m2
sed -i -e "/$ESC_m2_last/a\\
\n\n${ESC_m2_res%?}\n" ${sys}.top
