#!/bin/bash
# Goal: Add ions to the GROMACS TOP file by editing in-place
## Note: The pre-written blocks came from converting a test
## prmtop/inpcrd with the specific FF that had those ions

top_file=system.top
gro_file=system.gro

#-------- Add Atom Types
## Add after last defined
## AT is last defined, ESC_AT is safe for future sed commands
AT='HW             1   1.008000  0.00000000  A              0              0'
ESC_AT=$(printf '%s\n' "$AT" | sed -e 's/[]\/$*.^[]/\\&/g');

## ION_AT is new AT to add
ION_AT='K+            19  39.100000  0.00000000  A     0.28384033      1.7978874
Cl-           17  35.450000  0.00000000  A     0.48304528     0.05349244'

## ESC_ION_AT is safe for future sed commands
ESC_ION_AT=$(printf '%s\n' "$ION_AT" |
  sed 's/\\/&&/g;s/^[[:blank:]]/\\&/;s/$/\\/')

#------- Define ions as molecule types
## Add in a specific location, then clean that location for sed
add_after="    some line to add after, probably a dihedral"
ESC_add_after=$(printf '%s\n' "$add_after" | sed -e 's/[]\/$*.^[]/\\&/g');

## Add K+ as K
k_block='[ moleculetype ]
; Name            nrexcl
K           3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB      massB
; residue    1 K+ rtp K+ q 1.0
    1         K+      1     K      K       1 1.00000000  39.100000   ; qtot 1.000000'

## Clean K+
ESC_k_block=$(printf '%s\n' "$k_block" |
  sed 's/\\/&&/g;s/^[[:blank:]]/\\&/;s/$/\\/')

## Add Cl-
cl_block='[ moleculetype ]
; Name            nrexcl
CL           3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB      massB
; residue    1 Cl- rtp Cl- q -1.0
    1        Cl-      1    CL     CL       1 -1.00000000  35.450000   ; qtot -1.000000'

## Clean Cl-
ESC_cl_block=$(printf '%s\n' "$cl_block" |
  sed 's/\\/&&/g;s/^[[:blank:]]/\\&/;s/$/\\/')

## Add them in. The %? is crucial for executing the variable
sed -i -e "/$ESC_add_after/a\\
\n\n${ESC_k_block%?}\n\n\n${ESC_cl_block%?}\n" ${top_file}

#---- Additional commands if you want to get funky with WAT

## Deal with GRO
#sed -i -e 's/WAT/SOL/g' ${gro_file}

## Add the WAT and SOL residues (last 2 listings in top)
#NUM=$(tail -n 2 ${top_file} | awk '{ sum += $2 } END {print sum}')
## Remove last 2 lines
#sed -i -e :a -e '$d;N;2,2ba' -e 'P;D' ${top_file}
## Add new last line
#echo "SOL              ${NUM}" >> ${top_file}

## Change WAT molecule to SOL
#sed -i -e 's/WAT/SOL/g' ${top_file}

