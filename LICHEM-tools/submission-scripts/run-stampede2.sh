#!/bin/bash
#SBATCH -J LICHEM_Job                   ## Job name
#SBATCH -o lichem.o%j                   ## Name of stdout output file
##SBATCH -e lichem.e%j                  ## Name of stderr error file
#SBATCH -p skx-normal                   ## Queue (partition) name
#SBATCH -N 7                            ## Total # of nodes (QSM use beads-2)
#SBATCH -n 7                            ## Total # of mpi tasks (match -N)
#SBATCH -A my_TACC_alloc                ## Allocation ID
#SBATCH -t 47:59:59                     ## Wallclock time (hh:mm:ss)

## Special Instructions for QSM
## QSM requires BurstStruct.xyz which can be made through 2 steps:
##     lichem -path -b 9 -r reactant.xyz -p product.xyz
## (where number of beads can change, here it's 9)
## which makes BeadStartStruct.xyz AND
##     lichem -splitpath -b 9 -p BeadStartStruct.xyz
## Then save a copy from overwrite with:
##     cp BurstStruct.xyz initial-path.xyz
## The reactant XYZ is the XYZ listed with the command
## The number of nodes and MPI tasks specified above should be
## the number of beads -2 (so here, 9-2 = 7)

## Name of relevant files for LICHEM
## out_tag is used to name the output XYZ and logfile (before .xyz and .log)
## The TINKER key (if present) MUST be tinker.key
## If using Gaussian with QM_basis: GEN keyword, there must be a BASIS file
xyz=reactant.xyz
connect=connect.inp
reg=regions.inp
out_tag=my_output_structure

## Load LICHEM module
module load lichem/2020-tinker7.1

## Print helpful info
pwd
date

## -n number of processors
## -x input xyz
## -c connectivity file
## -r regions file
## -o output xyz file
## -l output log file

#--------------------------#
#          SP/DFP          #
#--------------------------#

## Select OpenMP threads
#export OMP_NUM_THREADS=44

#lichem -n 44 -x ${xyz} -c ${connect} -r ${reg} -o ${out_tag}.xyz \
# -l ${out_tag}.log

#--------------------------#
#   QSM or otherwise MPI   #
#--------------------------#

## Select OpenMP threads
export OMP_NUM_THREADS=44

ibrun \
   lichem.MPI -n 44 \
   -x ${xyz} -c ${connect} -r ${reg} \
   -o ${out_tag}.xyz -l ${out_tag}.log
