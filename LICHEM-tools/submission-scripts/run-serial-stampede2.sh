#!/bin/bash
#SBATCH -J LICHEM_Opt                    ## Job name
#SBATCH -o lichem.o%j                    ## Name of stdout output file
##SBATCH -e lichem.e%j                   ## Name of stderr error file
#SBATCH -p skx-normal                    ## Queue (partition) name
#SBATCH -N 1                             ## Total # of nodes
#SBATCH -n 1                             ## Total # of mpi tasks (match -N)
#SBATCH -A my_TACC_alloc                 ## Allocation ID
#SBATCH -t 47:59:59                      ## Wallclock time (hh:mm:ss)

## Name of relevant files for LICHEM
## out_tag is used to name the output XYZ and logfile (before .xyz and .log)
## The TINKER key (if present) MUST be tinker.key
## If using Gaussian with QM_basis: GEN keyword, there must be a BASIS file
xyz=xyzfile.xyz
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
export OMP_NUM_THREADS=44

lichem -n 44 -x ${xyz} -c ${connect} -r ${reg} -o ${out_tag}.xyz \
 -l ${out_tag}.log
