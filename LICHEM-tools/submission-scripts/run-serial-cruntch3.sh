#!/bin/bash

#PBS -N LICHEM_OPT
#PBS -j oe
#PBS -o qmmm.err
#PBS -q my_cpu_alloc
#PBS -l nodes=1:ppn=20,mem=80GB
#PBS -r n

## Name of relevant files for LICHEM
## The TINKER key (if present) MUST be tinker.key
## If using Gaussian with QM_basis: GEN keyword, there must be a BASIS file
xyz=xyzfile.xyz
connect=connect.inp
reg=regions.inp

## out is used to name the output XYZ and logfile (before .xyz and .log)
out=my_output_structure

cd $PBS_O_WORKDIR

## Load LICHEM
module load lichem/2020.12.10

## Set scratch
export GAUSS_SCRDIR="/tmp"

## -n number of processors
## -x input xyz
## -c connectivity file
## -r regions file
## -o output xyz file
## -l output log file
lichem -n 20 \
 -x ${xyz} \
 -c ${connect} \
 -r ${reg} \
 -o ${out}.xyz \
 -l ${out}.log
