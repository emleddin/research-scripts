#!/bin/bash
#PBS -N LICHEM_QSM
#PBS -j oe
#PBS -o qmmm.err
#PBS -q my_mpi_alloc
#PBS -l nodes=5:ppn=20,mem=100GB
#PBS -r n

#------------ Preparation for QSM -------------#
## The nodes requested above should ideally be "beads - 2"
## So for 7 beads, use 5 nodes.
## Try to use resources effectively, so if you want 11 beads,
## You could still go for 5 nodes. R1: 5 cores, R2: 4 cores

## Create BeadStartStruct.xyz with number of beads (-b)
## Can add possible TS state (-t)
#lichem -path -b 7 -r reactant.xyz  -p product.xyz
#lichem -path -b 7 -r reactant.xyz -t ts.xyz -p product.xyz

## Create BurstStruct.xyz
#lichem -splitpath -b 7 -p BeadStartStruct.xyz

## Copy BeadStartStruct.xyz so you remember what it is
#cp BeadStartStruct.xyz initial-path.xyz
#----------------------------------------------#

## Name of relevant files for LICHEM
## The TINKER key (if present) MUST be tinker.key
## If using Gaussian with QM_basis: GEN keyword, there must be a BASIS file
## It will look for BeadStartStruct.xyz and BurstStruct.xyz
react_xyz=reactant.xyz
connect=connect.inp
reg=regions.inp

## Output XYZ/Log
out=my_output_structure

## Access job directory
cd $PBS_O_WORKDIR

## Load LICHEM
module load lichem/2020.12.10

## Set scratch
export GAUSS_SCRDIR="/tmp"

# cat $PBS_NODEFILE > all_hosts

#     replace YYY below,
#     with the given ppn value above
#     i.e. if
#          PBS -l nodes=5:ppn=20,mem=MMM
#          then
#          /bin/sed -n '1~20p' all_hosts > host_list
/bin/sed -e 's/$/-ib/' ${PBS_NODEFILE} > all_hosts
/bin/sed -n '1~20p' all_hosts > host_list

## -n number of processors
## -x input xyz
## -c connectivity file
## -r regions file
## -o output xyz file
## -l output log file
mpirun -n 5 --hostfile host_list \
lichem.MPI -n 20 \
 -x ${react_xyz} \
 -c ${connect} \
 -r ${reg} \
 -o ${out}.xyz \
 -l ${out}.log


