#!/bin/bash
## Begin LSF Directives
#BSUB -P my_proj_id        ## Project ID
#BSUB -W 23:59             ## Wallclock time
#BSUB -nnodes 1            ## Number of CPU nodes
#BSUB -J WT_protein        ## name of job for queue
#BSUB -o WT_prot.out.%J    ## output file name
#BSUB -e WT_prot.err.%J    ## error file name

## Load modules
module purge
module load gcc/6.4.0 spectrum-mpi cuda \
 cmake readline zlib bzip2 boost python \
 netcdf netcdf-cxx4 netcdf-fortran parallel-netcdf \
 openblas netlib-lapack fftw

## PROJID: Summit allocation name
## sys: the file names used (ex: WT-prot-sy-wat)
## prm: the prmtop name (ex: ${sys}.prmtop == WT-prot-sys-wat.prmtop)
## STRUCT: the folder name for the files in MemberWork

PROJID=my_proj_id
sys=WT_protein_system_wat
prm=${sys}.prmtop
STRUCT='WT_protein'

## Folder with files in MemberWork
SCRATCHDIR=$MEMBERWORK/$PROJID/$STRUCT

cd $SCRATCHDIR

## j is last, num is next
j=0
num=1

while [ $j -lt 16 ]; do
jsrun -n1 -a24 -c24 -g0 $AMBERHOME/bin/pmemd.MPI -O -i mdin.$num \
 -o ${sys}_init${num}.out \
 -p ${prm} \
 -c ${sys}_init${j}.rst \
 -r ${sys}_init${num}.rst \
 -x ${sys}_init${num}.nc \
 -ref ${sys}_init${j}.rst

num=$[$num+1]
j=$[$j+1]
done

exit 0
