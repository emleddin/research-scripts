#!/bin/bash
## Begin LSF Directives
#BSUB -P my_proj_id         ## Project ID
#BSUB -W 23:59              ## Wallclock time
#BSUB -nnodes 4             ## number of GPU nodes
#BSUB -alloc_flags gpumps   ## process multiple GPUs
#BSUB -J WT_protein         ## name of job for queue
#BSUB -o WT_prot.out.%J     ## output file name
#BSUB -e WT_prot.err.%J     ## error file name

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

#cd $LSB_SUBCWD

## Folder with files in MemberWork
SCRATCHDIR=$MEMBERWORK/$PROJID/$STRUCT

cd $SCRATCHDIR

## j is last, num is new
j=0
num=1

while [ $num -lt 16 ]; do
## Average 15 ns/day, so submit in chunks of 15
# 16 31 46 61 76 91 106 121 136 151 166 181 196 211 226
# 241 256 271 286 301 316 331 346 361 376

## Replace /path/to with location of executable
## Make sure that the `md.mdin` is your mdin file!
jsrun -n12 -a1 -c4 -g1 -b packed:4 -d packed /path/to/pmemd.cuda.MPI -O -i md.mdin \
 -o ${sys}_md${num}.out \
 -p ${prm} \
 -c ${sys}_md${j}.rst \
 -r ${sys}_md${num}.rst \
 -x ${sys}_md${num}.nc \
 -ref ${sys}_md${j}.rst

num=$[$num+1]
j=$[$j+1]
done

exit 0
