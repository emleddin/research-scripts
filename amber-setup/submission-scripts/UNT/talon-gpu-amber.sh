#!/bin/bash

#SBATCH -J WT_protein_system    # name of job in queue
#SBATCH -o Sample_job.o%j       # output file (%j appends job name)
#SBATCH -p public               # partition
#SBATCH --qos general           # quality of service
#SBATCH --ntasks=1              # Number of nodes
#SBATCH --gres=gpu:2            # 2 GPUs
#SBATCH -t 12:00:00             # Wallclock time

## Set up a variable for file naming
sys=WT_protein_system_wat

### Loading modules
module load amber/16-cuda-mpi

e=0
f=1

## Loop for 15 files (change to files+1)
while [ $f -lt 16 ]; do

## Loop for 500 files (change to files+1)
$AMBERHOME/bin/pmemd.cuda -O -i md.mdin \
-o ${sys}_md$f.out \
-p ${sys}.prmtop \
-c ${sys}_md$e.rst \
-r ${sys}_md$f.rst \
-x ${sys}_md$f.nc \
-ref ${sys}_md$e.rst

e=$[$e+1]
f=$[$f+1]
done
