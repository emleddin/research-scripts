#!/bin/bash
#PBS -q my_gpu_alloc      ## Queue allocation
#PBS -l nodes=n01-02-03   ## Specify this GPU node to run on
#PBS -j oe                ## Combine standard output & standard error files
#PBS -r n                 ## Says the job is not rerunnable
#PBS -o err.error         ## Write printed errors to a file titled err.error
#PBS -N WT_prot_sys       ## Name of the job to appear in queue

## Copy final _init*.rst7 as _md0.rst7
sys=WT_protein_system_wat
prm=${sys}.prmtop

## What GPU card to run on
export CUDA_VISIBLE_DEVICES=5

cd $PBS_O_WORKDIR

#module load amber/19-cuda_mvapich2
module load amber/19-cuda_serial
#export MV2_ENABLE_AFFINITY=0

e=0
f=1

#276 for 275 ns NVT-NPT Langevin
while [ $f -lt 276 ]; do

#nohup mpirun --bind-to none -np 4 \
#-hostfile $PWD/PBS_NODEFILE
$AMBERHOME/bin/pmemd.cuda -O \
-i mdin.4 \
-p ${prm} \
-c ${sys}_md$e.rst7 \
-ref ${sys}_md$e.rst7 \
-o ${sys}_md$f.out \
-r ${sys}_md$f.rst7 \
-x ${sys}_md$f.nc


e=$[$e+1]
f=$[$f+1]
done
