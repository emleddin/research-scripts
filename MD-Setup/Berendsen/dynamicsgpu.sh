#!/bin/bash
#PBS -q my_gpu_alloc      ## Queue allocation
#PBS -l nodes=n01-02-03   ## Specify this GPU node to run on
#PBS -j oe                ## Combine standard output & standard error files
#PBS -r n                 ## Says the job is not rerunnable
#PBS -o err.error         ## Write printed errors to a file titled err.error
#PBS -N WT_prot_sys       ## Name of the job to appear in queue

export CUDA_VISIBLE_DEVICES=4

cd $PBS_O_WORKDIR

#module load amber/19-cuda_mvapich2
module load amber/19-cuda_serial
#export MV2_ENABLE_AFFINITY=0

e=0
f=1

while [ $f -lt 201 ]; do

#nohup mpirun --bind-to none -np 4 \
#-hostfile $PWD/PBS_NODEFILE 
$AMBERHOME/bin/pmemd.cuda -O -i mdin.11 \
-o WT_protein_system_wat_md$f.out \
-p WT_protein_system_wat.prmtop \
-c WT_protein_system_wat_md$e.rst \
-r WT_protein_system_wat_md$f.rst \
-x WT_protein_system_wat_md$f.mdcrd \
-ref WT_protein_system_wat_md$e.rst

e=$[$e+1]
f=$[$f+1]
done
