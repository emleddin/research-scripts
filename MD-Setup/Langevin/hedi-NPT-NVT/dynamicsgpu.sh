#!/bin/bash
#PBS -q my_gpu_alloc
#PBS -l nodes=n01-02-03
#PBS -j oe
#PBS -r n
#PBS -o err.error
#PBS -N WT_protein_system

export CUDA_VISIBLE_DEVICES=5

cd $PBS_O_WORKDIR

#module load amber/19-cuda_mvapich2
module load amber/19-cuda_serial
#export MV2_ENABLE_AFFINITY=0

e=0
f=1

while [ $f -lt 101 ]; do

#101 for 100 ns NVT-NPT Langevin
#nohup mpirun --bind-to none -np 4 \
#-hostfile $PWD/PBS_NODEFILE 
$AMBERHOME/bin/pmemd.cuda -O -i mdin.4 \
-o WT_protein_system_wat_md$f.out \
-p WT_protein_system_wat.prmtop \
-c WT_protein_system_wat_md$e.rst \
-r WT_protein_system_wat_md$f.rst \
-x WT_protein_system_wat_md$f.mdcrd \
-ref WT_protein_system_wat_md$e.rst

e=$[$e+1]
f=$[$f+1]
done
