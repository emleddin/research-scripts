#!/bin/bash
#PBS -q my_gpu_alloc       ## queue allocation
#PBS -l nodes=n11-12-13    ## use this GPU node
#PBS -j oe                 ## same output and error file
#PBS -r n                  ## Job not re-runnable
#PBS -o err.error          ## name of error file
#PBS -N WT_protein         ## name of job for queue

sys=WT_protein_system_wat
prm=${sys}.prmtop

## Use pairs like 0,1 | 2,3 | 4,5 | 6,7
export CUDA_VISIBLE_DEVICES=0,1

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE  > $PWD/PBS_NODEFILE

module load amber/18-cuda_mvapich2
export MV2_ENABLE_AFFINITY=0

e=0
f=1
while [ $f -lt 501 ]; do

## Make sure that the `md.mdin` is your mdin file!
mpirun -np 2 -hostfile $PWD/PBS_NODEFILE \
$AMBERHOME/bin/pmemd.cuda.MPI -O -i md.mdin \
-o ${sys}_md$f.out \
-p ${prm} \
-c ${sys}_md$e.rst \
-r ${sys}_md$f.rst \
-x ${sys}_md$f.nc \
-ref ${sys}_md$e.rst

e=$[$e+1]
f=$[$f+1]

done
