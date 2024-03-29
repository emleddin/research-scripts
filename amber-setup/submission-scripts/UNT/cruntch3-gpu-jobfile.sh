#!/bin/bash
#PBS -q my_gpu_alloc       ## queue allocation
#PBS -l nodes=n11-12-13    ## use this GPU node
#PBS -j oe                 ## same output and error file
#PBS -r n                  ## Job not re-runnable
#PBS -o err.error          ## name of error file
#PBS -N WT_protein         ## name of job for queue

## Copy inpcrd to file named like WT_protein_system_wat_md0.rst
sys=WT_protein_system_wat
prm=${sys}.prmtop

## The specific GPU you're running on
export CUDA_VISIBLE_DEVICES=0

cd $PBS_O_WORKDIR

## You need a cuda_serial version for GPU runs
module load amber/19-cuda_serial

e=0
f=1

## Loop for 500 files (change to files+1)
while [ $f -lt 501 ]; do

## Make sure that the `md.mdin` is your mdin file!
$AMBERHOME/bin/pmemd.cuda -O -i md.mdin \
-o ${sys}_prod$f.out \
-p ${prm} \
-c ${sys}_prod$e.rst \
-r ${sys}_prod$f.rst \
-x ${sys}_prod$f.nc \
-ref ${sys}_prod$e.rst

## Check that rst was made, if not break loop
if [ -f "${sys}_prod${f}.rst" ]; then
    :
else
    break
fi

e=$[$e+1]
f=$[$f+1]
done
