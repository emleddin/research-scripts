#!/bin/bash
#PBS -q my_cpu_alloc              ## queue allocation
#PBS -l nodes=1:ppn=20,mem=20GB   ## 20 processors, 80 GB memory
#PBS -j oe                        ## same output and error file
#PBS -r n                         ## Job not re-runnable
#PBS -o err.error                 ## name of error file
#PBS -N WT_protein                ## name of job for queue

sys=WT_protein_system_wat
prm=${sys}.prmtop

cd $PBS_O_WORKDIR

cat $PBS_NODEFILE  > $PWD/PBS_NODEFILE

#module load amber/16.7-mvapich2
#module load amber/18-mvapich2
module load amber/19-mvapich2

e=0
f=1
while [ $f -lt 4 ]; do

mpirun -np 20 -hostfile $PWD/PBS_NODEFILE $AMBERHOME/bin/pmemd.MPI \
-O -i mdin.$f \
-o ${sys}_init$f.out \
-p ${prm} \
-c ${sys}_init$e.rst \
-r ${sys}_init$f.rst \
-x ${sys}_init$f.nc \
-ref ${sys}_init$e.rst

e=$[$e+1]
f=$[$f+1]
done
