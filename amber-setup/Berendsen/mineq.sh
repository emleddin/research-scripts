#!/bin/bash
#PBS -q my_cpu_alloc              ## Use CPUs to run the job
#PBS -l nodes=1:ppn=20,mem=20GB   ## Use 1 node, with 20 processors per node
#PBS -j oe                        ## Combine standard output & std error files
#PBS -r n                         ## Says the jobs is not rerunnable
#PBS -o err.error                 ## Write printed errors to a err.error file
#PBS -N WT_prot_sys               ## Name of the job to appear in queue

## Copy inpcrd to file named like WT_protein_system_wat_init0.rst7
sys=WT_protein_system_wat
prm=${sys}.prmtop

cd $PBS_O_WORKDIR

cat $PBS_NODEFILE  > $PWD/PBS_NODEFILE

module load amber/19-mvapich2

e=0
f=1

while [ $f -lt 11 ]; do

mpirun -np 20 -hostfile $PWD/PBS_NODEFILE $AMBERHOME/bin/pmemd.MPI -O \
-i mdin.$f \
-p ${prm} \
-c ${sys}_init$e.rst7 \
-ref ${sys}_init$e.rst7
-o ${sys}_init$f.out \
-r ${sys}_init$f.rst7 \
-x ${sys}_init$f.nc

e=$[$e+1]
f=$[$f+1]
done
