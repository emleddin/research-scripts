#/bin/bash
#PBS -q my_cpu_alloc
#PBS -l nodes=1:ppn=20,mem=20GB
#PBS -j oe
#PBS -r n
#PBS -o err.error
#PBS -N WT_protein_system

cd $PBS_O_WORKDIR

cat $PBS_NODEFILE  > $PWD/PBS_NODEFILE

module load amber/19-mvapich2

e=0
f=1
while [ $f -lt 4 ]; do

mpirun -np 20 -hostfile $PWD/PBS_NODEFILE $AMBERHOME/bin/pmemd.MPI -O -i mdin.$f \
-o WT_protein_system_wat_init$f.out \
-p WT_protein_system_wat.prmtop \
-c WT_protein_system_wat_init$e.rst \
-r WT_protein_system_wat_init$f.rst \
-x WT_protein_system_wat_init$f.mdcrd \
-ref WT_protein_system_wat_init$e.rst

e=$[$e+1]
f=$[$f+1]
done

