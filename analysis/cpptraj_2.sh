#/bin/bash
#PBS -q my_CPU_alloc
#PBS -l nodes=1:ppn=20,mem=20GB
#PBS -j oe
#PBS -r n
#PBS -o err.error
#PBS -N WT_protein_sys_cpp2

cd $PBS_O_WORKDIR

cat $PBS_NODEFILE  > $PWD/PBS_NODEFILE

module load amber/19-mvapich2

mpirun -np 20 -hostfile $PWD/PBS_NODEFILE $AMBERHOME/bin/cpptraj.MPI -p strip.WT_protein_system_wat.prmtop -i cpptraj_analysis.in
