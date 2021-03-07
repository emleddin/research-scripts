#!/bin/bash

#SBATCH -J Tutorial1            ## name of job in queue
#SBATCH -o Sample_job.o%j       ## output file (%j appends job name)
#SBATCH -p public               ## partition
#SBATCH --qos general           ## quality of service
#SBATCH -N 1                    ## Number of nodes
#SBATCH -n 16                   ## Number of processors
#SBATCH -t 12:00:00             ## Wallclock time

sys=polyAT_vac

## Loading modules
module load amber/16-gen

## Run the job -- ${sys}_init_min.in is the mdin file
mpirun -np 16 sander.MPI \
 -O -i ${sys}_init_min.in \
 -o ${sys}_init_min.out \
 -c ${sys}.inpcrd \
 -p ${sys}.prmtop \
 -r ${sys}_init_min.rst
