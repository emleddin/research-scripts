#!/bin/bash

#SBATCH --nodes=1          ## request 1 node
#SBATCH -p gpu-shared      ## queue (partition) -- normal, development, etc.
#SBATCH --gres=gpu:2       ## resources you want to use (2 GPUs)
#SBATCH --tasks-per-node=2 ## number GPUs
#SBATCH --export=ALL       ## Keep the current environment stuff
#SBATCH -J WT_protein_sys  ## job name
#SBATCH -A my_comet_alloc  ## PI allocation
#SBATCH -o amber.out       ## output and error file name (%j expands to jobID)
#SBATCH -t 23:59:59        ## run time (hh:mm:ss)

## if you want to submit to gpu rather than gpu-shared
## have number of gpus, tasks per node, and
## OMP_NUM_THREADS equal to 4

## The name of the directory that these files are in
## (used to copy mdinfo to your comet home directory)
prefix=WT_protein_system
sys=WT_protein_system_wat
parm=${sys}.prmtop

## Set up the job environment
module unload intel
module load amber/18
module load cuda

## Set number of threads, should equal number of GPUs
#export OMP_NUM_THREADS=2

## Loop variables to restart calculation
## e=input, f=output
e=0
f=1

## All files should be located in the Lustre filesystem
## So, place them in:
## /oasis/scratch/comet/$USER/temp_project/$prefix

## Copy the necessary files from the submission location
## to the place the job will run
cp $SLURM_SUBMIT_DIR/${parm} /scratch/$USER/$SLURM_JOBID
cp $SLURM_SUBMIT_DIR/*md$e.rst /scratch/$USER/$SLURM_JOBID
cp $SLURM_SUBMIT_DIR/mdin* /scratch/$USER/$SLURM_JOBID

## Access the place to run the job
cd /scratch/$USER/$SLURM_JOBID

## Loop for 100 files (change to files+1)
while [ $f -lt 101 ]; do

## Make sure that the `md.mdin` is your mdin file!
ibrun $AMBERHOME/bin/pmemd.cuda.MPI -O -i md.mdin \
-o ${sys}_md$f.out \
-p ${parm} \
-c ${sys}_md$e.rst \
-r ${sys}_md$f.rst \
-x ${sys}_md$f.nc \
-ref ${sys}_md$e.rst

## Puts time info in home directory
cp mdinfo $HOME/mdinfo.$prefix

## Puts output files into directory accessible outside of job
## Environment--MUST BE IN LOOP
cp /scratch/$USER/$SLURM_JOBID/*md$f.out $SLURM_SUBMIT_DIR
cp -R /scratch/$USER/$SLURM_JOBID/*md$f.rst $SLURM_SUBMIT_DIR
cp -R /scratch/$USER/$SLURM_JOBID/*md$f.nc $SLURM_SUBMIT_DIR

e=$[$e+1]
f=$[$f+1]
done
