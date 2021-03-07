#!/bin/bash

#SBATCH -A my_comet_alloc    ## PI allocation
#SBATCH --nodes=2            ## request 2 nodes
#SBATCH --tasks-per-node=24  ## number CPU
#SBATCH -J WT_protein_system ## job name
#SBATCH -o amber.out         ## output and error file name (%j expands to jobID)
#SBATCH -t 11:59:59          ## run time (hh:mm:ss)
#SBATCH --export=ALL

## mkdir /oasis/scratch/comet/$USER/temp_project/
## NOTE: You should submit this from a directory off of the above path

## The name of the directory that these files are in
## (used to copy mdinfo to your comet home directory)
prefix=WT_protein_system
sys=WT_protein_system_wat
parm=${sys}.prmtop

#Set up the amber environment
module load amber/18

## Copy the necessary files from the submission location
## to the place the job will run
cp $SLURM_SUBMIT_DIR/${parm} /scratch/$USER/$SLURM_JOBID
cp $SLURM_SUBMIT_DIR/*init0.rst /scratch/$USER/$SLURM_JOBID
cp $SLURM_SUBMIT_DIR/mdin* /scratch/$USER/$SLURM_JOBID

## Access the place to run the job
cd /scratch/$USER/$SLURM_JOBID

## Loop variables to restart calculation
## e=input, f=output
e=0
f=1

while [ $f -lt 4 ]; do

ibrun $AMBERHOME/bin/pmemd.MPI -O -i mdin.$f \
-o ${sys}_init$f.out \
-p ${parm} \
-c ${sys}_init$e.rst \
-r ${sys}_init$f.rst \
-x ${sys}_init$f.nc \
-ref ${sys}_init$e.rst

## if calculation will not finish within 48 hours, make sure to
## copy calculation so far to permanent scratch dir INSIDE loop
#cp -R /scratch/$USER/$SLURM_JOBID/* $SLURM_SUBMIT_DIR

## Puts time info in home directory
cp mdinfo $HOME/mdinfo.$prefix

e=$[$e+1]
f=$[$f+1]
done

## these lines copy the files into the submission directory
## after the calculation has finished--make sure to be within
## the wallclock time!
cp -R /scratch/$USER/$SLURM_JOBID/*md$f.out $SLURM_SUBMIT_DIR
cp -R /scratch/$USER/$SLURM_JOBID/*md$f.rst $SLURM_SUBMIT_DIR
cp -R /scratch/$USER/$SLURM_JOBID/*md$f.nc $SLURM_SUBMIT_DIR
