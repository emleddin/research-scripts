#!/bin/bash
#SBATCH --partition=gpu-shared    ## queue (partition) -- normal, development, etc.
#SBATCH --nodes=1                 ## request 1 node (don't change this)
#SBATCH --gpus=2                  ## resources you want to use (1-4 GPUs)
#SBATCH --ntasks-per-node=2
# --cpus-per-task=1         ## How many CPUs per GPU
#SBATCH --mem=120G                ## How much memory (up to ~90G/GPU)
#SBATCH --export=none             ## Don't keep current environment stuff
#SBATCH -A my_expanse_alloc       ## PI allocaton
#SBATCH -t 47:59:59               ## run time (hh:mm:ss)
#SBATCH -J WT_protein_sys         ## job name
#SBATCH -o amber.out              ## output & error file name (%j expands to jobID)

## All files should be located in the Lustre filesystem
## So, place them in:
## /expanse/lustre/scratch/$USER/temp_project/$prefix

## The name of the directory that these files are in
## (used to copy mdinfo to your expanse home directory)
prefix=WT_protein_system
sys=WT_protein_system_wat
parm=${sys}.prmtop

## Loop variables to restart calculation
## e=input, f=output
e=0
f=1

## Set up the job environment
module purge
module load slurm
module load gpu
module load openmpi
module load amber

## Copy the necessary files from the submission location
## to the place the job will run
cp $SLURM_SUBMIT_DIR/${parm} /scratch/$USER/job_$SLURM_JOBID
cp $SLURM_SUBMIT_DIR/*md$e.rst /scratch/$USER/job_$SLURM_JOBID
cp $SLURM_SUBMIT_DIR/md.mdin /scratch/$USER/job_$SLURM_JOBID

## Access the place to run the job
cd /scratch/$USER/job_$SLURM_JOBID

while [ $f -lt 276 ]; do

## Make sure that the `md.mdin` is your mdin file!
pmemd.cuda -O -i md.mdin \
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
cp /scratch/$USER/job_$SLURM_JOBID/*md$f.out $SLURM_SUBMIT_DIR
cp -R /scratch/$USER/job_$SLURM_JOBID/*md$f.rst $SLURM_SUBMIT_DIR
cp -R /scratch/$USER/job_$SLURM_JOBID/*md$f.nc $SLURM_SUBMIT_DIR

e=$[$e+1]
f=$[$f+1]
done
