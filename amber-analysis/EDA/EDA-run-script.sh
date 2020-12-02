#/bin/bash
#PBS -q my_CPU_alloc
#PBS -j oe
#PBS -r n
#PBS -o EDA.error
#PBS -N EDA-run-script

##Load in the Intel compiler
module load intel/17.0

##Access the folder where the files are
cd $PBS_O_WORKDIR

##Compile the EDA program
ifort Residue_E_Decomp_07_15.f90 -o Residue_E_Decomp_07_15.x

##Sleep for 5 seconds, ensuring that the program was compiled
sleep 5

##Run the program; read in the prompt answers [Line 1: Name of input; Line 2: Name of prmtop]
./Residue_E_Decomp_07_15.x < ans.txt

##Acquire the process ID for the program execution
proc_PID=$!

##Wait until the program execution is over before ending the script
wait $proc_PID

echo "All done!"
