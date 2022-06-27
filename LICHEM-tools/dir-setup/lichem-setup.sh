#!/bin/bash

# Set up a number of frames for QM/MM for 2 systems
# run with
# `bash lichem-setup.sh`
# or
# `bash lichem-setup.sh > setup.log` to save a copy of the output
#
# This script creates a number of subdirectories for `BASE_DIR`,
#  copying and modifying scripts held within another directory (`base-files`).
#
#--------------- Example Directory Tree for BASE_DIR at Start ---------------#
## `$ tree`
# .
# ├── A
# │   ├── base-files
# │   │   ├── met1_protein_system_FF.prm
# │   │   ├── create-reg.py
# │   │   ├── pdbxyz4amber-pmd-params.py
# │   │   ├── run-serial-stampede2.sh
# │   │   ├── tinker.key
# │   │   └── vmd-regions.py
# │   └── write-frames
# │       ├── met1_protein_system_frame_1234.pdb
# │       └── met1_protein_system_frame_9876.pdb
# └── B
# │   ├── base-files
# │   │   ├── met2_protein_system_FF.prm
# │   │   ├── create-reg.py
# │   │   ├── pdbxyz4amber-pmd-params.py
# │   │   ├── run-serial-stampede2.sh
# │   │   ├── tinker.key
# │   │   └── vmd-regions.py
# │   └── write-frames
# │       ├── met2_protein_system_frame_1234.pdb
# │       └── met2_protein_system_frame_9876.pdb
# └── setup-all.sh
#
#--------------- Python Script Mods Allowing Arguments ---------------#
# The Python scripts need to accept system arguments!
# Make sure the following are in the scripts:
#
## In `create-reg.py`:
# import sys
#
# # The script itself is sys.argv[0]
# orig_pdb = sys.argv[1]
# tink_xyz = sys.argv[2]
#
## In `pdbxyz4amber-pmd-params.py`
# import sys
#
# # The script itself is sys.argv[0]
# infile = sys.argv[1]
# outfile = sys.argv[2]

#--------------- Set-Up Variables ---------------#
# Base directory with A_DIR and B_DIR as subdirectories
BASE_DIR=/path/to/qmmm

# Path to each metal's files
A_DIR=${BASE_DIR}/MET1
B_DIR=${BASE_DIR}/MET2

# Tinker parameter files for each system
A_PARAM=met1_protein_system_FF.prm
B_PARAM=met2_protein_system_FF.prm

# Subdirectory names within each metal (Ex: ${A_DIR}/${BF_DIR_NAME})
## Where all the scripts and common files are located (e.g., `tinker.key`)
BF_DIR_NAME=base-files
## Where all the PDBs are located
PDB_DIR_NAME=write-frames

# Name of the first folder within the FRAME folder (for converting PDB!)
ONE_DIR_NAME=1-xyz-conversion
# The PDBXYZ script (with sys.argv!)
PDBXYZ=pdbxyz4amber-pmd-params.py

# Name of the second folder within the FRAME folder (for creating regions file!)
TWO_DIR_NAME=2-lich-conversion
# The create-regions script tailored to the specific system (with sys.argv!)
CREATE_REG=create-reg.py
# The vmd-regions script
VMD_REG=vmd-regions.py

# Name of the third folder within the FRAME folder (for single point test!)
THREE_DIR_NAME=3-SP

# List of frames for system A - in parentheses!
list_A_frames=(1234 9876)

# List of frames for system B - in parentheses!
list_B_frames=(1234 9876)

#---- File names!

# This script assumes that the start PDBs for system A are ${A_TAG}_${FRAME}.pdb
#  and start PDBs for system B are ${B_TAG}_${FRAME}.pdb
#  Ex: met1_protein_system_frame_1234.pdb
A_TAG=met1_protein_system_frame
B_TAG=met2_protein_system_frame

# Tinker key file
T_KEY=tinker.key

# HPC script to use as a base
HPC_SCRIPT=run-serial-stampede2.sh

# Short (4-characters or less) name for each system for the HPC script job
A_NAME=MET1
B_NAME=MET2

#----------------------------------------------#
#--------------- Run the Script ---------------#
#----------------------------------------------#

# Load the LICHEM module
module load lichem/2021.12.03

#--- Typical users will not need to modify past this point!
# Atypical cases:
# - Not using SLURM (change the `sed` search)
# - The default output names were changed in `create-reg.py`
# - Wanting to clear a previous attempt

################
# Clean Priors #
################
# rm -rf ${A_DIR}/f_* || true
# rm -rf ${B_DIR}/f_* || true

######
# A #
######

# Define outside loop because same for entire loop
BF_DIR=${A_DIR}/${BF_DIR_NAME}
PDB_DIR=${A_DIR}/${PDB_DIR_NAME}

for FRAME in "${list_A_frames[@]}"
do
  FRAME_DIR=${A_DIR}/f_${FRAME}
  ONE_DIR=${FRAME_DIR}/${ONE_DIR_NAME}
  TWO_DIR=${FRAME_DIR}/${TWO_DIR_NAME}
  THREE_DIR=${FRAME_DIR}/${THREE_DIR_NAME}

  # Where we are in the loop
  echo $FRAME_DIR

  # Create subdirectories
  mkdir -p ${ONE_DIR}
  mkdir -p ${TWO_DIR}
  mkdir -p ${THREE_DIR}

  cd ${A_DIR}

  # Copy relevant files to first
  cp ${PDB_DIR}/${A_TAG}_${FRAME}.pdb ${ONE_DIR}/
  cp ${BF_DIR}/${A_PARAM} ${ONE_DIR}/.
  cp ${BF_DIR}/${PDBXYZ} ${ONE_DIR}/.

  # Enter subdirectory and run script
  cd ${ONE_DIR}
  # python3 script infile outfile
  python3 ${PDBXYZ} ${A_TAG}_${FRAME}.pdb ${A_TAG}_${FRAME}.xyz

  # Status report
  echo "Completed PDBXYZ for ${FRAME}"

  # Backup to frame
  cd ${FRAME_DIR}

  # Copy relevant files to second
  cp ${BF_DIR}/${CREATE_REG} ${TWO_DIR}/.
  cp ${BF_DIR}/${VMD_REG} ${TWO_DIR}/.
  cp ${BF_DIR}/${A_PARAM} ${TWO_DIR}/.
  cp ${BF_DIR}/${T_KEY} ${TWO_DIR}/.
  cp ${ONE_DIR}/${A_TAG}_${FRAME}.pdb ${TWO_DIR}/.
  cp ${ONE_DIR}/${A_TAG}_${FRAME}.xyz ${TWO_DIR}/.

  cd ${TWO_DIR}
  python3 ${CREATE_REG} ${A_TAG}_${FRAME}.pdb ${A_TAG}_${FRAME}.xyz

  # Convert to LICHEM format
  lichem -convert -t ${A_TAG}_${FRAME}.xyz -k ${T_KEY}
  # Overwrite blank regions file
  cp regions.inp_backup regions.inp

  # Create a VMD-selections file
  python3 ${VMD_REG}

  # Status report
  echo "Completed LICHEM conversion for ${FRAME}"

  # Backup to frame
  cd ${FRAME_DIR}

  # Copy relevant files to third
  cp ${BF_DIR}/${A_PARAM} ${THREE_DIR}/.
  cp ${BF_DIR}/${HPC_SCRIPT} ${THREE_DIR}/.
  cp ${BF_DIR}/${T_KEY} ${THREE_DIR}/.
  cp ${TWO_DIR}/BASIS ${THREE_DIR}/.
  cp ${TWO_DIR}/connect.inp ${THREE_DIR}/.
  cp ${TWO_DIR}/regions.inp ${THREE_DIR}/.
  cp ${TWO_DIR}/xyzfile.xyz ${THREE_DIR}/.
  cp ${TWO_DIR}/vmd-selections.vmd ${THREE_DIR}/.

  cd ${THREE_DIR}

  # Fix the queue name and  out_tag in script
  sed -i -e "s/#SBATCH -J LICHEM_Opt/#SBATCH -J ${A_NAME%?}_${FRAME%?}_SP/g" \
     ${HPC_SCRIPT}
  sed -i -e "s/out_tag=my_output_structure/out_tag=${A_TAG%?}_${FRAME%?}_out/g" \
     ${HPC_SCRIPT}

  echo "Completed SP set-up for ${FRAME}"

done

######
# B #
######

# Define outside loop because same for entire loop
BF_DIR=${B_DIR}/${BF_DIR_NAME}
PDB_DIR=${B_DIR}/${PDB_DIR_NAME}

for FRAME in "${list_B_frames[@]}"
do
  FRAME_DIR=${B_DIR}/f_${FRAME}
  ONE_DIR=${FRAME_DIR}/${ONE_DIR_NAME}
  TWO_DIR=${FRAME_DIR}/${TWO_DIR_NAME}
  THREE_DIR=${FRAME_DIR}/${THREE_DIR_NAME}

  # Where we are in the loop
  echo $FRAME_DIR

  # Create subdirectories
  mkdir -p ${ONE_DIR}
  mkdir -p ${TWO_DIR}
  mkdir -p ${THREE_DIR}

  cd ${B_DIR}

  # Copy relevant files to first
  cp ${PDB_DIR}/${B_TAG}_${FRAME}.pdb ${ONE_DIR}/.
  cp ${BF_DIR}/${B_PARAM} ${ONE_DIR}/.
  cp ${BF_DIR}/${PDBXYZ} ${ONE_DIR}/.

  # Enter subdirectory and run script
  cd ${ONE_DIR}
  # python3 script infile outfile
  python3 ${PDBXYZ} ${B_TAG}_${FRAME}.pdb ${B_TAG}_${FRAME}.xyz

  # Status report
  echo "Completed PDBXYZ for ${FRAME}"

  # Backup to frame
  cd ${FRAME_DIR}

  # Copy relevant files to second
  cp ${BF_DIR}/${CREATE_REG} ${TWO_DIR}/.
  cp ${BF_DIR}/${VMD_REG} ${TWO_DIR}/.
  cp ${BF_DIR}/${B_PARAM} ${TWO_DIR}/.
  cp ${BF_DIR}/${T_KEY} ${TWO_DIR}/.
  cp ${ONE_DIR}/${B_TAG}_${FRAME}.pdb ${TWO_DIR}/.
  cp ${ONE_DIR}/${B_TAG}_${FRAME}.xyz ${TWO_DIR}/.

  cd ${TWO_DIR}
  python3 ${CREATE_REG} ${B_TAG}_${FRAME}.pdb ${B_TAG}_${FRAME}.xyz

  # Convert to LICHEM format
  lichem -convert -t ${B_TAG}_${FRAME}.xyz -k ${T_KEY}
  # Overwrite blank regions file
  cp regions.inp_backup regions.inp

  # Create a VMD-selections file
  python3 ${VMD_REG}

  # Status report
  echo "Completed LICHEM conversion for ${FRAME}"

  # Backup to frame
  cd ${FRAME_DIR}

  # Copy relevant files to third
  cp ${BF_DIR}/${B_PARAM} ${THREE_DIR}/.
  cp ${BF_DIR}/${HPC_SCRIPT} ${THREE_DIR}/.
  cp ${BF_DIR}/${T_KEY} ${THREE_DIR}/.
  cp ${TWO_DIR}/BASIS ${THREE_DIR}/.
  cp ${TWO_DIR}/connect.inp ${THREE_DIR}/.
  cp ${TWO_DIR}/regions.inp ${THREE_DIR}/.
  cp ${TWO_DIR}/xyzfile.xyz ${THREE_DIR}/.
  cp ${TWO_DIR}/vmd-selections.vmd ${THREE_DIR}/.

  cd ${THREE_DIR}

  # Fix the queue name and  out_tag in script
  sed -i -e "s/#SBATCH -J LICHEM_Opt/#SBATCH -J ${B_NAME%?}_${FRAME%?}_SP/g" \
     ${HPC_SCRIPT}
  sed -i -e "s/out_tag=my_output_structure/out_tag=${B_TAG%?}_${FRAME%?}_out/g" \
     ${HPC_SCRIPT}

  echo "Completed SP set-up for ${FRAME}"

done

# Function to print a statement to stderr
echoerr() { echo "$@" 1>&2; }

echoerr I did not check the Tinker XYZs with analyze!!!
echoerr If you run into issues running a single point, start there!
