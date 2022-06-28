#!/bin/bash

#------------------#
# Define Variables #
#------------------#

# A is previous DFP opt, B is cont opt
A_DIR=/path/to/initial/dfp/opt/files
A_XYZ=optimized_path.xyz
A_QUEUE=run-serial-stampede2.sh

B_DIR=/path/to/continued/dfp/opt/files

# TINKER Parameter File
PRM=amber_tinker_params.prm

# Number of atoms in the reactant.xyz
NATOM=12345

#------------------- No need to modify past the curtain --------------------#

# Get the number of lines in the XYZ file
#  (line w natom, blank, natom lines)
NLINE=$((${NATOM}+2))

#------------#
# Set-up DFP #
#------------#

# Copy the requisite files
cp ${A_DIR}/${PRM} \
   ${A_DIR}/BASIS \
   ${A_DIR}/connect.inp \
   ${A_DIR}/regions.inp \
   ${A_DIR}/tinker.key \
   ${A_DIR}/vmd-selections.vmd \
   ${A_DIR}/${A_QUEUE} \
   ${B_DIR}/.

# Don't re-copy regions.inp if it was already modified
#cp ${A_DIR}/${PRM} \
#   ${A_DIR}/BASIS \
#   ${A_DIR}/connect.inp \
#   ${A_DIR}/tinker.key \
#   ${A_DIR}/vmd-selections.vmd \
#   ${A_DIR}/${A_QUEUE} \
#   ${B_DIR}/.

# Copy the final optimization step for start
tail -n ${NLINE} ${A_DIR}/${A_XYZ} > ${B_DIR}/xyzfile.xyz

# Reminders
printf " -----Reminder: Modify your queue script! :) \n"

