#!/bin/bash

#------------------#
# Define Variables #
#------------------#

# A is previous QSM opt, B is new opt
A_DIR=/path/to/first/opt/qsm/files
A_XYZ=qsm_optimized_path.xyz

B_DIR=/path/to/second/opt/qsm/files

# TINKER Parameter File
PRM=amber_tinker_params.prm

# Number of atoms in the reactant.xyz
NATOM=12345

# Number of Beads
BEADS=13

#------------------- No need to modify past the curtain --------------------#

# Get the number of lines in the QSM XYZ file
#  (line w natom*beads, blank, natom*beads lines)
NLINE=$((${NATOM}*${BEADS}+2))

#------------#
# Set-up QSM #
#------------#

# Copy the requisite files
cp ${A_DIR}/${PRM} \
   ${A_DIR}/BASIS \
   ${A_DIR}/connect.inp \
   ${A_DIR}/regions.inp \
   ${A_DIR}/tinker.key \
   ${A_DIR}/vmd-selections.vmd \
   ${B_DIR}/.

# Don't re-copy regions.inp if it was already modified
#cp ${A_DIR}/${PRM} \
#   ${A_DIR}/BASIS \
#   ${A_DIR}/connect.inp \
#   ${A_DIR}/tinker.key \
#   ${A_DIR}/vmd-selections.vmd \
#   ${B_DIR}/.

# Copy the final optimization path as initial BeadStartStruct.xyz
tail -n ${NLINE} ${A_DIR}/${A_XYZ} > ${B_DIR}/BeadStartStruct.xyz

# Copy the final split path as initial split path
cp ${A_DIR}/BurstStruct_1.xyz ${B_DIR}/BurstStruct.xyz

# Copy reactant because of frozen ends
cp ${A_DIR}/reactant.xyz ${B_DIR}/reactant.xyz

# Reminders
printf " -----Reminder: Update the regions.inp file! \n"
printf " -----          Remove the force constant and change tols! \n"
printf " -----Reminder: You still need a queue script! :) \n"

