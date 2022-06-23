#!/bin/bash

#------------------#
# Define Variables #
#------------------#

# LICHEM Module
MODULE=lichem/2021.12.03

# Number of atoms in the XYZ
NATOM=12345

# TINKER Parameter File
PRM=amber_tinker_params.prm

# Reactant Directory and optimized XYZ
RXT_DIR=/path/to/reactant/optimization
RXT_XYZ=optimized_reactant.xyz

# Product Directory and optimized XYZ
PROD_DIR=/path/to/product/optimization
PROD_XYZ=optimized_product.xyz

# QSM Directory
QSM_DIR=/path/to/qsm/files

# Number of Beads (including reactant and product!)
BEADS=13

#------------------- No need to modify past the curtain --------------------#

# Get the number of lines in the XYZ file (line w natom, blank, natom lines)
NLINE=$((${NATOM}+2))

#------------#
# Set-up QSM #
#------------#

# Copy the requisite files
cp ${RXT_DIR}/amber_ff14SB_OL15_polk_mg.prm \
   ${RXT_DIR}/BASIS \
   ${RXT_DIR}/connect.inp \
   ${RXT_DIR}/regions.inp \
   ${RXT_DIR}/tinker.key \
   ${RXT_DIR}/vmd-selections.vmd \
   ${QSM_DIR}/.

# Don't copy regions.inp to QSM if it was already modified
#cp ${RXT_DIR}/amber_ff14SB_OL15_polk_mg.prm \
#   ${RXT_DIR}/BASIS \
#   ${RXT_DIR}/connect.inp \
#   ${RXT_DIR}/tinker.key \
#   ${RXT_DIR}/vmd-selections.vmd \
#   ${QSM_DIR}/.

# Copy the final optimization frame as:
# - starting reactant
tail -n ${NLINE} ${RXT_DIR}/${RXT_XYZ} > ${QSM_DIR}/reactant.xyz

# - starting product
tail -n ${NLINE} ${PROD_DIR}/${PROD_XYZ} > ${QSM_DIR}/product.xyz

# Go to QSM directory
cd ${QSM_DIR}

# Load the LICHEM module
printf "\n-----Loading the Module-----\n"
module load ${MODULE}

# Create a path of X beads (X-2 CPUs)
printf "\n-----Creating the Path-----\n"
lichem -path -b ${BEADS} -r reactant.xyz -p product.xyz

# Create BurstStruct.xyz
printf "-----Splitting the Path-----\n"
lichem -splitpath -b ${BEADS} -p BeadStartStruct.xyz

# Save the original path
cp BeadStartStruct.xyz initial-path.xyz

# Reminder
printf " -----Reminder: View BurstStruct.xyz before optimizing!----- \n"
printf " -----Reminder: Change your regions.inp file for QSM!!!----- \n"
printf " -----Reminder: You still need a queue script! :) ----- \n"

