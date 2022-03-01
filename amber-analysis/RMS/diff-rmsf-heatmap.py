#!/usr/env/python

## This script reads in MUT and WT RMSF data from cpptraj, averages it, and
##  then creates a defined attribute file for Chimera.

import numpy as np
import pandas as pd

save_name = "MUT-WT-RMSF-map.txt"
data_type = "rmsf"
sys_tag = "MUT-WT-RMSF"

## Read in data
## header = 0 reads header in first row because Python starts at 0
## -------- System 1
d1_a = pd.read_csv('r1/MUT-rmsf-byres.dat', \
 delim_whitespace=True, header=0)
d1_b = pd.read_csv('r2/MUT-rmsf-byres.dat', \
 delim_whitespace=True, header=0)
d1_c = pd.read_csv('r3/MUT-rmsf-byres.dat', \
 delim_whitespace=True, header=0)

## Concatenate together into 1 dataframe
avg_MUT = pd.concat([d1_a, d1_b, d1_c]).groupby(level=0).mean()

## -------- System 2
d2_a = pd.read_csv('r1/WT-rmsf-byres.dat', \
 delim_whitespace=True, header=0)
d2_b = pd.read_csv('r2/WT-rmsf-byres.dat', \
 delim_whitespace=True, header=0)
d2_c = pd.read_csv('r3/WT-rmsf-byres.dat', \
 delim_whitespace=True, header=0)

avg_WT = pd.concat([d2_a, d2_b, d2_c]).groupby(level=0).mean()

## Define the column names to use throughout
## This renames the columns in case something was weird
avg_MUT.columns = ['Frame', 'RMSF']
avg_WT.columns = ['Frame', 'RMSF']

## Create the Difference DataFrame
RMSF_diffs = pd.DataFrame()
#RMSF_diffs = avg_MUT[["Frame"]].copy()
RMSF_diffs["DiffRMSF"] = avg_MUT["RMSF"] - avg_WT["RMSF"]

## Convert to Numpy Array
#RMSF_diffs = RMSF_diffs.to_numpy()

def attribute_to_residues(array, filename, attribute, description,
    matchmode="1-to-1", recipient="residues"):
    """
    Writes an array of values to a chimera-formatted attribute file to allow
    users to color proteins by residue and corresponding value (RMSF, EDA,
    Correlated movement w.r.t specific residue, etc.)
    """
    f = open(filename,"w+")
    f.write("#"+str(description)+"\n")
    f.write("attribute: "+str(attribute)+"\n")
    f.write("match mode: "+matchmode+"\n")
    f.write("recipient: "+recipient+"\n")
    for i in range(len(array)):
        f.write("\t:"+str(i+1)+"\t"+str(array[i])+" \n")
    f.close()
    return

#attribute_to_residues(RMSF_diffs, save_name, data_type, sys_tag)
attribute_to_residues(RMSF_diffs["DiffRMSF"], save_name, data_type, sys_tag)
