#!/bin/python3
## 28 June 2021
## This takes `cluster.csv`

import numpy as np
import pandas as pd
from collections import OrderedDict

## Data file
## Header looks like: Node,NodePlot,Window,Cluster,resid
infile = "cluster.csv"

## File with each residue and whether it changes clusters or is unchanged
outfile = "residue-continuity.dat"

## File with cluster information for each window per residue that changes
changed_out = "cluster-changers.dat"

#---------- No need to modify behind the curtain ----------#

## Read in the by-frame clusters
df_c = pd.read_csv(infile, header=0)

## Force column names
df_c.columns = ["Node", "NodePlot", "Window", "Cluster", "ResID"]

## Save the columns you care about
df = df_c[["Window", "Cluster", "ResID"]]

## Create a list of residue names
residues = []
for resname in df.ResID:
    residues.append(resname)

## Keep only unique values in order
res_list = list(dict.fromkeys(residues))

## Create a dictionary of residues
res_dict = OrderedDict()
for residue in res_list:
    res_dict[residue] = df.loc[df['ResID'] == residue ]

## Create an empty dataframe for writing
place_df = pd.DataFrame(columns=['Residue', "Clusters"])

## Test each residue for cluster uniqueness and write the residue to file
##  if multiple clusters
f = open(changed_out, "w")

for residue in res_dict:
    ## Match rows based on the residue name in the dictionary
    test = df.loc[df['ResID'] == residue ]
    ## Check the cluster columns in those rows for uniqueness
    unq = test["Cluster"].nunique()
    ## If more than 1 cluster for the residue
    if unq > 1:
        place_df = place_df.append({'Residue': residue, 'Clusters': 'ChangedClust'}, ignore_index=True)
        f.write(test.to_string())
        f.write("\n\n")
    ## If only 1 cluster for the residue
    else:
        place_df = place_df.append({'Residue': residue, 'Clusters': 'UnchangedClust'}, ignore_index=True)

## Close the opened file
f.close()

## Write the Changed/Unchanged DataFrame to a new CSV
place_df.to_csv(outfile, index=False)
