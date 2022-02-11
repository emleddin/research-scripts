#!/usr/env/python3
## Provide robust information for a QM PDB (created using mda-qm-part1.py)

import pandas as pd

## Read in the space-delimited BASIS_verification file
df = pd.read_csv('BASIS_verification.txt', sep=r'\s{1,}', engine='python', \
 header=0)

## Remove the PB rows
df = df[df.Type != 'PB']

## Reindex, then force index to start at '1' instead of '0'
df = df.reset_index()
df.index += 1

## Add a column with the new QM PDB indexing
df['QM_PDB'] = df.index

## Write a new file with the QM_PDB_ID indexing
with open("PDB_verification.txt", "w+") as bv_out:
    bv_out.write("QM_PDB_ID AtomName ResName ResNum Regions_ID TINKER_ID BASIS_ID Type\n")
    for r in df.itertuples(index=True, name='Pandas'):
        bv_out.write("{:<9} {:<8} {:<7} {:<6} {:<10} {:<9} {:<8} {:<5}\n".format(\
        r.QM_PDB, r.AtomName, r.ResName, r.ResNum, r.Regions_ID, r.TINKER_ID,\
         r.BASIS_ID, r.Type))
    bv_out.close()
