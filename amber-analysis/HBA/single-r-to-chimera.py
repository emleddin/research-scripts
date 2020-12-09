#!/bin/python3
## This takes output from running `rmagic-hbond-avg.r` separately on 2 systems

import sys
import numpy as np
import pandas as pd

command_file = "chimera-commands.dat"

## Read in the WT
df_WT = pd.read_csv('WT_res130_total_hbond_avg.dat', \
 delim_whitespace=True, header=0)
df_MUT = pd.read_csv('MUTA_res130_total_hbond_avg.dat', \
 delim_whitespace=True, header=0)

## Which residue was changed compared to WT per system (SNP location)
## Use the PDB number, not biology
changed_resnum = 130

## The name of the system being compared to WT
system = 'MutA'

## Difference cutoff (ex: 20% different)
cutoff = 20

## Syst 1 Column names:
## tag1a and tag1b (AbsFrac of each system)
## tag1c (AbsDiff of tag1a and tag1b)
tag1a = 'AF_WT'
tag1b = 'AF_MUTA'
tag1c = 'AD_WT_MUTA'

#----------------------------------------------------------------------#
#---------Behind the Curtain: No Need to Modify Past This Line---------#
#----------------------------------------------------------------------#

def clean_dfs(df_WT, df_MUT, tag1a, tag1b, tag1c):
    ## Rename columns -- index is read correctly
    ## Originally ['Acceptor', 'Donor', 'AvgFrac']
    df_WT[tag1a] = df_WT['AvgFrac']
    df_WT = df_WT.drop('AvgFrac', axis=1)
    df_MUT[tag1b] = df_MUT['AvgFrac']
    df_MUT = df_MUT.drop('AvgFrac', axis=1)

    ## Merge the tables together
    tab1ab = df_WT.merge(df_MUT, on=['Acceptor','Donor'])

    ## Explicitly set all NA values as "0"
    tab1ab = tab1ab.fillna(0)

    ## Get the absolute difference between columns
    tab1ab[tag1c] = abs(tab1ab[tag1a] - tab1ab[tag1b])

    ## Change to percentages
    tab1ab[tag1a] = tab1ab[tag1a]*100
    tab1ab[tag1b] = tab1ab[tag1b]*100
    tab1ab[tag1c] = tab1ab[tag1c]*100

    ## Get values > cutoff
    tab1ab_cut = tab1ab[(tab1ab[tag1c] > cutoff)]

    ## Limit to 4 sig figs after decimal
    tab1ab_cleancut = tab1ab_cut.round(decimals=4)
    return tab1ab_cleancut

def chimera_setup(df):
    df.columns = ['Acceptor','Donor','AF_1A','AF_1B','AD_1']
    ## https://stackoverflow.com/questions/50140131/multiple-logical-comparisons-in-pandas-df
    ## Create a clean df to work from with just Acceptors and Donors
    print_df = pd.concat([df.Acceptor,df.Donor], axis=1)
    ## Keep only residue number -- use \' because nucleic acids are mean
    print_df.Acceptor = print_df.Acceptor.str.split('\'').str[0].str\
    .split('_').str[1].str.split('@').str[0]
    print_df.Donor = print_df.Donor.str.split('\'').str[0].str.split('_')\
    .str[1].str.split('@').str[0]
    ## Cycle through for those favored in 1A and those favored in 1B, then
    ## denote that in a new column. If they don't match, label 'NA'
    conditions = [(df.AF_1A > df.AF_1B),
              (df.AF_1A < df.AF_1B)]
    values = ['A', 'B']
    print_df['AF1_favored'] = np.select(conditions, values, 'NA')
    #
    ## Save a dataframe of just the As
    only_1A = print_df[(print_df.AF1_favored == 'A')]
    ## Write a list of donors and a list of acceptors
    A_acc_list = only_1A['Acceptor'].to_list()
    A_don_list = only_1A['Donor'].to_list()
    ## Add these lists into one big list
    only_1A_list = A_acc_list + A_don_list
    ## Remove 'Res' (for Changed_Res) if present
    only_1A_list = [ x for x in only_1A_list if x != 'Res' ]
    ## Turn that list into a set, which removes duplicates
    only_1A_list = set(only_1A_list)
    #
    ## Save a dataframe of just the Bs
    only_1B = print_df[(print_df.AF1_favored == 'B')]
    A_acc_list = only_1B['Acceptor'].to_list()
    A_don_list = only_1B['Donor'].to_list()
    only_1B_list = A_acc_list + A_don_list
    ## Remove 'Res' (for Changed_Res) if present
    only_1B_list = [ x for x in only_1B_list if x != 'Res' ]
    ## Turn that list into a set, which removes duplicates
    only_1B_list = set(only_1B_list)
    #
    ## Get the ones matching only 1A, only 1B, and both using set logic
    only1A = only_1A_list - only_1B_list
    only1B = only_1B_list - only_1A_list
    both1 = only_1A_list.intersection(only_1B_list)
    return only1A, only1B, both1

def write_chimera(command_file, changed_resnum, system, only1A, only1B, both1):
    original_stdout = sys.stdout
    with open(command_file,'w+') as f:
        sys.stdout = f
        f.write("#--------- System: {} ---------#\n".format(system))
        f.write("Favors Sys 1A (WT)\n")
        print('select',*sorted(only1A), sep=" :")
        f.write("\nFavors Sys 1B (SNP)\n")
        print('select',*sorted(only1B), sep=" :")
        f.write("\nIn Both 1A and 1B\n")
        print('select',*sorted(both1), sep=" :")
        f.write("\nSNP Position\n")
        print('select :{}'.format(changed_resnum))
        f.close()
        sys.stdout = original_stdout

## Run the system
tab1ab_cleancut = clean_dfs(df_WT, df_MUT, tag1a, tag1b, tag1c)

only1A, only1B, both1 = chimera_setup(tab1ab_cleancut)

write_chimera(command_file, changed_resnum, system, only1A, only1B, both1)
