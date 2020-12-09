#!/bin/python3
## Use with a table of 5 systems from system-hbond-table.r

import sys
import numpy as np
import pandas as pd
from tabulate import tabulate

outfile = "filtered-table.dat"
command_file = "chimera-commands.dat"

df = pd.read_csv('attempt.dat', delim_whitespace=True, header=0)

## Which residue was changed compared to WT per system based
## Use the PDB number, not biology
changed_resnum = [130, 130, 145, 275, 375]

## The name of the system being compared to WT
systems = ['SysA','SysB','SysC','SysD','SysE']

## Define some functions

## https://stackoverflow.com/questions/16490261/python-pandas-write-dataframe-to-fixed-width-file-to-fwf
def write_fwf(df):
    content = tabulate(df.values.tolist(), list(df.columns), tablefmt="plain")
    return content

def filter_345(df):
    saved_colnames = df.columns
    ## Rename the columns for script
    df.columns = ['Acceptor','Donor','AF_1A','AF_1B','AD_1','AF_2A','AF_2B',\
    'AD_2','AF_3A','AF_3B','AD_3','AF_4A','AF_4B','AD_4','AF_5A','AF_5B','AD_5']
    #
    ## Keep those without '-' in first group
    #d2 = df[(df.AF_1A != '-') & (df.AF_1B != '-') & (df.AD_1 != '-')]
    all5 = df[(df.AF_1A != '-') & (df.AF_1B != '-') & (df.AD_1 != '-') \
    & (df.AF_2A != '-') & (df.AF_2B != '-') & (df.AD_2 != '-') \
    & (df.AF_3A != '-') & (df.AF_3B != '-') & (df.AD_3 != '-') \
    & (df.AF_4A != '-') & (df.AF_4B != '-') & (df.AD_4 != '-') \
    & (df.AF_5A != '-') & (df.AF_5B != '-') & (df.AD_5 != '-')]
    ## Reset index
    all5 = all5.reset_index(drop=True)
    #
    ##########################
    ##     Get 4 Groups     ##
    ##########################
    #
    ## | is or , & is and -- pandas gets angry with the actual words
    ## Get the not types, dfs where only 1 rep has - - -
    not1 = df[(df.AF_2A != '-') & (df.AF_2B != '-') & (df.AD_2 != '-') \
    & (df.AF_3A != '-') & (df.AF_3B != '-') & (df.AD_3 != '-') \
    & (df.AF_4A != '-') & (df.AF_4B != '-') & (df.AD_4 != '-') \
    & (df.AF_5A != '-') & (df.AF_5B != '-') & (df.AD_5 != '-')]
    ## Reset index
    not1 = not1.reset_index(drop=True)
    #
    not2 = df[(df.AF_1A != '-') & (df.AF_1B != '-') & (df.AD_1 != '-') \
    & (df.AF_3A != '-') & (df.AF_3B != '-') & (df.AD_3 != '-') \
    & (df.AF_4A != '-') & (df.AF_4B != '-') & (df.AD_4 != '-') \
    & (df.AF_5A != '-') & (df.AF_5B != '-') & (df.AD_5 != '-')]
    ## Reset index
    not2 = not2.reset_index(drop=True)
    #
    not3 = df[(df.AF_1A != '-') & (df.AF_1B != '-') & (df.AD_1 != '-') \
    & (df.AF_2A != '-') & (df.AF_2B != '-') & (df.AD_2 != '-') \
    & (df.AF_4A != '-') & (df.AF_4B != '-') & (df.AD_4 != '-') \
    & (df.AF_5A != '-') & (df.AF_5B != '-') & (df.AD_5 != '-')]
    ## Reset index
    not3 = not3.reset_index(drop=True)
    #
    not4 = df[(df.AF_1A != '-') & (df.AF_1B != '-') & (df.AD_1 != '-') \
    & (df.AF_2A != '-') & (df.AF_2B != '-') & (df.AD_2 != '-') \
    & (df.AF_3A != '-') & (df.AF_3B != '-') & (df.AD_3 != '-') \
    & (df.AF_5A != '-') & (df.AF_5B != '-') & (df.AD_5 != '-')]
    ## Reset index
    not4 = not4.reset_index(drop=True)
    #
    not5 = df[(df.AF_1A != '-') & (df.AF_1B != '-') & (df.AD_1 != '-') \
    & (df.AF_2A != '-') & (df.AF_2B != '-') & (df.AD_2 != '-') \
    & (df.AF_3A != '-') & (df.AF_3B != '-') & (df.AD_3 != '-') \
    & (df.AF_4A != '-') & (df.AF_4B != '-') & (df.AD_4 != '-')]
    ## Reset index
    not5 = not5.reset_index(drop=True)
    #
    ## Append all 4-groups together
    only4match = not1.append(not2,ignore_index=True).append(not3,ignore_index=True)\
    .append(not4,ignore_index=True).append(not5,ignore_index=True)
    #
    ## Drop duplicates and reset the index
    only4match = only4match.drop_duplicates().reset_index(drop=True)
    #
    ## Filter out the 5-matches
    only4 = pd.concat([only4match, all5]).drop_duplicates(keep=False).reset_index(drop=True)
    #
    ##########################
    ##     Get 3 Groups     ##
    ##########################
    #
    ## | is or , & is and -- pandas gets angry with the actual words
    ## Get the not types, dfs where 2 reps have - - -
    not12 = df[(df.AF_3A != '-') & (df.AF_3B != '-') & (df.AD_3 != '-') \
    & (df.AF_4A != '-') & (df.AF_4B != '-') & (df.AD_4 != '-') \
    & (df.AF_5A != '-') & (df.AF_5B != '-') & (df.AD_5 != '-')]
    ## Reset index
    not12 = not12.reset_index(drop=True)
    #
    not13 = df[(df.AF_2A != '-') & (df.AF_2B != '-') & (df.AD_2 != '-') \
    & (df.AF_4A != '-') & (df.AF_4B != '-') & (df.AD_4 != '-') \
    & (df.AF_5A != '-') & (df.AF_5B != '-') & (df.AD_5 != '-')]
    ## Reset index
    not13 = not13.reset_index(drop=True)
    #
    not14 = df[(df.AF_2A != '-') & (df.AF_2B != '-') & (df.AD_2 != '-') \
    & (df.AF_3A != '-') & (df.AF_3B != '-') & (df.AD_3 != '-') \
    & (df.AF_5A != '-') & (df.AF_5B != '-') & (df.AD_5 != '-')]
    ## Reset index
    not14 = not14.reset_index(drop=True)
    #
    not15 = df[(df.AF_2A != '-') & (df.AF_2B != '-') & (df.AD_2 != '-') \
    & (df.AF_3A != '-') & (df.AF_3B != '-') & (df.AD_3 != '-') \
    & (df.AF_4A != '-') & (df.AF_4B != '-') & (df.AD_4 != '-')]
    ## Reset index
    not15 = not15.reset_index(drop=True)
    #
    not23 = df[(df.AF_1A != '-') & (df.AF_1B != '-') & (df.AD_1 != '-') \
    & (df.AF_4A != '-') & (df.AF_4B != '-') & (df.AD_4 != '-') \
    & (df.AF_5A != '-') & (df.AF_5B != '-') & (df.AD_5 != '-')]
    ## Reset index
    not23 = not23.reset_index(drop=True)
    #
    not24 = df[(df.AF_1A != '-') & (df.AF_1B != '-') & (df.AD_1 != '-') \
    & (df.AF_3A != '-') & (df.AF_3B != '-') & (df.AD_3 != '-') \
    & (df.AF_5A != '-') & (df.AF_5B != '-') & (df.AD_5 != '-')]
    ## Reset index
    not24 = not24.reset_index(drop=True)
    #
    not25 = df[(df.AF_1A != '-') & (df.AF_1B != '-') & (df.AD_1 != '-') \
    & (df.AF_3A != '-') & (df.AF_3B != '-') & (df.AD_3 != '-') \
    & (df.AF_4A != '-') & (df.AF_4B != '-') & (df.AD_4 != '-')]
    ## Reset index
    not25 = not25.reset_index(drop=True)
    #
    not34 = df[(df.AF_1A != '-') & (df.AF_1B != '-') & (df.AD_1 != '-') \
    & (df.AF_2A != '-') & (df.AF_2B != '-') & (df.AD_2 != '-') \
    & (df.AF_5A != '-') & (df.AF_5B != '-') & (df.AD_5 != '-')]
    ## Reset index
    not34 = not34.reset_index(drop=True)
    #
    not35 = df[(df.AF_1A != '-') & (df.AF_1B != '-') & (df.AD_1 != '-') \
    & (df.AF_2A != '-') & (df.AF_2B != '-') & (df.AD_2 != '-') \
    & (df.AF_4A != '-') & (df.AF_4B != '-') & (df.AD_4 != '-')]
    ## Reset index
    not35 = not35.reset_index(drop=True)
    #
    not45 = df[(df.AF_1A != '-') & (df.AF_1B != '-') & (df.AD_1 != '-') \
    & (df.AF_2A != '-') & (df.AF_2B != '-') & (df.AD_2 != '-') \
    & (df.AF_3A != '-') & (df.AF_3B != '-') & (df.AD_3 != '-')]
    ## Reset index
    not45 = not45.reset_index(drop=True)
    #
    ## Append all 3-groups together
    only3match = not12.append(not13,ignore_index=True).append(not14,\
        ignore_index=True).append(not15,ignore_index=True).append(not23,\
        ignore_index=True).append(not24,ignore_index=True).append(not25,\
        ignore_index=True).append(not34,ignore_index=True).append(not35,\
        ignore_index=True).append(not45,ignore_index=True)
    #
    ## Drop duplicates and reset the index
    only3match = only3match.drop_duplicates().reset_index(drop=True)
    #
    ## Filter out the 5-matches and 4-matches
    only3 = pd.concat([only3match, all5]).drop_duplicates(keep=False).reset_index(drop=True)
    only3 = pd.concat([only3, only4]).drop_duplicates(keep=False).reset_index(drop=True)
    #
    return saved_colnames, df, only3, only4, all5

def write_filtered(only3, only4, all5, outfile):
    with open(outfile, "w+") as f:
        f.write("Different across 5 systems\n\n")
        f.write(write_fwf(all5))
        f.write("\n\n\nDifferent across 4 systems\n\n")
        f.write(write_fwf(only4))
        f.write("\n\n\nDifferent across 3 systems\n\n")
        f.write(write_fwf(only3))

def chimera_setup(df, conditions):
    ## https://stackoverflow.com/questions/50140131/multiple-logical-comparisons-in-pandas-df
    ## Create a clean df to work from with just Acceptors and Donors
    print_df = pd.concat([df.Acceptor,df.Donor], axis=1)
    ## Keep only residue number -- use \' because nucleic acids are mean
    print_df.Acceptor = print_df.Acceptor.str.split('\'').str[0].str.split('_').str[1].str.split('@').str[0]
    print_df.Donor = print_df.Donor.str.split('\'').str[0].str.split('_').str[1].str.split('@').str[0]
    ## Cycle through for those favored in 1A and those favored in 1B, then
    ## denote that in a new column. If they don't match, label 'NA'
    ## conditions = [(df.AF_1A != '-') & (df.AF_1A > df.AF_1B),
    ##           (df.AF_1A != '-') & (df.AF_1A < df.AF_1B)]
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

def write_chimera(command_file, changed_resnum, systems, only1A, only1B, both1,\
 only2A, only2B, both2,only3A, only3B, both3, only4A, only4B, both4,\
 only5A, only5B, both5):
    original_stdout = sys.stdout
    with open(command_file,'w+') as f:
        sys.stdout = f
        f.write("#--------- System 1: {} ---------#\n".format(systems[0]))
        f.write("Favors Sys 1A (WT)\n")
        print('select',*sorted(only1A), sep=" :")
        f.write("\nFavors Sys 1B (SNP)\n")
        print('select',*sorted(only1B), sep=" :")
        f.write("\nIn Both 1A and 1B\n")
        print('select',*sorted(both1), sep=" :")
        f.write("\nSNP Position\n")
        print('select :{}'.format(changed_resnum[0]))
        # 2
        f.write("\n\n#--------- System 2: {} ---------#\n".format(systems[1]))
        f.write("Favors Sys 2A (WT)\n")
        print('select',*sorted(only2A), sep=" :")
        f.write("\nFavors Sys 2B (SNP)\n")
        print('select',*sorted(only2B), sep=" :")
        f.write("\nIn Both 2A and 2B\n")
        print('select',*sorted(both2), sep=" :")
        f.write("\nSNP Position\n")
        print('select :{}'.format(changed_resnum[1]))
        # 3
        f.write("\n\n#--------- System 3: {} ---------#\n".format(systems[2]))
        f.write("Favors Sys 3A (WT)\n")
        print('select',*sorted(only3A), sep=" :")
        f.write("\nFavors Sys 3B (SNP)\n")
        print('select',*sorted(only3B), sep=" :")
        f.write("\nIn Both 3A and 3B\n")
        print('select',*sorted(both3), sep=" :")
        f.write("\nSNP Position\n")
        print('select :{}'.format(changed_resnum[2]))
        # 4
        f.write("\n\n#--------- System 4: {} ---------#\n".format(systems[3]))
        f.write("Favors Sys 4A (WT)\n")
        print('select',*sorted(only4A), sep=" :")
        f.write("\nFavors Sys 4B (SNP)\n")
        print('select',*sorted(only4B), sep=" :")
        f.write("\nIn Both 4A and 4B\n")
        print('select',*sorted(both4), sep=" :")
        f.write("\nSNP Position\n")
        print('select :{}'.format(changed_resnum[3]))
        # 5
        f.write("\n\n#--------- System 5: {} ---------#\n".format(systems[4]))
        f.write("Favors Sys 5A (WT)\n")
        print('select',*sorted(only5A), sep=" :")
        f.write("\nFavors Sys 5B (SNP)\n")
        print('select',*sorted(only5B), sep=" :")
        f.write("\nIn Both 5A and 5B\n")
        print('select',*sorted(both5), sep=" :")
        f.write("\nSNP Position\n")
        print('select :{}'.format(changed_resnum[4]))
        f.close()
        sys.stdout = original_stdout

## Run the system
saved_colnames, df, only3, only4, all5 = filter_345(df)
write_filtered(only3, only4, all5, outfile)

## For 1
conditions = [(df.AF_1A != '-') & (df.AF_1A > df.AF_1B),
          (df.AF_1A != '-') & (df.AF_1A < df.AF_1B)]
only1A, only1B, both1 = chimera_setup(df, conditions)

## For 2
conditions = [(df.AF_2A != '-') & (df.AF_2A > df.AF_2B),
          (df.AF_2A != '-') & (df.AF_2A < df.AF_2B)]
only2A, only2B, both2 = chimera_setup(df, conditions)

## For 3
conditions = [(df.AF_3A != '-') & (df.AF_3A > df.AF_3B),
          (df.AF_3A != '-') & (df.AF_3A < df.AF_3B)]
only3A, only3B, both3 = chimera_setup(df, conditions)

## For 4
conditions = [(df.AF_4A != '-') & (df.AF_4A > df.AF_4B),
          (df.AF_4A != '-') & (df.AF_4A < df.AF_4B)]
only4A, only4B, both4 = chimera_setup(df, conditions)

## For 5
conditions = [(df.AF_5A != '-') & (df.AF_5A > df.AF_5B),
          (df.AF_5A != '-') & (df.AF_5A < df.AF_5B)]
only5A, only5B, both5 = chimera_setup(df, conditions)

## Make the file for Chimera copy/paste
write_chimera(command_file, changed_resnum, systems, only1A, only1B, both1,\
 only2A, only2B, both2,only3A, only3B, both3, only4A, only4B, both4,\
 only5A, only5B, both5)
