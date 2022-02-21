#!/usr/env/python3
import numpy as np
import pandas as pd

## Input Coulomb
df1Ac = pd.read_csv('../WT/r1/fort_coulomb_interaction.dat', \
 delim_whitespace=True, header=1)
df2Ac = pd.read_csv('../WT/r2/fort_coulomb_interaction.dat', \
 delim_whitespace=True, header=1)
df3Ac = pd.read_csv('../WT/r3/fort_coulomb_interaction.dat', \
 delim_whitespace=True, header=1)

## Input VdW
df1Av = pd.read_csv('../WT/r1/fort_vdw_interaction.dat', \
 delim_whitespace=True, header=0)
df2Av = pd.read_csv('../WT/r2/fort_vdw_interaction.dat', \
 delim_whitespace=True, header=0)
df3Av = pd.read_csv('../WT/r3/fort_vdw_interaction.dat', \
 delim_whitespace=True, header=0)

## Outfile
outfile = "WT-protein-DNA-int.txt"

## Do you want to look for a residue of interest?
## If not, the script will find the protein:DNA interaction
ROI_vers = False

## Set this if ROI_vers = True
## ROI is the residue of interest to get interactions with reference to.
# ROI = 123

## Set these if ROI_vers = False
## protein residues: [first,last]
protein = [1, 300]
## DNA residues: [first,last]
DNA = [301, 310]
## A tag for the system
sys = 'WT'

##--------------------------- Behind the Curtain ---------------------------##

## Find the protein-DNA interactions
def find_matches(df):
    """Find all interactions for the protein-DNA pairs.
    This assumes residues are indexed continuously with all protein residues
    in the first chunk.
    """
    global DNA
    ## Remember how the program writes out numbers and only write 1 case
    ## Aka 1 & 2 never 2 & 1. You can verify by checking that your line count =
    ##  len(range(protein[0], protein[1]+1)) * len(range(DNA[0], DNA[1]+1)) +1
    save_df = df.loc[(df['ResA']<DNA[0]) & (df['ResB']>=DNA[0])]
    ## Reset index and drop copy of old
    save_df.reset_index(inplace=True, drop=True)
    return save_df

def ROI_matches(df, ROI):
    """Find all interactions against a given residue of interest.
    """
    ## Remember how the program writes out numbers and only write 1 case
    save_df = df.loc[(df['ResA']==ROI) | (df['ResB']==ROI)]
    ## Reset index and drop copy of old
    save_df.reset_index(inplace=True, drop=True)
    return save_df

def average_coul(df1, df2, df3):
    """Average three sets of Coulomb data."""
    df = pd.concat([df1, df2, df3])
    df_coul = df.groupby(by=['ResA','ResB'], as_index=False, sort=False).mean()
    df_coul = df_coul.rename(columns={"Index": "Index", "ResA": "ResA", \
     "ResB": "ResB", "CoulEnergy": "AvgCoulomb", "StdErr": "AvgCoulSD"})
    return df_coul

def average_vdw(df1, df2, df3):
    """Average three sets of vdW data."""
    df = pd.concat([df1, df2, df3])
    df_vdw = df.groupby(by=['ResA','ResB'], as_index=False, sort=False).mean()
    df_vdw = df_vdw.rename(columns={"Index": "Index", "ResA": "ResA", \
     "ResB": "ResB", "VdWEnergy": "AvgVdW", "StdErr": "AvgVdWSD"})
    return df_vdw

def average_tot(df_coul, df_vdw, remove_ROI=False, ROI=None):
    """Sum the averaged Coulomb and vdW data. Forcibly set values on either
    side of the residue of interest to 0, since their bond and angle terms
    cannot be separated from Coul/vdW."""
    df_coul = df_coul.rename(columns={"Index": "Index", "ResA": "ResA", \
     "ResB": "ResB", "AvgCoulomb": "AvgIntTot", "AvgCoulSD": "AvgSD"})
    df_vdw = df_vdw.rename(columns={"Index": "Index", "ResA": "ResA", \
     "ResB": "ResB", "AvgVdW": "AvgIntTot", "AvgVdWSD": "AvgSD"})
    df = pd.concat([df_coul, df_vdw])
    df_tot = df.groupby(by=['ResA','ResB'],as_index=False, sort=False).sum()
    if remove_ROI == True:
        float(ROI)
        ## Force set values to zero
        df_tot['AvgIntTot'].mask(df_tot['ResA'] == ROI+1, 0, inplace=True)
        df_tot['AvgIntTot'].mask(df_tot['ResA'] == ROI-1, 0, inplace=True)
        df_tot['AvgIntTot'].mask(df_tot['ResB'] == ROI+1, 0, inplace=True)
        df_tot['AvgIntTot'].mask(df_tot['ResB'] == ROI-1, 0, inplace=True)
        df_tot['AvgIntTot'].mask(df_tot['ResA'] == ROI+1, 0, inplace=True)
        df_tot['AvgIntTot'].mask(df_tot['ResA'] == ROI-1, 0, inplace=True)
        df_tot['AvgIntTot'].mask(df_tot['ResB'] == ROI+1, 0, inplace=True)
        df_tot['AvgIntTot'].mask(df_tot['ResB'] == ROI-1, 0, inplace=True)
        #
        df_tot['AvgSD'].mask(df_tot['ResA'] == ROI+1, 0, inplace=True)
        df_tot['AvgSD'].mask(df_tot['ResA'] == ROI-1, 0, inplace=True)
        df_tot['AvgSD'].mask(df_tot['ResB'] == ROI+1, 0, inplace=True)
        df_tot['AvgSD'].mask(df_tot['ResB'] == ROI-1, 0, inplace=True)
        df_tot['AvgSD'].mask(df_tot['ResA'] == ROI+1, 0, inplace=True)
        df_tot['AvgSD'].mask(df_tot['ResA'] == ROI-1, 0, inplace=True)
        df_tot['AvgSD'].mask(df_tot['ResB'] == ROI+1, 0, inplace=True)
        df_tot['AvgSD'].mask(df_tot['ResB'] == ROI-1, 0, inplace=True)
    return df_tot

def write_tot(df_tot, outfile):
    """Write out the matched list of residues."""
    with open(outfile, "w+") as f:
        f.write("ResA  ResB  AvgIntTot  AvgSD\n")
        for c in df_tot.itertuples(index=True, name='Pandas'):
            f.write("{:>4}  {:>4}  {:>9.2f}  {:>5.2f}\n".format(\
                c.ResA, c.ResB, c.AvgIntTot, c.AvgSD))
        f.close()

def ROI_totals(df1Ac, df1Av, df2Ac, df2Av, df3Ac, df3Av, ROI, outfile):
    """Search for all pairs with a residue of interest, then average the
    Coulomb and vdW interactions. Sum the total, and write it to an output file.
    """
    r1c = ROI_matches(df1Ac, ROI)
    r2c = ROI_matches(df2Ac, ROI)
    r3c = ROI_matches(df3Ac, ROI)
    df_coul = average_coul(r1c, r2c, r3c)

    r1v = ROI_matches(df1Av, ROI)
    r2v = ROI_matches(df2Av, ROI)
    r3v = ROI_matches(df3Av, ROI)
    df_vdw = average_vdw(r1v, r2v, r3v)

    remove_ROI = True
    df_tot = average_tot(df_coul, df_vdw, remove_ROI, ROI)
    write_tot(df_tot, outfile)
    return df_tot

def protein_int_totals(df1Ac, df1Av, df2Ac, df2Av, df3Ac, df3Av, DNA, sys, outfile):
    """Search for all protein-DNA pairs, then average the Coulomb and vdW
    interactions. Sum the total, and write it to an output file.
    Print out a sum total of the protein:DNA interactions.
    """
    r1c = find_matches(df1Ac)
    r2c = find_matches(df2Ac)
    r3c = find_matches(df3Ac)
    df_coul = average_coul(r1c, r2c, r3c)

    r1v = find_matches(df1Av)
    r2v = find_matches(df2Av)
    r3v = find_matches(df3Av)
    df_vdw = average_vdw(r1v, r2v, r3v)

    remove_ROI = False
    df_tot = average_tot(df_coul, df_vdw, remove_ROI)
    write_tot(df_tot, outfile)

    prot_sum = df_tot["AvgIntTot"].sum()
    print(f"The protein:DNA interaction for {sys:<6} is: {prot_sum:.2f} kcal/mol")
    return df_tot

### Run the script
if ROI_vers == True:
    df_tot = ROI_totals(df1Ac, df1Av, df2Ac, df2Av, df3Ac, df3Av, ROI, outfile)
else:
    df_tot = protein_int_totals(df1Ac, df1Av, df2Ac, df2Av, df3Ac, df3Av, DNA, sys, outfile)
