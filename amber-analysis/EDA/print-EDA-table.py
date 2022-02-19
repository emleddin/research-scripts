#!/usr/env/python3
import numpy as np
import pandas as pd

## Threshold value (kcal/mol) -- only values above +thresh and
##  below -thresh will written out. Default: 1.
thresh=1

## Read in data
## header = 0 reads header in first row because Python starts at 0
## Note: if you want outputs in multiple formats, read in the dataset multiple
##  times. Ex: df1 = pd.read_csv('file1'); df2 = pd.read_csv('file1')
df1 = pd.read_csv('WT-MUTA_total_interaction_res220_avg.dat', \
 delim_whitespace=True, header=0)
df2 = pd.read_csv('WT-MUTB_total_interaction_res220_avg.dat', \
 delim_whitespace=True, header=0)
df3 = pd.read_csv('WT-MUTC_total_interaction_res220_avg.dat', \
 delim_whitespace=True, header=0)
df4 = pd.read_csv('WT-MUTD_total_interaction_res220_avg.dat', \
delim_whitespace=True, header=0)

## Are you providing a map file for the residues? It expects output from
##  research-scripts/amber-setup/system-preparation/pdb-numbering.py
## True/False required!
map_values = True

## Protein map file created with
##  research-scripts/amber-setup/system-preparation/pdb-numbering.py
map = pd.read_csv("PDB_numbering.txt", delim_whitespace=True, header=1)

#------------------#
# Setting Up Files #
#------------------#
## Define a list of tuples with (data, outfile, output_type) for the output
##  files as `datasets`.
##  Options for outut_type: 'latex', 'latex-siunitx', 'csv', or 'txt'
##    'latex' writes out table row format. Ex: A123 1.02 \pm 0.22 \\
##    'latex-siunitx' writes out numbers formatted with siunitx
##     Ex: A123 \num{1.02} \pm \num{0.22} \\
datasets = [
    (df1, "resthresh_A.tex", 'latex'),
    (df2, "resthresh_B.csv", 'latex-siunitx'),
    (df3, "resthresh_C.csv", 'csv'),
    (df4, "resthresh_D.txt", 'txt'),
]

##---------------------------------------------------------------------------##
##---------------------------- Behind the Curtain ---------------------------##
##---------------------------------------------------------------------------##

def get_threshold_values(thresh, df, map_values, outfile=None, map=None):
    """Select values above +thresh and below -thresh to write out.
    If a map of biological protein numbering is a available, that information
    will be saved instead of the exisiting numbering.
    """
    thresh=np.positive(thresh)
    #
    ## Define the column names to use throughout
    ## This renames the columns in case something was weird
    df.columns = ['Residue', 'DiffE', 'AvgSTDEV']
    #
    if map_values == True:
        ## Clean up map data to remove any blank lines in string columns
        ## Start by making sure a map was given
        if isinstance(map, type(type)):
            print("\nmap_values=True, but no map given.")
            print(f"Using residues instead for {outfile}.")
            df["Map"] = df["Residue"]
        else:
            map['LEAP_ResName'].replace('', np.nan, inplace=True)
            map['RCSB_ResName'].replace('', np.nan, inplace=True)
            map['1LC'].replace('', np.nan, inplace=True)
            map.dropna(inplace=True)
            ## Add protein map to EDA data
            df["Map"] = map["RCSB_ResName"]+map["RCSB_ResID"].astype(int).astype(str)
    #
    ## Print the rows > thresh
    gt_rows = df[df.DiffE >= thresh]
    lt_rows = df[df.DiffE <= -thresh]
    #
    ## Create a complete dataframe of those exceeding threshold
    thresh_rows = pd.concat([gt_rows, lt_rows])
    thresh_rows.sort_values(by = 'Residue', inplace=True)
    return thresh_rows

def write_out_latex(outfile, thresh_rows, map_values):
    """Saves a TEX file with the Residue, DiffE, and AvgSTDEV as rows of a
    LaTeX table.
    """
    with open(outfile, "w+") as f:
        if map_values == True:
            ## Use {{ and }} as the escape for a single brace
            for c in thresh_rows.itertuples(index=True, name='Pandas'):
                f.write("{:<8} & {:>6.2f} \pm {:>6.2f} \\\\\n".format(\
                c.Map, c.DiffE, c.AvgSTDEV))
            f.close()
        else:
            for c in thresh_rows.itertuples(index=True, name='Pandas'):
                f.write("{:<8} & {:>6.2f} \pm {:>6.2f} \\\\\n".format(\
                c.Residue, c.DiffE, c.AvgSTDEV))
            f.close()

def write_out_latex_siunitx(outfile, thresh_rows, map_values):
    """Saves a TEX file with the Residue, DiffE, and AvgSTDEV as rows of a
    LaTeX table.
    """
    with open(outfile, "w+") as f:
        if map_values == True:
            ## Use {{ and }} as the escape for a single brace
            for c in thresh_rows.itertuples(index=True, name='Pandas'):
                f.write("{:<8} & \\num{{{:>6.2f}  \pm {:>6.2f}}} \\\\\n".format(\
                c.Map, c.DiffE, c.AvgSTDEV))
            f.close()
        else:
            for c in thresh_rows.itertuples(index=True, name='Pandas'):
                f.write("{:<8} & \\num{{{:>6.2f}  \pm {:>6.2f}}}\\\\\n".format(\
                c.Residue, c.DiffE, c.AvgSTDEV))
            f.close()

def write_out_csv(outfile, thresh_rows, map_values):
    """Saves a CSV file with the Residue, DiffE, and AvgSTDEV.
    """
    with open(outfile, "w+") as f:
        if map_values == True:
            for c in thresh_rows.itertuples(index=True, name='Pandas'):
                f.write("{},{:.2f},{:.2f}\n".format(\
                c.Map, c.DiffE, c.AvgSTDEV))
            f.close()
        else:
            for c in thresh_rows.itertuples(index=True, name='Pandas'):
                f.write("{},{:.2f},{:.2f}\n".format(\
                c.Residue, c.DiffE, c.AvgSTDEV))
            f.close()

def write_out_txt(outfile, thresh_rows, map_values):
    """Saves a TXT file with the Residue, DiffE, and AvgSTDEV.
    """
    with open(outfile, "w+") as f:
        if map_values == True:
            for c in thresh_rows.itertuples(index=True, name='Pandas'):
                f.write("{:<8} {:>6.2f} {:>6.2f}\n".format(\
                c.Map, c.DiffE, c.AvgSTDEV))
            f.close()
        else:
            for c in thresh_rows.itertuples(index=True, name='Pandas'):
                f.write("{:<8} {:>6.2f} {:>6.2f}\n".format(\
                c.Residue, c.DiffE, c.AvgSTDEV))
            f.close()

def run_all(df, outfile, output_type, map_values, thresh=1., map=None):
    """
    Get the threshold values and write out the requested file formats.
    """
    thresh_rows = get_threshold_values(thresh, df, map_values, outfile, map)
    #
    ## Write out the appropriate file
    if  output_type.lower() == 'latex':
        write_out_latex(outfile, thresh_rows, map_values)
    elif  output_type.lower() == 'latex-siunitx':
        print("\nMake sure to add \\usepackage{siunitx} to your LaTeX preamble!\n")
        write_out_latex_siunitx(outfile, thresh_rows, map_values)
    elif output_type.lower() == 'csv':
        write_out_csv(outfile, thresh_rows, map_values)
    elif output_type.lower() == 'txt':
        write_out_txt(outfile, thresh_rows, map_values)
    else:
        print(f"I can't write {output_type}.")
        print("Valid options: 'latex', 'latex-siunitx', 'csv', & 'txt'.")
        print(f"Writing a txt file as {outfile}.txt.")
        write_out_txt(outfile+".txt", thresh_rows, map_values)

#---------------------#
# Create Output Files #
#---------------------#
for df, outfile, output_type in datasets:
    thresh = float(thresh)
    run_all(df, outfile, output_type, map_values, thresh, map)
