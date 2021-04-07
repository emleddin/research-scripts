import pandas as pd
import math

## Read in data from rmagic-EDA-avg.R
## header = 0 reads header in first row because Python starts at 0
d1 = "MUT-A"

d1C = pd.read_csv('MUT-A/EDA/MUT-protein-system-A-EDA-resXXX-coul-avg.dat', \
 delim_whitespace=True, header=0)

d1V = pd.read_csv('MUT-A/EDA/MUT-protein-system-A-EDA-resXXX-vdw-avg.dat', \
 delim_whitespace=True, header=0)

d1T = pd.read_csv('MUT-A/EDA/MUT-protein-system-A-EDA-resXXX-tot-avg.dat', \
 delim_whitespace=True, header=0)

# D2
d2 = "MUT-B"

d2C = pd.read_csv('MUT-B/EDA/MUT-protein-system-B-EDA-resXXX-coul-avg.dat', \
 delim_whitespace=True, header=0)

d2V = pd.read_csv('MUT-B/EDA/MUT-protein-system-B-EDA-resXXX-vdw-avg.dat', \
 delim_whitespace=True, header=0)

d2T = pd.read_csv('MUT-B/EDA/MUT-protein-system-B-EDA-resXXX-tot-avg.dat', \
 delim_whitespace=True, header=0)

# D3
d3 = "MUT-C"

d3C = pd.read_csv('MUT-C/EDA/MUT-protein-system-C-EDA-resXXX-coul-avg.dat', \
 delim_whitespace=True, header=0)

d3V = pd.read_csv('MUT-C/EDA/MUT-protein-system-C-EDA-resXXX-vdw-avg.dat', \
 delim_whitespace=True, header=0)

d3T = pd.read_csv('MUT-C/EDA/MUT-protein-system-C-EDA-resXXX-tot-avg.dat', \
 delim_whitespace=True, header=0)

# D4
d4 = "WT-A"

d4C = pd.read_csv('WT-A/EDA/WT-protein-system-A-EDA-resXXX-coul-avg.dat', \
 delim_whitespace=True, header=0)

d4V = pd.read_csv('WT-A/EDA/WT-protein-system-A-EDA-resXXX-vdw-avg.dat', \
 delim_whitespace=True, header=0)

d4T = pd.read_csv('WT-A/EDA/WT-protein-system-A-EDA-resXXX-tot-avg.dat', \
 delim_whitespace=True, header=0)

# D5
d5 = "WT-B"

d5C = pd.read_csv('WT-B/EDA/WT-protein-system-B-EDA-resXXX-coul-avg.dat', \
 delim_whitespace=True, header=0)

d5V = pd.read_csv('WT-B/EDA/WT-protein-system-B-EDA-resXXX-vdw-avg.dat', \
 delim_whitespace=True, header=0)

d5T = pd.read_csv('WT-B/EDA/WT-protein-system-B-EDA-resXXX-tot-avg.dat', \
 delim_whitespace=True, header=0)

# D6
d6 = "WT-C"

d6C = pd.read_csv('WT-C/EDA/WT-protein-system-C-EDA-resXXX-coul-avg.dat', \
 delim_whitespace=True, header=0)

d6V = pd.read_csv('WT-C/EDA/WT-protein-system-C-EDA-resXXX-vdw-avg.dat', \
 delim_whitespace=True, header=0)

d6T = pd.read_csv('WT-C/EDA/WT-protein-system-C-EDA-resXXX-tot-avg.dat', \
 delim_whitespace=True, header=0)

#------ Calculate
coul_d1 = d1C.loc[:,"AvgCoulomb"].sum()
coulsd_d1 = math.sqrt(d1C.loc[:, "AvgCoulombSD"].var())

vdw_d1 = d1V.loc[:,"AvgVdW"].sum()
vdwsd_d1 = math.sqrt(d1V.loc[:, "AvgVdWSD"].var())

tot_d1 = d1T.loc[:,"AvgIntTot"].sum()
totsd_d1 = math.sqrt(d1T.loc[:, "AvgStdDev"].var())

# D2
coul_d2 = d2C.loc[:,"AvgCoulomb"].sum()
coulsd_d2 = math.sqrt(d2C.loc[:, "AvgCoulombSD"].var())

vdw_d2 = d2V.loc[:,"AvgVdW"].sum()
vdwsd_d2 = math.sqrt(d2V.loc[:, "AvgVdWSD"].var())

tot_d2 = d2T.loc[:,"AvgIntTot"].sum()
totsd_d2 = math.sqrt(d2T.loc[:, "AvgStdDev"].var())

# D3
coul_d3 = d3C.loc[:,"AvgCoulomb"].sum()
coulsd_d3 = math.sqrt(d3C.loc[:, "AvgCoulombSD"].var())

vdw_d3 = d3V.loc[:,"AvgVdW"].sum()
vdwsd_d3 = math.sqrt(d3V.loc[:, "AvgVdWSD"].var())

tot_d3 = d3T.loc[:,"AvgIntTot"].sum()
totsd_d3 = math.sqrt(d3T.loc[:, "AvgStdDev"].var())

# D4
coul_d4 = d4C.loc[:,"AvgCoulomb"].sum()
coulsd_d4 = math.sqrt(d4C.loc[:, "AvgCoulombSD"].var())

vdw_d4 = d4V.loc[:,"AvgVdW"].sum()
vdwsd_d4 = math.sqrt(d4V.loc[:, "AvgVdWSD"].var())

tot_d4 = d4T.loc[:,"AvgIntTot"].sum()
totsd_d4 = math.sqrt(d4T.loc[:, "AvgStdDev"].var())

# D5
coul_d5 = d5C.loc[:,"AvgCoulomb"].sum()
coulsd_d5 = math.sqrt(d5C.loc[:, "AvgCoulombSD"].var())

vdw_d5 = d5V.loc[:,"AvgVdW"].sum()
vdwsd_d5 = math.sqrt(d5V.loc[:, "AvgVdWSD"].var())

tot_d5 = d5T.loc[:,"AvgIntTot"].sum()
totsd_d5 = math.sqrt(d5T.loc[:, "AvgStdDev"].var())

# D6
coul_d6 = d6C.loc[:,"AvgCoulomb"].sum()
coulsd_d6 = math.sqrt(d6C.loc[:, "AvgCoulombSD"].var())

vdw_d6 = d6V.loc[:,"AvgVdW"].sum()
vdwsd_d6 = math.sqrt(d6V.loc[:, "AvgVdWSD"].var())

tot_d6 = d6T.loc[:,"AvgIntTot"].sum()
totsd_d6 = math.sqrt(d6T.loc[:, "AvgStdDev"].var())

#------ Print
with open("sum-report-XXX.txt", "w+") as f:
    f.write("                   All Protein\n")
    f.write(f"{d1}: Coul {coul_d1:.2f} ± {coulsd_d1:.2f}\n")
    f.write(f"{d1}: VdW {vdw_d1:.2f} ± {vdwsd_d1:.2f}\n")
    f.write(f"{d1}: tot {tot_d1:.2f} ± {totsd_d1:.2f}\n")
    f.write("\n")

    f.write(f"{d2}: Coul {coul_d2:.2f} ± {coulsd_d2:.2f}\n")
    f.write(f"{d2}: VdW {vdw_d2:.2f} ± {vdwsd_d2:.2f}\n")
    f.write(f"{d2}: tot {tot_d2:.2f} ± {totsd_d2:.2f}\n")
    f.write("\n")

    f.write(f"{d3}: Coul {coul_d3:.2f} ± {coulsd_d3:.2f}\n")
    f.write(f"{d3}: VdW {vdw_d3:.2f} ± {vdwsd_d3:.2f}\n")
    f.write(f"{d3}: tot {tot_d3:.2f} ± {totsd_d3:.2f}\n")
    f.write("\n")

    f.write(f"{d4}: Coul {coul_d4:.2f} ± {coulsd_d4:.2f}\n")
    f.write(f"{d4}: VdW {vdw_d4:.2f} ± {vdwsd_d4:.2f}\n")
    f.write(f"{d4}: tot {tot_d4:.2f} ± {totsd_d4:.2f}\n")
    f.write("\n")

    f.write(f"{d5}: Coul {coul_d5:.2f} ± {coulsd_d5:.2f}\n")
    f.write(f"{d5}: VdW {vdw_d5:.2f} ± {vdwsd_d5:.2f}\n")
    f.write(f"{d5}: tot {tot_d5:.2f} ± {totsd_d5:.2f}\n")
    f.write("\n")

    f.write(f"{d6}: Coul {coul_d6:.2f} ± {coulsd_d6:.2f}\n")
    f.write(f"{d6}: VdW {vdw_d6:.2f} ± {vdwsd_d6:.2f}\n")
    f.write(f"{d6}: tot {tot_d6:.2f} ± {totsd_d6:.2f}\n")
    f.write("\n")
