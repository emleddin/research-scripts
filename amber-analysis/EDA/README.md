# Energy Decomposition Analysis (EDA)
These are a number of scripts for performing an energy decomposition analysis
for AMBER MD trajectories. They use the `Residue_E_Decomp_07_15.f90` program.
An OpenMP version of this program can be found on the
[CisnerosRes GitHub](https://github.com/CisnerosResearch/AMBER-EDA).

## Running the FORTRAN Program
The program will have 3 output files: `fort.803`, `fort.804`, and `fort.806`.
`fort.803` contains the Coulomb energies, `fort.804` contains a sanity check,
and `fort.806` contains the van der Waals energies.
If you're using the OpenMP version, the naming will be
`fort_coulomb_interaction.dat` for the Coulomb energies,
`fort_sanity_check.txt` for the sanity check, and
`fort_vdw_interaction.dat` for the van der Waals energies.

A critical note for the FORTRAN program is that it will not work on NetCDF files,
only the ASCII mdcrd.
Thus, they must be converted with `cpptraj` to mdcrd.
While converting, do not use autoimage or stripping, as these will affect the
results and give absurdly high energies.
If they're originally written to an mdcrd, you're good to go.

### `EDA_script.sh`
This is a PBS run script.

### `ans.txt`
This contains the answers to the FORTRAN program prompts so it can run through
a queuing system.

### `EDA_new.inp`
The actual input file for the program that contains information on the number
of atoms and frames.

## Processing the FORTRAN Data
These scripts are used to process through the `fort.803` and `fort.806` data
(from replicate runs).

### `rmagic-EDA-avg.r`
The script relies on the `data.table` and `tidyverse` packages to run.
The `ROI` (residue of interest) selects every other residue with that residue.
If a SNP is being investigated, the `ROI` should likely be the SNP position.
Because of the inability to truly separated bonded and non-bonded interactions,
the two residues directly next to the `ROI` are set to zero.

### `rmagic-EDA-avg-diffs.r`
This script finds the difference between two files of averages obtained through
`rmagic-EDA-avg.r`. The `X_val` variable should be the same as the `ROI` in
that script.

### `rmagic-EDA-single-run.r`
The script relies on the `data.table` and `tidyverse` packages to run.
The `ROI` (residue of interest) selects every other residue with that residue.
If a SNP is being investigated, the `ROI` should likely be the SNP position.
Because of the inability to truly separated bonded and non-bonded interactions,
the two residues directly next to the `ROI` are set to zero.
This script is for a single run (so only 1 set of `fort.803` and `fort.806`
files, NOT for replicates that need to be averaged).

### `chimeraprint.sh`
This script can be used to map differences in the Coulomb and van der Waals
energies onto the protein structure using Chimera.
It calls the file created using `rmagic-EDA-avg-diffs.r`.

## Plotting the Processed Data

### `EDA-bar-comp.py`
This script works with the file generated from `rmagic-EDA-avg.r` or
`rmagic-EDA-single-run.r` to plot 2 systems side-by-side in a barplot
(i.e., the x-axis would look like `A1 B1 A2 B2` for systems A and B.)

### `EDA-diffs.py`
This script will plot barplots of 5 sets of EDA difference data.
The difference data files can be made by running `rmagic-EDA-avg-diffs.r`
multiple times [i.e., (1) WT-MutA, (2) WT-MutB, (3) WT-MutC, (4) WT-MutD, and
(5) WT-MutD].
The resulting plot has one x-axis for residues and individual y-axes for
each dataset.

### `EDA-single-diff.py`
This script will create a barplot of 1 set of EDA difference data.
The difference data file can be made by running `rmagic-EDA-avg-diffs.r`.

### `EDA-candlestick.gnu`
This script will create a box and whisker plot of 1 set of EDA difference data
using gnuplot.
The difference data file can be made by running `rmagic-EDA-avg-diffs.r`.
For the best auto-arrangement, first save the image as an EPS file and then
convert it to PNG.
Conversion can be done with ImageMagick (`-flatten` gives it a white background
instead of transparent).
```
convert -flatten image.eps -rotate 90 image.png
```

### Rendered Examples
This directory has examples of the Python and Gnuplot plots made from the above
scripts. Toy data were created through commands like:
```python
import pandas as pd
import numpy as np

## Create a sample data frame with 450 values
d1 = pd.DataFrame(np.arange(1,451))
## Name that column "Residues"
d1.columns = ['Residues']
## Add a normalized (Gaussian) dataset centered on 0 with StDev 0.5
d1['DiffE'] = np.random.normal(loc=0,scale=0.5,size=450)
## Create a uniform dataset for StDev from 0 to 0.1
d1['AvgSTDEV'] = np.random.uniform(low=0., high=0.1, size=450)

## Save a file for Gnuplot examples
with open("toy_data.dat", 'w+') as f:
    f.write("  Residue   TotIntE    AvgStDev\n")
    np.savetxt(f, d1.values, fmt='%10.5f')
```

A secondary way to create toy data, particularly with large outliers, is
based off a
[StackExchange answer](https://stackoverflow.com/questions/55351782/how-should-i-generate-outliers-randomly).
```python
import pandas as pd
import numpy as np

def generate(median=0, std_dev=0.5, max_err=38, size=420, outlier_size=15):
    errs = np.random.normal(loc=median, scale=std_dev, size=size)
    data = median + errs
    #
    #
    lower_errs = np.random.randint(1, max_err + 1) * \
     np.random.rand(outlier_size)
    lower_outliers = median - lower_errs
    #
    for i in range(len(lower_outliers)):
        if lower_outliers[i] < -19. and lower_outliers[i] > -40.:
            lower_outliers[i] += -20.
    #
    upper_errs = np.random.randint(1, max_err + 1) * \
     np.random.rand(outlier_size)
    upper_outliers = median + upper_errs
    #
    data = np.concatenate((data, lower_outliers, upper_outliers))
    np.random.shuffle(data)
    #
    return data

data = generate()

d1 = pd.DataFrame(np.arange(1,451))
## Name that column "Residues"
d1.columns = ['Residues']
## Add generated data with outliers
d1['DiffE'] = generate()
## Create a uniform dataset for StDev from 0 to 1
d1['AvgSTDEV'] = np.random.uniform(low=0., high=1, size=450)
```
