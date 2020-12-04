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
