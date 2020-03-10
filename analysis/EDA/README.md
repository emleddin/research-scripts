# Energy Decomposition Analysis (EDA)
These are a number of scripts for performing an energy decomposition analysis
for AMBER MD trajectories. They use the `Residue_E_Decomp_07_15.f90` program.

## Running the FORTRAN Program
The program will have 3 output files: `fort.803`, `fort.804`, and `fort.806`.
`fort.803` contains the Coulomb energies, `fort.804` contains the van der Waals
energies, and `fort.806` contains a sanity check.

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
These scripts are used to process through the `fort.803` and `fort.804` data
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

### `chimeraprint.sh`
This script can be used to map differences in the Coulomb and van der Waals
energies onto the protein structure using Chimera.
It calls the file created using `rmagic-EDA-avg-diffs.r`.
