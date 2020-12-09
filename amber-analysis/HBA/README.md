# Hydrogen Bond Analysis (HBA)

These are a collection of scripts for processing hydrogen bond data.

Each script assumes that `../cpptraj_2.sh` was run with `cpptraj_analysis.in`,
as this generates a specific `.dat` file for each replicate/system.

The `cpptraj` data needs to go through some additional processing in order to
be particularly helpful. Why?
* Multiple lines that have the same hydrogen bond for different lengths of time, which need to be summed.
* Difficult to average across simulations, since they may or may not have the exact same bonds.

## Bash Way (non-ideal)
These scripts use `bash` and `awk` to parse through the hbond data from `cpptraj`.
These are non-ideal because if a match isn't found between two files, it prints
the information to a no-match list, instead of having a file with any hydrogen
bond not meeting the specified cutoff.

### `hbond-analysis-prep.sh`
This averages 3 replicates together, so that they can be compared with another
system using `hbond-analysis.sh`.

### `hbond-analysis.sh`
This script works for
- Comparing single trials of 2 systems/replicates
- Comparing 2 systems previously averaged using `hbond-analysis-prep.sh`

## R way (Good!)
These script use `R` to get around some of the issues with the `bash` approach.
The scripts rely on the `data.table` and `tidyverse` packages to run.

### `rmagic-hbond-avg.r`
This script calls the data files generated with `cpptraj`.
It uses `fix_string` to find a residue number that changed between the two
systems you eventually want to compare.
For instance, if you have SNP I123L, you need to set `fix_string=_123` for both
your wild type and mutant data file processing, so that residue can be compared
as either I123 (ILE) or I123L (LEU).

### `rmagic-hbond-avg-2res.r`
When comparing distinct SNP systems, you need to go through the residue name
change process for both residues.
This script builds on `rmagic-hbond-avg.r` by adding in a second string to find
and replace.
The residue numbers are changed through `fix_stringA` and `fix_stringB`.
The name they are changed to is set with `fixed_stringA` and `fixed_stringB`.

### `system-hbond-table.r`
This script uses the files generated with either `rmagic-hbond-avg.r` or
`rmagic-hbond-avg-2res.r` to create a table of relevant hydrogen bonds, based
on a specific cutoff.

As written, the script is written to perform up to 5 comparisons (3 by default).
Each comparison is matched:
  - infile1a and infile1b --> A_diff
  - infile2a and infile2b --> B_diff
  - infile3a and infile3b --> C_diff
  - infile4a and infile4b --> D_diff
  - infile5a and infile5b --> E_diff

The number of comparisons is set through `sets`, which must be a float.
The tags (e.g., `tag1a`, `tag1b`, and `tag1c`) set the column names for the
output files.
`cutoff` sets the percentage of total simulation time to compare between.
While it is based off total time, that time will be system dependent.
So, if you're comparing A-system that ran for 20 ns and B-system that ran for
200 ns, a 4 ns presence in A-system would be equivalent to a 40 ns presence in
B-system.

### `df-to-chimera.py`
This script uses the file generated from `system-hbond-table.r` to print
filtered tables that show an Acceptor/Donor pair as important in 3, 4, or 5
systems.
Independent of these filtered groups, it will also print a file of Chimera
command line commands for selecting residues across each WT/MUT pair that are
- (a) present longer in the WT
- (b) present for longer in the MUT system
- (c) multiple bonds present, some for longer in WT, others longer in MUT, and
- (d) the MUT position itself.

This will miss bonds that were marked as "No Match", but there's typically
only 2-3 per system, so it's not as much of a hassle to add them by hand.
