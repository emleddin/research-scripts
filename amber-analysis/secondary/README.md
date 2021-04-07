# Secondary
These scripts (well, this script) pertains to analyzing the secondary structure
of a protein.

## `secstruct-gnuplot-fix.sh`

This script modifies the secondary structure gnuplot script that `cpptraj`
generates.
The unmodified version from `cpptraj` may crash `gnuplot` (depending on the
`cpptraj` release), and prints axis labels for each residue.
While this may not be an issue for a small system, the data from proteins with
more than about 20 residues become unreadable.

You can use this script to process 1+ files by defining filenames (`$gnu_files`)
and image save files (`$out_pngs`) on their own line.

```
gnu_files="just_1_file"

gnu_files="file1
file2
file3
...
fileX"
```

## `2SA-plot-fix.py`
This is a Python version of the script to modify the secondary structure
gnuplot script that `cpptraj` generates.
It is much faster than the bash version.
You will want to remember to modify the `ytics` of the residue range.
