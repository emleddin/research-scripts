# TINKER Params

## param_files
This directory contains combined leaprc files and their dependencies.

## `generate_TINKER_parameters.py`
This script uses ParmEd to read through AMBER FF files and write them in the
TINKER FF format.

### User Notes
- Only one leaprc file can be read at a time, and they cannot be added together.
Thus, you'll need a "summary" leaprc file with all the relevant information.
- The `solvents.lib` file that ships with AMBER has **major issues** with
the ParmEd parser, specifically because of the `BOX` types within it. That's
why there's a `param_files/solvents_hack.lib` file, which doesn't contain the
`BOX` types, and thus won't print hundreds of warnings. The warnings aren't
the problem, but it can't create the residue template, and the `atom` line
generation will throw an error. Fixing the `solvents.lib` connections would
then create thousands of lines that are not pertinent to the simulations I
personally run, so I don't have the incentive to fix it. :smiley:

### Bugs/Things to Address
- TINKER requires connectivity information, and without it, they're interpreted
as zeroes. So, you know, valence issues for everything!
- Need to add a way to loop through the torsions with Xs and write them out.
Concern with a straight loop is that it will overwrite important parameters.
Starting with the X loops, however, should work, because TINKER overwrites
new information.
- The ions hate me. Figure that stuff out.

### Things to Test
- Check mol2/frcmod read/write
- Do I actually need biotypes?!?
