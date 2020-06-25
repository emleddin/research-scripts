# TINKER Params

## param_files
This directory contains combined leaprc files and their dependencies.

## `generate_TINKER_parameters.py`
This script uses ParmEd to read through prmtop files or AMBER force field
files and write them in the TINKER force field format.

### User Notes

- Make sure that you test any parameters using TINKER `analyze`!!!!! REALLY!!!!
Do NOT forget this part!!!!
- A `pdbxyz`-like script that works directly with parameter files generated
through this method exists (`pdbxyz4amber-pmd-params.py`).
You can access it (and the reverse `xyzpdb` script, `xyzpdb4amber-pmd-params.py`)
[here](https://github.com/emleddin/pdbxyz-xyzpdb).
- This script assumes that `pandas v1.0` is installed.
Previous versions do not have the `ignore_index` functionality for
`drop_duplicates`, which *should* only be an issue when working from a prmtop.

#### Building from a `prmtop` file
- This is the **recommended** approach for building a parameter file. It will
build the parameters based on everything in the `prmtop` file corresponding
to your system. Use an unstripped prmtop file so that the solvent and ion
parameters are built.
- This approach seems to drop an angle in the solvent.
If you're using a force field with TIP3P water built from `parm10.dat`,
you're likely missing
`angle        HW   OW   HW     100.00     104.52` (where `HW` and `OW`
correspond to their atom types in your parameter file). This line can go
anywhere in your new `.prm` file -- the sectioning is more for organization.

#### Building from a `leaprc` file
- Only one `leaprc` file can be read at a time, and they cannot be added together.
Thus, you'll need a "summary" leaprc file with all the relevant information.
- The `solvents.lib` file that ships with AMBER has **major issues** with
the ParmEd parser, specifically because of the `BOX` types within it. That's
why there's a `param_files/solvents_hack.lib` file, which doesn't contain the
`BOX` types, and thus won't print hundreds of warnings. The warnings aren't
the problem, but it can't create the residue template, and the `atom` line
generation will throw an error. Fixing the `solvents.lib` connections would
then create thousands of lines that are not pertinent to the simulations I
personally run, so I don't have the incentive to fix it. :smiley:
- This approach will build the residues defined in the `leaprc` file, so
any non-standard residues specific to your system will need to be incorporated
into the leaprc.
- There's an issue with getting the `atom` and `charge` lines for ions.
It builds an atom type number for them, but they don't appear as
possible residues (probably because they don't have a template).

#### Torsions with `X`
- The `leave_as_X = True` will leave any general parameter terms as `X`.
The parameters will be printed using `999` as the atom type, so if the atoms
that are missing follow those general types, you can add them by hand without
too much digging.
If `leave_as_X = False`, then the generalized dihedrals will be made
into "real" dihedrals for any atom types that are likely relevant.
What is considered relevant is anything that is __NOT__ in the `ions` list of the
`build_X_dihedrals` function.
- Looping through torsions with Xs is absolutely horrific, and you will likely
run into `MAXPRM` issues with pre-compiled or as-is compilations of TINKER.
Because the torsions are resorted after all of the torsions have been built,
important parameters may be overwritten.

### Bugs/Things to Be Aware Of
- Every residue will have a definition for every atom. For example, there are
separate entries (and associated atom types) for the `H1` and `H2` in water,
instead of a single `H` definition.
- This script will NOT write `biotype` information.
- Ions don't exist when built using a `leaprc`.
- If the connectivity information wasn't guessed based off the atom type
(as defined in the `guess_connectivity` function), then it will be listed
as `0` with the comment `!! GUESSED CONNECTION`. No warning will print.
- A warning will print if the rmin14 and eps14 (VDW) parameters weren't found
for a particular atom. They will be set to `-0.0000` for easy finding.

### Things to Test
- Check mol2/frcmod read/write
- Check/update to work with CHARMM/AMOEBA prmtop files (likely an if-statement
change)
- Do I actually need biotypes?!?
