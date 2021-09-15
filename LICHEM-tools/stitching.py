## stitching.py
## Create a complete BeadStartStruct.xyz from two path components.

## XYZ: The name for the reactant-to-intermediate BeadStartStruct.xyz
## Beads: The number of beads given to LICHEM to create BeadStartStruct.xyz
##      lichem -path -b 9 -r reactant.xyz -p intermediate.xyz
##      mv BeadStartStruct.xyz react-int-BSS.xyz
r_i_xyz="react-int-BSS.xyz"
r_i_beads=9

## XYZ: The name for the intermediate-to-product BeadStartStruct.xyz
## Beads: The number of beads given to LICHEM to create BeadStartStruct.xyz
##      lichem -path -b 9 -r intermediate.xyz -p product.xyz
##      mv BeadStartStruct.xyz int-prod-BSS.xyz
i_p_xyz="int-prod-BSS.xyz"
i_p_beads=9

#------------------
def get_natom(xyz):
    """
    Returns the number of atoms from an XYZ file.
    Parameters
    ----------
    xyz : XYZ file
        A standard XYZ file.
    Returns
    -------
    natom : int
        The number of atoms in the file.
    """
    with open(xyz) as f:
        natom = f.readlines()[0].strip()
    f.close()
    natom = int(natom)
    return natom

def num_beads(r_i_beads, i_p_beads):
    """
    Determine the total number of beads based on each step.
    Parameters
    ----------
    r_i_beads : int
        Number of beads represented in the r_i_xyz file.
    i_p_beads : int
        Number of beads represented in the i_p_xyz file.
    Returns
    -------
    tot_beads : int
        The total number of beads, including reactant and product.
    cpu_beads : int
        The number of CPU nodes to use for beads. Does not include reactant and
        product as beads.
    """
    tot_beads = r_i_beads + i_p_beads - 1
    cpu_beads = tot_beads - 2
    #
    return tot_beads, cpu_beads

def stitch_xyzs(r_i_xyz, r_i_natom, r_i_beads, i_p_xyz, i_p_natom, i_p_beads, \
 tot_beads):
    """
    Stitch together the BeadStartStruct.xyz files into one, complete
    BeadStartStruct.xyz.
    Parameters
    ----------
    r_i_xyz : str
        The filename for the BeadStartStruct.xyz for the path from the reactant
        to the intermediate structure.
    r_i_natom : int
        The number of atoms in the r_i_xyz file.
    r_i_beads : int
        Number of beads represented in the r_i_xyz file.
    i_p_xyz : str
        The filename for the BeadStartStruct.xyz for the path from the
        intermediate to the product structure.
    i_p_natom : int
        The number of atoms in the i_p_xyz file.
    i_p_beads : int
        Number of beads represented in the i_p_xyz file.
    tot_beads : int
        The total number of beads, including reactant and product.
    Returns
    -------
    BeadStartStruct.xyz : xyz
        The combined XYZ file.
    """
    r_i_atoms_per_bead = r_i_natom / r_i_beads
    i_p_atoms_per_bead = i_p_natom / i_p_beads
    #
    if r_i_atoms_per_bead != i_p_atoms_per_bead:
        raise ValueError(f"r_i has {r_i_atoms_per_bead} and i_p has {i_p_atoms_per_bead} atoms per bead.\n"+
        "For some reason, your r_i structure and your i_p structure\n" +
        "do not have the same number of atoms in each bead. Try remaking them.\n")
    #
    ## Ignore the intermediate when writing the first time to ensure it's only
    ## written once!
    r_atoms_write = r_i_natom - r_i_atoms_per_bead
    r_atoms_write = int(r_atoms_write)
    #
    ## Sanity check: Do the number being written make sense based on beads?
    overall_natoms = r_atoms_write + i_p_natom
    tot_bead_natoms = r_i_atoms_per_bead * tot_beads
    if overall_natoms != tot_bead_natoms:
        raise ValueError("Either the number of beads specified to the script\n" +
        "or the number of atoms in the structure may be incorrect.\n"+
        f"Total beads: {tot_beads}, Overall natoms: {overall_natoms}\n")
    #
    with open("BeadStartStruct.xyz", "w+") as bss_out:
        ## Write the new natom line
        bss_out.write(f"{int(overall_natoms)}\n\n")
        ## Write the reactant to intermediate XYZ
        with open(r_i_xyz) as file:
            ## Skip the natom and blank line
            for _ in range(2):
                next(file)
            ## Write out all other lines
            count = 0
            for line in file.readlines():
                if count < r_atoms_write:
                    count += 1
                    bss_out.write(line)
                if count == r_atoms_write:
                    count += 2
                    print(f"Check line {count} in BeadStartStruct.xyz if you " +
                    "encounter issues.")
                    print("There may be a blank line or missing atom.")
        file.close()
        ## Write the intermediate to product XYZ
        with open(i_p_xyz) as file:
            ## Skip the natom and blank line
            for _ in range(2):
                next(file)
            ## Write out all other lines
            for line in file.readlines():
                bss_out.write(line)
        file.close()

## Run the script!
r_i_natom = get_natom(r_i_xyz)
i_p_natom = get_natom(i_p_xyz)

tot_beads, cpu_beads = num_beads(r_i_beads, i_p_beads)

stitch_xyzs(r_i_xyz, r_i_natom, r_i_beads, i_p_xyz, i_p_natom, i_p_beads, \
 tot_beads)

print(f"\nVerification: this system is set up for {tot_beads} total beads.\n")
print(f"You should use {cpu_beads} CPU nodes, or a nice divisible number.\n")
print("If you get the 'Restart file does not have the correct format' error,")
print(f"first check that {tot_beads} beads are listed in your regions.inp file.\n")
