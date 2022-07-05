import MDAnalysis as mda
import numpy as np
import parmed as pmd
import pandas as pd

orig_pdb="WT_protein_system_frame_23456.pdb"
tink_xyz="WT_protein_system_frame_23456_convert.xyz"

## Alternative: using system arguments
## The script itself is sys.argv[0]
# orig_pdb = sys.argv[1]
# tink_xyz = sys.argv[2]

## Atom number for center of active atom shell
shell_center=1234
## Did you use the index from VMD for the shell_center? If yes, set True
VMD_index_shell=False

## Specify AMBER (aka "CHARGES" for LICHEM) or AMOEBA
electro="AMBER"
## Starting criteria -- loose, medium, or tight
criteria="tight"
## Set the QM parameters for the resulting regions file!
method="B3LYP"
mem="80 GB"
charge="-2"
spin="1"

## Set the regular basis level and the higher basis level for atoms defined
## in the select_higher_basis function
reg_basis_level="6-31G*"
high_basis_level="6-31+G(d,p)"
## The information for the pseudobond atoms (assuming all are the same type)
## The linebreaks are important! :)
#PB1="\n[PB basis set]"
#PB2="[PB pseudopotentials]"
PB1="STO-2G\nSP 2 1.00\n0.9034 1.00 1.00\n0.21310 1.90904 0.57864"
PB2="try1 1 2\nS Component\n1\n1 7.75 16.49\nP\n1\n1 1.0 0.0"

##------------- Functions requiring modification!!!!
def select_QM(system, shell_center):
    """
    Select the QM atoms using the atom selection language of MDAnalysis.

    Parameters
    ---------
    system : MDAnalysis.core.universe.Universe
        The Tinker XYZ information mapped onto a PDB topology.

    Returns
    -------
    all_QM : MDAnalysis.core.groups.AtomGroup
        An atom group of all of the selected atoms for the QM region.
    all_PB : MDAnalysis.core.groups.AtomGroup
        An atom group of all of the selected pseudobond atoms.
    all_BA : MDAnalysis.core.groups.AtomGroup
        An atom group of all of the selected boundary atoms.
    all_FR : list
        A list of all of the frozen atoms.

    Examples
    ---------
    Here are some select examples. Learn more in the MDAnalysis documentation.

    To use the values from VMD indices, use slices of `system.atoms[]`.
    This option is not inclusive. So, the following pulls atoms that have a VMD
    index of [194, 195, ..., 214, 215], which can be checked by listing the
    atoms.

        >>> QM_LYS_15 = system.atoms[194:216]
        >>> list(QM_LYS_15[:])
        [<Atom 195: N of type N of resname LYS, resid 15 and segid and altLoc >,
        <Atom 196: H of type H of resname LYS, resid 15 and segid and altLoc >,
        ...,
        <Atom 215: C of type C of resname LYS, resid 15 and segid and altLoc >,
        <Atom 216: O of type O of resname LYS, resid 15 and segid and altLoc >]

    To use the atom numbers from the TINKER XYZ, you can use the `bynum`
    keyword with `system.select_atoms()`. This option is inclusive.

        >>> QM_LYS_15 = system.select_atoms("bynum 195:216")
        >>> list(QM_LYS_15[:])
        [<Atom 195: N of type N of resname LYS, resid 15 and segid and altLoc >,
        <Atom 196: H of type H of resname LYS, resid 15 and segid and altLoc >,
        ...,
        <Atom 215: C of type C of resname LYS, resid 15 and segid and altLoc >,
        <Atom 216: O of type O of resname LYS, resid 15 and segid and altLoc >]

    Entire residues can be selected with either `resnum` or `resid`. These
    match the residue numbers found in the PDB.

        >>> QM_LYS_15 = system.select_atoms("resnum 15")
        >>> list(QM_LYS_15[:])
        [<Atom 195: N of type N of resname LYS, resid 15 and segid and altLoc >,
        <Atom 196: H of type H of resname LYS, resid 15 and segid and altLoc >,
        ...,
        <Atom 215: C of type C of resname LYS, resid 15 and segid and altLoc >,
        <Atom 216: O of type O of resname LYS, resid 15 and segid and altLoc >]

        >>> QM_LYS_15 = system.select_atoms("resid 15")
        >>> list(QM_LYS_15[:])
        [<Atom 195: N of type N of resname LYS, resid 15 and segid and altLoc >,
        <Atom 196: H of type H of resname LYS, resid 15 and segid and altLoc >,
        ...,
        <Atom 215: C of type C of resname LYS, resid 15 and segid and altLoc >,
        <Atom 216: O of type O of resname LYS, resid 15 and segid and altLoc >]

    Atom masks can select atoms that match (or don't match) from a residue
    number.

        >>> QM_LYS_15 = system.select_atoms("resnum 15 and (name C or name O)")
        >>> list(QM_LYS_15[:])
        [<Atom 215: C of type C of resname LYS, resid 15 and segid and altLoc >,
        <Atom 216: O of type O of resname LYS, resid 15 and segid and altLoc >]

    Distance masks can be used to select atoms within in a certain distance
    cutoff (in angstroms). This example would select any atoms from water
    residues within 4 Å of two specified residues.

        >>> QM_WAT = system.select_atoms("(around 4 resnum 475) or \
         (around 4 resnum 477) and (resname WAT)")
        >>> QM_WAT = QM_WAT.residues.atoms
        >>> list(QM_WAT[:])
        [<Atom 22150: O of type O of resname WAT, resid 5154 and segid and altLoc >,
        ...,
        <Atom 49182: H2 of type H of resname WAT, resid 14164 and segid and altLoc >]
    """
    ## My QM Atoms
    QM_AAA_20 = system.select_atoms("resnum 20 and (name C or name O)")
    QM_BBB_22 = system.select_atoms("bynum 50:55")
    QM_CCC_24 = system.select_atoms("resnum 24")
    QM_DDD_26 = system.select_atoms("resnum 26")
    QM_EEE_28 = system.select_atoms("resnum 28 and (name N or name H)")
    QM_FFF_30 = system.select_atoms("bynum 288:293")
    QM_GGG_32 = system.select_atoms("bynum 303:308")
    QM_DX_100 = system.select_atoms("resnum 100 and name O3\'")
    QM_DX3_101 = system.select_atoms("resnum 101")
    QM_ME_200 = system.select_atoms("resnum 200")
    QM_NTP_225 = system.select_atoms("resnum 225")
    QM_ME_250 = system.select_atoms("resnum 250")
    QM_WAT = system.select_atoms("(around 4 resnum 200) or (around 4 resnum 250) and (resname WAT)")
    QM_WAT = QM_WAT.residues.atoms
    #
    ## Combine the QM atoms. Consider using `|` instead of '+' to make `all_QM`
    ## ordered with a single copy of an atom.
    all_QM = QM_AAA_20 | QM_BBB_22 | QM_CCC_24 | QM_DDD_26 | QM_EEE_28 | QM_FFF_30 | \
    QM_GGG_32 | QM_DX_100 | QM_DX3_101 | QM_ME_200 | QM_NTP_225 | QM_ME_250 | \
    QM_WAT
    #
    ## My Pseudo Atoms
    PB_AAA_20 = system.select_atoms("resnum 20 and name CA")
    PB_BBB_22 = system.select_atoms("resnum 22 and name CA")
    PB_EEE_28 = system.select_atoms("resnum 28 and name CA")
    PB_FFF_30 = system.select_atoms("resnum 30 and name CA")
    PB_GGG_32 = system.select_atoms("resnum 32 and name CB")
    PB_DX_100 = system.select_atoms("resnum 100 and name C3\'")
    #
    ## Combine the PB atoms. Consider using `|` instead of '+' to make `all_PB`
    ## ordered with a single copy of an atom.
    all_PB = PB_AAA_20 | PB_BBB_22 | PB_EEE_28 | PB_FFF_30 | PB_GGG_32 | PB_DX_100
    #
    ## My Boundary Atoms
    BA_AAA_20 = system.select_atoms("resnum 20 and (name N or name H or name HA \
     or name CB or name HB or name CG2 or name HG21 or name HG22 or name HG23 \
     or name CG1 or name HG12 or name HG13 or name CD1 or name HD11 or name HD12 \
     or name HD13)")
    BA_BBB_22 = system.select_atoms("resnum 22 and (name N or name H or name HA \
     or name C or name O)")
    BA_EEE_28 = system.select_atoms("resnum 28 and (name HA or name C or name O or \
     name CB or name HB1 or name HB2 or name HB3)")
    BA_FFF_30 = system.select_atoms("resnum 30 and (name N or name H or name HA \
     or name C or name O)")
    BA_GGG_32 = system.select_atoms("resnum 32 and (name CA or name HB2 or name \
     HB3)")
    BA_DX_100 = system.select_atoms("resnum 100 and (name C4\' or name O4\' or \
     name C1\' or name C2\' or name H4\' or name H3\' or name H2\' or name H2\'\' \
     or name H1\')")
    #
    ## Combine the BA atoms. Consider using `|` instead of '+' to make `all_BA`
    ## ordered with a single copy of an atom.
    all_BA = BA_AAA_20 | BA_BBB_22 | BA_EEE_28 | BA_FFF_30 | BA_GGG_32 | BA_DX_100
    #
    print("There are {} QM, {} pseudobond, and {} boundary atoms.\n".format( \
     len(all_QM), len(all_PB), len(all_BA)))
    #
    ## Redo the sphere for the unfrozen list
    my_sphere = system.select_atoms("sphzone 15.0 (index %d)" %shell_center)
    #
    print_SC = system.atoms[shell_center]
    #
    print("You used residue {} {} at atom {} {} for the sphere center.\n"\
     .format(print_SC.residue.resname, print_SC.residue.resnum, print_SC.id, \
     print_SC.name))
    print("There are {} active atoms.\n".format(len(my_sphere)))
    #
    ## Get a list of all the unfrozen atoms
    unfrozen = np.concatenate((my_sphere.ix, all_QM.atoms.ix, all_PB.atoms.ix,
     all_BA.atoms.ix))
    #
    print("There are {} unfrozen atoms in my array.\n".format(len(unfrozen)))
    #
    tot = len(my_sphere.ix) + len(all_QM.atoms.ix) + len(all_PB.atoms.ix) + \
     len(all_BA.atoms.ix)
    #
    print("There should be {} total unfrozen atoms.\n".format(tot))
    #
    ## Get array of all atoms
    all_ix = system.atoms.ix
    #
    ## len(all_ix) - len(unfrozen) doesn't need to match the total number of
    ## frozen atoms because the sphere will overlap with QM, BA, and PB
    #
    all_FR = [i for i in all_ix if i not in unfrozen]
    #
    print("There are {} total frozen atoms.\n".format(len(all_FR)))
    return all_QM, all_BA, all_PB, all_FR


def select_higher_basis(system):
    """
    Select the QM atoms using the atom selection language of MDAnalysis.

    Parameters
    ---------
    system : MDAnalysis.core.universe.Universe
        The Tinker XYZ information mapped onto a PDB topology.

    Returns
    -------
    all_HB : MDAnalysis.core.groups.AtomGroup
        An atom group of all of the atoms in the QM region selected for a
        higher basis.
    """
    HB_DX3_101 = system.select_atoms("resnum 101 and (name O3' or name HO3')")
    HB_ME_200 = system.select_atoms("resnum 200")
    HB_NTP_225 = system.select_atoms("resnum 225 and (name O3A or name PA or name O5')")
    HB_ME_250 = system.select_atoms("resnum 250")
    #
    ## Combine the higher basis atoms. Consider using `|` instead of '+'
    ## to make `all_HB` ordered with a single copy of an atom.
    all_HB = HB_DX3_101 | HB_ME_200 | HB_NTP_225 | HB_ME_250
    return all_HB


##------------- Standard Function definitions (no modification required)
## You may want to modify `make_regions` if you know you want add'l keywords
## in your `regions.inp_backup` file

def get_box(orig_pdb):
    """
    Determines the box size from the PDB for the regions file.

    Parameters
    ----------
    orig_pdb : str
        The path to the PDB file used to generate the QM/MM TINKER XYZ.

    Returns
    -------
    x_box, y_box, z_box : float
        X, Y, and Z box size coordinates from the input PDB.
    """
    pdb = pmd.load_file(orig_pdb)
    save_box = pdb.get_box()
    if save_box.any() == None:
        print("Oopsie! This PDB didn't contain box information. "
              "Setting coordinates to 0\nfor right now. "
              "Please find this information and update your regions file.\n")
        x_box = 0.
        y_box = 0.
        z_box = 0.
    else:
        x_box = save_box.item(0)
        y_box = save_box.item(1)
        z_box = save_box.item(2)
    #
    print("\nBox size: {} {} {}\n".format(x_box, y_box, z_box))
    return x_box, y_box, z_box


def load_XYZ(orig_pdb, tink_xyz):
    """
    Load in the Tinker XYZ using the PDB as a topology.

    Parameters
    ----------
    orig_pdb : str
        The path to the PDB file used to generate the QM/MM TINKER XYZ.
    tink_xyz: str
        The path to the TINKER XYZ file for QM/MM.

    Returns
    -------
    system : MDAnalysis.core.universe.Universe
        The Tinker XYZ information mapped onto a PDB topology.
    """
    system = mda.Universe(orig_pdb, tink_xyz, format="TXYZ", dt=1.0, 
                          in_memory=True)
    #
    ## Remove the segment IDs (aka `SYSTEM`) for prettier AtomGroup printing
    for atom in system.atoms:
        atom.segment.segid = ''
    #
    return system

def check_shell(shell_center, VMD_index_shell):
    """
    Selects the correct atom for shell center, based on VMD or TINKER indexing.

    Parameters
    ----------
    shell_center : int
        The atom index.
    VMD_index_shell : bool
        `True` if shell_center given as VMD index, `False` if given as TINKER
        index.

    Returns
    -------
    shell_center : int
        The atom index using the VMD index.
    """
    if VMD_index_shell == False:
        shell_center -= 1
    return shell_center

def make_regions(x_box, y_box, z_box, all_QM, all_PB, all_BA, all_FR, electro,\
     criteria, method, mem, charge, spin):
    """
    Generates the regions.inp file. The quantum, pseudobond, boundary, and
    frozen atom lists are formmated for printing in 10 columns.

    Parameters
    ----------
    x_box, y_box, z_box : float
        X, Y, and Z box size coordinates from the input PDB.
    all_QM : MDAnalysis.core.groups.AtomGroup
        An atom group of all of the selected atoms for the QM region.
    all_PB : MDAnalysis.core.groups.AtomGroup
        An atom group of all of the selected pseudobond atoms.
    all_BA : MDAnalysis.core.groups.AtomGroup
        An atom group of all of the selected boundary atoms.
    all_FR : list
        A list of all of the frozen atoms.
    """
    ##
    electro = electro.upper()
    criteria = criteria.lower()
    #
    nc = 10 ## Number of columns to print
    #
    with open("regions.inp_backup", "w+") as reg_out:
        reg_out.write("Potential_type: QMMM\n")
        reg_out.write("QM_type: g16\n")
        # to match current
        reg_out.write("!QM_type:Gaussian\n")
        reg_out.write("QM_method: {}\n".format(method))
        reg_out.write("QM_basis: GEN\n")
        reg_out.write("QM_memory: {}\n".format(mem))
        reg_out.write("QM_charge: {}\n".format(charge))
        reg_out.write("QM_spin: {}\n".format(spin))
        reg_out.write("MM_type: TINKER\n")
        if electro == "AMBER" or "CHARGES":
            reg_out.write("Electrostatics: CHARGES\n")
        elif electro == "AMOEBA":
            reg_out.write("Electrostatics: AMOEBA\n")
        reg_out.write("Calculation_type: SP\n")
        reg_out.write("Opt_stepsize: 0.50\n")
        reg_out.write("Max_stepsize: 0.10\n")
        if criteria == "loose":
            reg_out.write("qm_opt_tolerance: 0.15\n")
            reg_out.write("qm_rms_force_tol: 0.10\n")
            reg_out.write("qm_max_force_tol: 0.020\n")
            reg_out.write("mm_opt_tolerance: 0.20\n")
        elif criteria == "medium":
            reg_out.write("qm_opt_tolerance: 0.05\n")
            reg_out.write("qm_rms_force_tol: 0.010\n")
            reg_out.write("qm_max_force_tol: 0.015\n")
            reg_out.write("mm_opt_tolerance: 0.05\n")
        elif criteria == "tight":
            reg_out.write("qm_opt_tolerance: 0.001\n")
            reg_out.write("qm_rms_force_tol: 0.005\n")
            reg_out.write("qm_max_force_tol: 0.015\n")
            reg_out.write("mm_opt_tolerance: 0.01\n")
        reg_out.write("max_opt_steps: 30\n")
        reg_out.write("max_qm_steps: 15\n")
        reg_out.write("PBC: Yes\n")
        reg_out.write(f"Box_size: {x_box:.6f} {y_box:.6f} {z_box:.6f}\n")
        reg_out.write("Use_LREC: Yes\n")
        reg_out.write("LREC_cut: 25.0\n")
        if electro == "AMBER":
            reg_out.write("LREC_exponent: 2\n")
        elif electro == "AMOEBA":
            reg_out.write("LREC_exponent: 3\n")
        reg_out.write("Use_Ewald: Yes\n")
        reg_out.write("Keep_files: Yes\n")
        reg_out.write(f"QM_atoms: {len(all_QM)}\n")
        # These are atom groups
        for i, item in enumerate(all_QM, start=1):
            # If i is divisible by number of columns or last item
            if (i) % nc == 0 or item == all_QM[-1]:
                reg_out.write(f"{item.ix}\n")
            else:
                reg_out.write(f"{item.ix} ")
        # Write PB atoms
        reg_out.write(f"Pseudobond_atoms: {len(all_PB)}\n")
        for i, item in enumerate(all_PB, start=1):
            # If i is divisible by number of columns or last item
            if (i) % nc == 0 or item == all_PB[-1]:
                reg_out.write(f"{item.ix}\n")
            else:
                reg_out.write(f"{item.ix} ")
        # Write BA
        reg_out.write(f"Boundary_atoms: {len(all_BA)}\n")
        for i, item in enumerate(all_BA, start=1):
            # If i is divisible by number of columns or last item
            if (i) % nc == 0 or item == all_BA[-1]:
                reg_out.write(f"{item.ix}\n")
            else:
                reg_out.write(f"{item.ix} ")
        # Write frozen
        # This is a list
        reg_out.write(f"Frozen_atoms: {len(all_FR)}\n")
        for i, item in enumerate(all_FR, start=1):
            # If last item
            if item == all_FR[-1]:
                reg_out.write(f"{item}")
            # Else if i is divisible by number of columns
            elif (i) % nc == 0:
                reg_out.write(f"{item}\n")
            else:
                reg_out.write(f"{item} ")
        reg_out.close()

def map_BASIS(all_QM, all_PB):
    """
    Generates a reference DataFrame map between the atoms in the `regions.inp`
    file and their numbering in the Gaussian BASIS file.

    Parameters
    ----------
    all_QM : MDAnalysis.core.groups.AtomGroup
        An atom group of all of the selected atoms for the QM region.
    all_PB : MDAnalysis.core.groups.AtomGroup
        An atom group of all of the selected pseudobond atoms.

    Returns
    -------
    basis_df : pandas.core.frame.DataFrame
        A DataFrame with the BASIS file mapping.
    """
    ## Make a group of all the atoms in the BASIS file
    ## Do a union (| instead of +) because you want them to be reordered
    BASIS = all_QM.atoms | all_PB.atoms
    basis_df = pd.DataFrame([atom.ix, atom.id, atom.resname, atom.resnum,
     atom.name] for atom in BASIS.atoms)
    ## Add ColNames
    basis_df.columns = ["Regions_ID", "TINKER_ID", "ResName", "ResNum", "AtomName"]
    ## Add index as column (BASIS listing)
    basis_df.insert(0, "BASIS_ID", basis_df.index+1)
    # ## Write a file with the map items
    # with open("BASIS_list.txt", "w+") as bl_out:
    #     bl_out.write("BASIS_ID Regions_ID TINKER_ID ResName ResNum AtomName\n")
    #     for r in basis_df.itertuples(index=True, name='Pandas'):
    #         bl_out.write("{:<8} {:<10} {:<9} {:<7} {:<6} {:<8}\n".format(\
    #         r.BASIS_ID, r.Regions_ID, r.TINKER_ID, r.ResName, r.ResNum, \
    #         r.AtomName))
    #     bl_out.close()
    return basis_df

def make_BASIS(basis_df, all_QM, all_HB, all_PB, reg_basis_level,
 high_basis_level, PB1, PB2):
    """
    Generates the Gaussian BASIS file.

    Parameters
    ----------
    basis_df : pandas.core.frame.DataFrame
        A DataFrame with the BASIS file mapping.

    Returns
    -------
    BASIS : txt
        The Gaussian file with basis set information.
    """
    ## Get the "regular BASIS" atoms
    reg_B = all_QM.subtract(all_HB)
    ## Create an empty column
    basis_df["Type"] = np.nan
    #
    ## Assign "REG" type to atoms with regular BASIS
    for index, row in basis_df.iterrows():
        for atom in reg_B.atoms:
            if atom.ix == row["Regions_ID"]:
                basis_df.loc[basis_df.index[index], 'Type'] = "REG"
    #
    ## Do the "higher BASIS" atoms
    for index, row in basis_df.iterrows():
        for atom in all_HB.atoms:
            if atom.ix == row["Regions_ID"]:
                basis_df.loc[basis_df.index[index], 'Type'] = "HIGH"
    #
    ## Do the Pseudobond atoms
    for index, row in basis_df.iterrows():
        for atom in all_PB.atoms:
            if atom.ix == row["Regions_ID"]:
                basis_df.loc[basis_df.index[index], 'Type'] = "PB"
    #
    ## Separate into different DataFrames
    reg = basis_df.loc[basis_df['Type'] == "REG"].reset_index(drop=True)
    high = basis_df.loc[basis_df['Type'] == "HIGH"].reset_index(drop=True)
    pb = basis_df.loc[basis_df['Type'] == "PB"].reset_index(drop=True)
    #
    ## Turn the BASIS_IDs for each DF into a list
    reg_list = reg["BASIS_ID"].to_list()
    high_list = high["BASIS_ID"].to_list()
    pb_list = pb["BASIS_ID"].to_list()
    #
    ## Write out the sanity check
    with open("BASIS_verification.txt", "w+") as bv_out:
        bv_out.write("BASIS_ID Regions_ID TINKER_ID ResName ResNum AtomName Type\n")
        for r in basis_df.itertuples(index=True, name='Pandas'):
            bv_out.write("{:<8} {:<10} {:<9} {:<7} {:<6} {:<8} {:<5}\n".format(\
            r.BASIS_ID, r.Regions_ID, r.TINKER_ID, r.ResName, r.ResNum, \
            r.AtomName, r.Type))
        bv_out.close()
    #
    ## Now write the BASIS file
    with open("BASIS", "w+") as b_out:
    #---- Regular Basis
        ## Write indices of those with regular BASIS level in groups of 8
        for i, item in enumerate(reg_list, start=1):
        # If i is divisible by 8 or last item
            if i % 8 == 0 or item == reg_list[-1]:
                b_out.write(f"{item}  0\n")
                b_out.write("{}\n".format(reg_basis_level))
                b_out.write("****\n")
            ## If there is a remainder:
            else:
                b_out.write(f"{item}")
    #---- Higher Basis
        # ## Write indices of those with higher BASIS level in groups of 8
        for i, item in enumerate(high_list, start=1):
        # If i is divisible by 8 or last item
            if i % 8 == 0 or item == high_list[-1]:
                b_out.write(f"{item}  0\n")
                b_out.write("{}\n".format(high_basis_level))
                b_out.write("****\n")
            ## If there is a remainder:
            else:
                b_out.write(f"{item}")
    #---- PB Atoms
        ## Part 1: write indices of PB atoms in groups of 12
        for i, item in enumerate(high_list, start=1):
        # If i is divisible by 12 or last item
            if i % 12 == 0 or item == high_list[-1]:
                b_out.write(f"{item}  0 {PB1.rstrip()}\n")
                b_out.write("****\n")
            ## If there is a remainder:
            else:
                b_out.write(f"{item}")
        ## Part 2: rewrite indices with additional info
        for i, item in enumerate(high_list, start=1):
        # If i is divisible by 12 or last item
            if i % 12 == 0 or item == high_list[-1]:
                b_out.write(f"{item}  0\n")
                b_out.write(f"{PB2.rstrip()}\n\n")
            ## If there is a remainder:
            else:
                b_out.write(f"{item}")
        b_out.close()
        if len(pb_list) > 12:
            print("You may need to check the PB section formatting in the "
                  "BASIS file.\n"
                  "I did the best I could, but more than 12 pseudobond atoms "
                  "gets funky.\n")
    return basis_df

#-------------- Run the program -------------#

## Get box information
x_box, y_box, z_box = get_box(orig_pdb)

## Load the system
system = load_XYZ(orig_pdb, tink_xyz)

## Correct the shell_center
shell_center = check_shell(shell_center, VMD_index_shell)

## Select QM, Boundary, and Pseudobond atoms
all_QM, all_BA, all_PB, all_FR = select_QM(system, shell_center)

## Make the regions file
make_regions(x_box, y_box, z_box, all_QM, all_PB, all_BA, all_FR, electro, \
 criteria, method, mem, charge, spin)

print("I have successfully generated the regions file.\n")

## Create a DataFrame relating the regions file IDs and the BASIS numbering
basis_df = map_BASIS(all_QM, all_PB)

## Select atoms using the higher basis
all_HB = select_higher_basis(system)

## Create the BASIS file
basis_df = make_BASIS(basis_df, all_QM, all_HB, all_PB, reg_basis_level, \
 high_basis_level, PB1, PB2)

print("I have successfully generated the BASIS file. Good luck on the next \n\
steps of the QM/MM process!!!")
