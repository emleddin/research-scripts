import MDAnalysis as mda
import numpy as np
import parmed as pmd

orig_pdb="WT_protein_system_frame_23456.pdb"
tink_xyz="WT_protein_system_frame_23456_convert.xyz"

## Atom number for center of active atom shell
shell_center=1234
## Did you use the index from VMD for the shell_center? If yes, set True
VMD_index_shell=False

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
        print("Oopsie! This PDB didn't contain box information. Setting coordinates to 0\n\
 for right now. Please find this information and update your regions file.\n")
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
    system = mda.Universe(orig_pdb, tink_xyz, format="TXYZ", dt=1.0, in_memory=True)

    ## Remove the segment IDs (aka `SYSTEM`) for prettier AtomGroup printing
    for atom in system.atoms:
        atom.segment.segid = ''

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

def select_QM(system, shell_center):
    """
    Select the QM atoms using the atom selection language of MDAnalysis.

    Paramters
    ---------
    system : MDAnalysis.core.universe.Universe
        The Tinker XYZ information mapped onto a PDB topology.

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

        >>> QM_WAT = system.select_atoms("(around 4 resnum 200) or \
         (around 4 resnum 250) and (resname WAT)")
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

def make_regions(x_box, y_box, z_box, all_QM, all_PB, all_BA, all_FR):
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
    nc = 10 ## Number of columns to print
    QM_len = (len(all_QM)-(len(all_QM)%(nc)))+(nc)
    QM_range = QM_len//(nc)
    print_QM = [all_QM.atoms.ix[i*(nc):i*(nc)+(nc)] for i in range(QM_range)]
    #
    PB_len = (len(all_PB)-(len(all_PB)%(nc)))+(nc)
    PB_range = PB_len//(nc)
    print_PB = [all_PB.atoms.ix[i*(nc):i*(nc)+(nc)] for i in range(PB_range)]
    #
    BA_len = (len(all_BA)-(len(all_BA)%(nc)))+(nc)
    BA_range = BA_len//(nc)
    print_BA = [all_BA.atoms.ix[i*(nc):i*(nc)+(nc)] for i in range(BA_range)]
    #
    FR_len = (len(all_FR)-(len(all_FR)%(nc)))+(nc)
    FR_range = FR_len//(nc)
    print_FR = [all_FR[i*(nc):i*(nc)+(nc)] for i in range(FR_range)]
    with open("regions.inp_backup", "w+") as reg_out:
        reg_out.write("Potential_type: QMMM\n")
        reg_out.write("QM_type: g16\n")
        # to match current
        reg_out.write("!QM_type:Gaussian\n")
        reg_out.write("QM_method: B3LYP\n")
        reg_out.write("QM_basis: GEN\n")
        reg_out.write("QM_memory: 80 GB\n")
        reg_out.write("QM_charge: -4\n")
        reg_out.write("QM_spin: 1\n")
        reg_out.write("MM_type: TINKER\n")
        reg_out.write("Electrostatics: CHARGES\n")
        reg_out.write("Calculation_type: SP\n")
        reg_out.write("Opt_stepsize: 0.50\n")
        reg_out.write("Max_stepsize: 0.10\n")
        reg_out.write("qm_opt_tolerance: 1e-3\n")
        reg_out.write("qm_rms_force_tol: 0.010\n")
        reg_out.write("qm_max_force_tol: 0.020\n")
        reg_out.write("mm_opt_tolerance: 1e-1\n")
        reg_out.write("max_opt_steps: 10\n")
        reg_out.write("max_qm_steps: 10\n")
        reg_out.write("PBC: Yes\n")
        reg_out.write("Box_size: {:.6f} {:.6f} {:.6f}\n".format(x_box, y_box, z_box))
        reg_out.write("Use_LREC: Yes\n")
        reg_out.write("LREC_cut: 25.0\n")
        reg_out.write("Use_Ewald: Yes\n")
        reg_out.write("Keep_files: Yes\n")
        reg_out.write("QM_atoms: {}\n".format(len(all_QM)))
        reg_out.write('\n'.join(' '.join(map(str,sl)) for sl in print_QM) + '\n')
        reg_out.write("Pseudobond_atoms: {}\n".format(len(all_PB)))
        reg_out.write('\n'.join(' '.join(map(str,sl)) for sl in print_PB) + '\n')
        reg_out.write("Boundary_atoms: {}\n".format(len(all_BA)))
        reg_out.write('\n'.join(' '.join(map(str,sl)) for sl in print_BA) + '\n')
        reg_out.write("Frozen_atoms: {}\n".format(len(all_FR)))
        reg_out.write('\n'.join(' '.join(map(str,sl)) for sl in print_FR))
        reg_out.close()

def map_BASIS(all_QM, all_PB):
    """
    Generates a reference between the atoms in the `regions.inp` file and their
    numbering in the Gaussian BASIS file.

    Parameters
    ----------
    all_QM : MDAnalysis.core.groups.AtomGroup
        An atom group of all of the selected atoms for the QM region.
    all_PB : MDAnalysis.core.groups.AtomGroup
        An atom group of all of the selected pseudobond atoms.
    """
    ## Make a group of all the atoms in the BASIS file
    ## Do a union (| instead of +) because you want them to be reordered
    BASIS = all_QM.atoms | all_PB.atoms
    count = 1
    with open("BASIS_list.txt", "w+") as bl_out:
        bl_out.write("BASIS_ID Regions_ID TINKER_ID ResName ResNum AtomName\n")
        for atom in BASIS.atoms:
            bl_out.write("{:<8} {:<10} {:<9} {:<7} {:<6} {:<8}\n".format(count, \
            atom.ix, atom.id, atom.resname, atom.resnum, atom.name))
            count += 1
        bl_out.close()

#-------------- Run the program -------------#

## Get box information
x_box, y_box, z_box = get_box(orig_pdb)

## Load the system
system = load_XYZ(orig_pdb, tink_xyz)

## Correct the shell_center
shell_center = check_shell(shell_center, VMD_index_shell)

## Select QM, Basis, and Pseudobond atoms
all_QM, all_BA, all_PB, all_FR = select_QM(system, shell_center)

## Make the regions file
make_regions(x_box, y_box, z_box, all_QM, all_PB, all_BA, all_FR)

## Create a key between the regions file IDs and the BASIS numbering
map_BASIS(all_QM, all_PB)

print("I have successfully generated the regions file. Good luck on the next \n\
steps of the QM/MM process!!!")
