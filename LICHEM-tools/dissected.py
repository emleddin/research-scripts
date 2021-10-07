## dissected.py
## Modify a BurstStruct.xyz and create the BeadStartStruct.xyz
##  maybe doing swapsies along the way
## Requires: BurstStruct.xyz, regions.inp

import MDAnalysis as mda
from shutil import copy2

burst_xyz = "BurstStruct.xyz"
regions_file="regions.inp"
n_beads = 13

## lichem_def_beads: "True" uses the default path MM enviornment
##   and "False" uses the reactant or product MM set-up below.
## Note: By default, LICHEM uses the reactant MM for all beads along the path!
## react_mm_beads and prod_mm_beads: how many beads to use the reactant
##   or product MM for.
##   The two values need to add up to n_bead - 2
##   The react_mm_beads will be used for the first X structures along the path,
##   and prod_mm_beads for the final X structures
##   The reactant and product will keep their current values (hence the `-2`)
## Ex: n_beads = 5, react_mm_beads = 2, prod_mm_beads = 1
## Ex: reactant QM/MM, 1 QM/rxt MM, 2 QM/rxt MM, 3 QM/prod MM, product QM/MM

lichem_def_beads = False
react_mm_beads = 5
prod_mm_beads = 6

##----------------------------- The Curtain :) -----------------------------##

def check_beads(react_mm_beads, prod_mm_beads, n_beads):
    """
    Verifies the input information.
    Parameters
    ----------
    react_mm_beads : int
        The number of beads along the path to substitute the reactant MM for.
    prod_mm_beads : int
        The number of beads along the path to substitute the product MM for.
    n_beads: int
        The number of beads represented by the Burst structure.
    """
    react_mm_beads = int(react_mm_beads)
    prod_mm_beads = int(prod_mm_beads)
    test_beads = react_mm_beads + prod_mm_beads + 2
    if test_beads == n_beads:
        print("Huzzah! Bead math is correct!")
        print(" Keeping the reactant and product the same and using")
        print(f" reactant MM for the first {react_mm_beads} beads and")
        print(f" product MM for the last {prod_mm_beads} beads.\n")
    else:
         raise ValueError("The number of beads specified to the script\n" +
        f"is not correct. My check gave {test_beads}. n_beads: {n_beads}\n"+
        f"Reactant MM beads: {react_mm_beads}, Product MM beads: {prod_mm_beads}\n")

def get_natom(xyz, n_beads):
    """
    Returns the number of atoms from an XYZ file.
    Parameters
    ----------
    xyz : XYZ file
        A standard XYZ file.
    n_beads: int
        The number of beads represented by the Burst structure.
    Returns
    -------
    natom : int
        The number of atoms in the file.
    start_atoms: int
        The total number of atoms to be in the BeadStartStruct.
    """
    with open(xyz) as f:
        natom = f.readlines()[0].strip()
    f.close()
    natom = int(natom)
    #
    start_atoms = natom * n_beads
    #
    return natom, start_atoms

def distribute_burst(burst_xyz, n_beads, natom):
    """
    Breaks the BurstStruct.xyz into multiple XYZ files.
    Parameters
    ----------
    burst_xyz : XYZ file
        A standard XYZ file with n_beads frames.
    n_beads: int
        The number of beads represented by the Burst structure.
    natom : int
        The number of atoms in the file.
    Saves
    -------
    BurstFrame_{}.xyz : XYZ file
        A single frame from the BurstStruct.xyz.
    """
    ## https://stackoverflow.com/questions/16289859/splitting-large-text-file-into-smaller-text-files-by-line-numbers-using-python
    lines_per_file = natom + 2
    struct = None
    bead_counter = 0
    with open(burst_xyz) as burst:
        for lineno, line in enumerate(burst):
            if lineno % lines_per_file == 0:
                if struct:
                    struct.close()
                struct_fname = f'BurstFrame_{bead_counter}.xyz'
                struct = open(struct_fname, "w")
                bead_counter += 1
            struct.write(line)
        if struct:
            struct.close()

def setup_swapsies(react_mm_beads, prod_mm_beads, n_beads):
    """
    Creates a list of which swapsies to use for each frame.
    Parameters
    ----------
    react_mm_beads : int
        The number of beads along the path to substitute the reactant MM for.
    prod_mm_beads : int
        The number of beads along the path to substitute the product MM for.
    n_beads: int
        The number of beads represented by the Burst structure.
    Returns
    -------
    swapsies_list : list
        A list with the type of swapies to perform for each bead in the
        BurstStruct.
    """
    ## Swapsies options: "regular" or "reverse"
    ## regular: put the product QM in the reactant MM
    ## reverse: put the reactant QM in the product MM
    swapsies_list = []
    test_react = 0
    for bead in range(n_beads):
        if bead == 0 or bead == n_beads-1:
            swapsies_list.append("no")
        elif test_react < react_mm_beads:
            swapsies_list.append("regular")
            test_react += 1
        else:
            swapsies_list.append("reverse")
    return swapsies_list

## Use a set: https://stackoverflow.com/questions/740287/how-to-check-if-one-of-the-following-items-is-in-a-list
## Between two strings: https://stackoverflow.com/questions/36559356/extract-values-between-two-strings-in-a-text-file
## Flatten list of lists: https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
def readreg(regions_file):
    '''
    Read the LICHEM regions file and save the QM atom indices. These indices
    correspond to VMD numbering (starting at 0), not TINKER or PDB indices.
    Parameters
    ----------
    regions_file : str
        The file name of the LICHEM regions file.
    Returns
    -------
    all_QM : list
        List of all the indices for atoms in the QM region.
    '''
    f = open(regions_file, 'r')
    reg_lines = f.readlines()
    #
    LICHEM_reg_keywords = ["beads:", "boundary_atoms:", "box_size:",
    "calculation_type:", "electrostatics:", "ensemble:", "eq_steps:",
    "force_constant:", "frozen_atoms:", "frozen_ends:", "init_path_chk:",
    "keep_files:", "lrec_cut:", "lrec_exponent:", "max_stepsize:",
    "max_opt_steps:", "max_qm_steps:", "mm_opt_cut:", "mm_opt_tolerance:",
    "mm_type:", "neb_atoms:", "opt_stepsize:", "pbc:", "potential_type:",
    "pressure:", "prod_steps:", "pseudobond_atoms:", "qm_basis:", "qm_charge:",
    "qm_max_force_tol:", "qm_memory:", "qm_method:", "qm_opt_tolerance:",
    "qm_rms_force_tol:", "qm_spin:", "qm_type:", "qm_units:", "restrain_mm:",
    "solv_model:", "spring_constant:", "ts_freqs:", "use_ewald:", "use_lrec:",
    "use_mm_cutoff:", "use_solvent:"]
    LICHEM_reg_keywords = set(LICHEM_reg_keywords)
    #
    save_QM_atoms = False
    all_QM = []
    #
    for line in reg_lines:
        ## Turns on the line saving for the line following qm_atoms
        ## THANKFULLY this eliminates the needs to remove number of QM atoms
        ## Also, save the number of QM atoms specified in the regions file
        ## the .lower() fixes potential case problems
        if "qm_atoms:" in line.strip().lower():
            save_QM_atoms = True
            junk,total_QM = line.split(":")
            total_QM = total_QM.strip("\n").strip()
        ## Check for keyword that follows
        elif [i for i in LICHEM_reg_keywords if i in line.strip().lower()]:
            save_QM_atoms = False
        ## Save the lines in between with the QM indices
        elif save_QM_atoms:
            all_QM.append(line.strip('\n'))
    f.close()
    #
    ## Clean up the QM list of lists by first splitting at spaces and then
    ## flattening the list of lists
    all_QM = [i.split(' ') for i in all_QM]
    all_QM = [item for sublist in all_QM for item in sublist]
    ## Check that all_QM has correct number of items
    if int(total_QM) == len(all_QM):
        print("It looks like I read the correct number of QM atoms.")
    else:
        print("Aw shucks. Your regions says there should be {} QM atoms, \
I read {}. You should check your regions file, because either \
I read it in wrong, or you have the wrong number of QM atoms listed \
for reading.".format(int(total_QM), len(all_QM)))
    #
    return all_QM

def load_XYZ(lichem_in_xyz):
    """
    Load in the LICHEM XYZ.
    Parameters
    ----------
    lichem_in_xyz: str
        The path to/name of the optimized LICHEM XYZ file.
    Returns
    -------
    system : MDAnalysis.core.universe.Universe
        The LICHEM XYZ.
    """
    system = mda.Universe(lichem_in_xyz, format="XYZ", dt=1.0, in_memory=True)
    # mda.coordinates.XYZ.XYZReader(
    #
    ## Remove the segment IDs (aka `SYSTEM`) for prettier AtomGroup printing
    for atom in system.atoms:
        atom.segment.segid = ''
    #
    ## Set to final XYZ frame
    system.universe.trajectory[-1]
    #
    return system

def write_qm_pdb(system, out_qm_pdb_mda, all_QM):
    """
    Write out a PDB file with the QM atoms.
    Parameters
    ----------
    system: MDAnalysis.core.universe.Universe
        The LICHEM XYZ.
    out_qm_pdb_mda : str
        The name of the PDB file to write out to.
    all_QM : list
        A list of the QM atom indices.
    Saves
    -----
    out_qm_pdb_mda : PDB file
        A PDB file with the QM atoms from the `system` final frame.
    Returns
    -------
    all_QM_ag : MDAnalysis.core.groups.AtomGroup
        An Atom Group consisting of the QM atoms.
    """
    all_QM_ag = mda.AtomGroup(all_QM, system)
    all_QM_ag.write(out_qm_pdb_mda, convert_units=False)
    return all_QM_ag

def read_adjusted_PDB(adj_qm_pdb_mda, all_QM):
    """
    Load in the PDB file with replacement QM atom coordinates.
    Parameters
    ----------
    adj_qm_pdb_mda : str
        A PDB file with the QM atoms from the `system` final frame.
    all_QM : list
        A list of the QM atom indices.
    Returns
    -----
    adj_pdb : MDAnalysis.core.universe.Universe
        The manipulated PDB file.
    all_QM_ag_adj : MDAnalysis.core.universe.Universe
        A copy of adj_pdb.
    """
    adj_pdb = mda.Universe(adj_qm_pdb_mda, format="PDB", dt=1.0, in_memory=True)
    all_QM_ag_adj = adj_pdb.copy()
    return adj_pdb, all_QM_ag_adj

def integrate_movements(system, all_QM, all_QM_ag_adj, lichem_out_xyz):
    """
    Generate the XYZ file with the replacement QM atom coordinates.
    Parameters
    ----------
    system : MDAnalysis.core.universe.Universe
        The LICHEM XYZ.
    all_QM : list
        List of all the indices for atoms in the QM region.
    all_QM_ag_adj : MDAnalysis.core.universe.Universe
        A copy of adj_pdb.
    lichem_out_xyz : str
    Saves
    -----
    lichem_out_xyz : XYZ file
        An XYZ file with the updated QM coordinates.
    """
    all_QM_ag = mda.AtomGroup(all_QM, system)
    all_QM_ag.atoms.positions = all_QM_ag_adj.atoms.positions
    system.atoms.write(lichem_out_xyz, remark='', convert_units=False)

def stitch_beads(lichem_def_beads, natom):
    """
    Stitches a new BurstStruct.xyz together for LICHEM.
    Parameters
    ----------
    lichem_def_beads : bool
        Whether the LICHEM MM region was used (True), or if the swapises
        protocol was followed (False).
    """
    if lichem_def_beads == True:
        find_files = "BurstFrame_"
    else:
        find_files = "FinalBurstFrame_"
    # Create the new BurstStruct.xyz
    with open("new_BurstStruct.xyz", "w+") as bs_out:
        for bnum in range(n_beads):
            with open(find_files + f"{bnum}.xyz") as file:
                lines = file.readlines()
                ## Check if the last line is blank (just "\n")
                last = lines[-1]
                for line in lines[:-1]:
                    bs_out.write(line)
                if last != "\n":
                    bs_out.write(last)
            file.close()
    bs_out.close()
    overall_natoms = natom * n_beads
    # Create the new BeadStartStruct.xyz
    with open("new_BeadStartStruct.xyz", "w+") as bss_out:
        for bnum in range(n_beads):
            if bnum == 0:
                bss_out.write(f"{int(overall_natoms)}\n\n")
            with open(find_files + f"{bnum}.xyz") as file:
                lines = file.readlines()
                ## Skip the natom and blank line
                coords = lines[2:-1]
                ## Check if the last line is blank (just "\n")
                last = lines[-1]
                for line in coords:
                    bss_out.write(line)
                if last != "\n":
                    bss_out.write(last)
            file.close()
    bss_out.close()

##----------------------------------------------------------------------------##
## Run the script!

## Don't run swapsies, just break apart and create new_BeadStartStruct.xyz
##  from the BurstStruct.xyz
if lichem_def_beads == True:
    natom, start_atoms = get_natom(burst_xyz, n_beads)
    distribute_burst(burst_xyz, n_beads, natom)
    stitch_beads(lichem_def_beads, natom)
## Run swapsies if not using the MM default from LICHEM.
elif lichem_def_beads == False:
    check_beads(react_mm_beads, prod_mm_beads, n_beads)
    natom, start_atoms = get_natom(burst_xyz, n_beads)
    distribute_burst(burst_xyz, n_beads, natom)
    swapsies_list = setup_swapsies(react_mm_beads, prod_mm_beads, n_beads)
    #
    all_QM = readreg(regions_file)
    #
    for bnum, swapsies_type in enumerate(swapsies_list):
        ## Use a specific MM environment for the beads
        ## Put the product QM in the reactant MM and reoptimize the product
        print(f"Reading frame {bnum} for {swapsies_type} swapsies.")
        if swapsies_type == "regular":
            print(f"Putting the frame {bnum} QM in the reactant MM.\n")
            ## Set filenames
            lichem_rxt_in_xyz = "BurstFrame_0.xyz"
            lichem_prod_in_xyz = f"BurstFrame_{bnum}.xyz"
            out_qm_prod_pdb_mda = f"product-{bnum}-qm-mda.pdb"
            lichem_out_prod_qm_rxt_mm_xyz = f"FinalBurstFrame_{bnum}.xyz"
            ## Read the product and get the QM atoms
            # all_QM = readreg(regions_file)
            product = load_XYZ(lichem_prod_in_xyz)
            all_QM_ag = write_qm_pdb(product, out_qm_prod_pdb_mda, all_QM)
            ## Put the product QM into the reactant
            reactant = load_XYZ(lichem_rxt_in_xyz)
            adj_pdb, all_QM_ag_adj = read_adjusted_PDB(out_qm_prod_pdb_mda, all_QM)
            integrate_movements(reactant, all_QM, all_QM_ag_adj, lichem_out_prod_qm_rxt_mm_xyz)
        ## Put reactant QM in the product MM and reopt the reactant
        elif swapsies_type == "reverse":
            print(f"Putting the frame {bnum} QM in the product MM.\n")
            ## Set filenames
            lichem_rxt_in_xyz = f"BurstFrame_{bnum}.xyz"
            lichem_prod_in_xyz = f"BurstFrame_{n_beads-1}.xyz"
            out_qm_rxt_pdb_mda = f"reactant-{bnum}-qm-mda.pdb"
            lichem_out_rxt_qm_prod_mm_xyz = f"FinalBurstFrame_{bnum}.xyz"
            ## Read the reactant and get the QM atoms
            # all_QM = readreg(regions_file)
            reactant = load_XYZ(lichem_rxt_in_xyz)
            all_QM_ag = write_qm_pdb(reactant, out_qm_rxt_pdb_mda, all_QM)
            ## Put the reactant QM into the product
            product = load_XYZ(lichem_prod_in_xyz)
            adj_pdb, all_QM_ag_adj = read_adjusted_PDB(out_qm_rxt_pdb_mda, all_QM)
            integrate_movements(product, all_QM, all_QM_ag_adj, lichem_out_rxt_qm_prod_mm_xyz)
        ## Reactant/product case ("no")
        else:
            copy2(f'BurstFrame_{bnum}.xyz', f'FinalBurstFrame_{bnum}.xyz')
    stitch_beads(lichem_def_beads, natom)
## Something broke along the way
else:
    print("You never set `lichem_def_beads` to a True or False value. Exiting...")
