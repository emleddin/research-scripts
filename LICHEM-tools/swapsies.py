import MDAnalysis as mda

## Swapsies: Put the product QM in the reactant MM and reoptimize the product
lichem_rxt_in_xyz="reactant_optimization_output.xyz"
lichem_prod_in_xyz="product_optimization_output.xyz"
regions_file="regions.inp"
out_qm_prod_pdb_mda="product-qm-mda.pdb"
lichem_out_prod_qm_rxt_mm_xyz="xyzfile.xyz"

## Reverse Swapsies: Put reactant QM in the product MM and reopt the reactant
## Switch the commenting at the end of the script for this option!!!
#lichem_rxt_in_xyz="reactant_optimization_output.xyz"
#lichem_prod_in_xyz="product_optimization_output.xyz"
#regions_file="regions.inp"
#out_qm_rxt_pdb_mda="reactant-qm-mda.pdb"
#lichem_out_rxt_qm_prod_mm_xyz="xyzfile.xyz"

## Use a set: https://stackoverflow.com/questions/740287/how-to-check-if-one-of-the-following-items-is-in-a-list
## Between two strings: https://stackoverflow.com/questions/36559356/extract-values-between-two-strings-in-a-text-file
## Flatten list of lists: https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
def readreg(regions_file):
    '''Read the LICHEM regions file and save the QM atom indices. These indices
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
    all_QM_ag.write(out_qm_pdb_mda)
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
    system.atoms.write(lichem_out_xyz, remark='')

## Run the script!

## -- For Swapsies
## Read the product and get the QM atoms
all_QM = readreg(regions_file)
product = load_XYZ(lichem_prod_in_xyz)
all_QM_ag = write_qm_pdb(product, out_qm_prod_pdb_mda, all_QM)

## Put the product QM into the reactant
reactant = load_XYZ(lichem_rxt_in_xyz)
adj_pdb, all_QM_ag_adj = read_adjusted_PDB(out_qm_prod_pdb_mda, all_QM)
integrate_movements(reactant, all_QM, all_QM_ag_adj, lichem_out_prod_qm_rxt_mm_xyz)

#-----------------------------------------------------------------------------#

#--- For Reverse Swapsies
## Read the reactant and get the QM atoms
#all_QM = readreg(regions_file)
#reactant = load_XYZ(lichem_rxt_in_xyz)
#all_QM_ag = write_qm_pdb(reactant, out_qm_rxt_pdb_mda, all_QM)
#
### Put the reactant QM into the product
#product = load_XYZ(lichem_prod_in_xyz)
#adj_pdb, all_QM_ag_adj = read_adjusted_PDB(out_qm_rxt_pdb_mda, all_QM)
#integrate_movements(product, all_QM, all_QM_ag_adj, lichem_out_rxt_qm_prod_mm_xyz)
