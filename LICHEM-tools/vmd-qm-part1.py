import MDAnalysis as mda

lichem_in_xyz="reactant_optimization_output.xyz"
lichem_out_xyz="xyzfile.xyz"
out_qm_pdb="reactant-qm-vmd.pdb"
regions_file="regions.inp"
vmd_file="view-qm.vmd"

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


def make_vmd(vmd_file, lichem_in_xyz, all_QM, out_qm_pdb):
    """
    Generates a VMD command file that will save a PDB of the QM region.

    Parameters
    ----------
    vmd_file : str
        The name of the VMD file to write to.
    lichem_in_xyz :  str
        The path to/name of the optimized LICHEM XYZ file.
    all_QM : list
        List of all the indices for atoms in the QM region.

    out_qm_pdb : str
        The name of the PDB file to write out to.

    Saves
    -----
    vmd_file : VMD command file
    """
    with open(vmd_file, "w+") as vmd_out:
        vmd_out.write("#!/usr/local/bin/vmd\n")
        vmd_out.write("# VMD Script generated for LICHEM QM region\n")
        vmd_out.write("mol new {} type xyz \n".format(lichem_in_xyz))
        vmd_out.write("set qm [atomselect top \"index ")
        vmd_out.write(" ".join([str(i) for i in all_QM]))
        vmd_out.write(" \"]\n")
        vmd_out.write("# Save the last frame of optimization\n")
        vmd_out.write("$qm frame last\n")
        vmd_out.write("$qm writepdb {}\n".format(out_qm_pdb))
        vmd_out.write("# Change to CPK\n")
        vmd_out.write("mol selection {index ")
        vmd_out.write(" ".join([str(i) for i in all_QM]))
        vmd_out.write(" }\n")
        vmd_out.write("mol representation CPK\n")
        vmd_out.write("mol addrep top\n")
        vmd_out.write("# Remove lines view\n")
        vmd_out.write("mol delrep 0 top\n")
        vmd_out.write("# Showcase QM\n")
        vmd_out.write("display resetview\n")
        vmd_out.close()

## Run the script!

all_QM = readreg(regions_file)
make_vmd(vmd_file, lichem_in_xyz, all_QM, out_qm_pdb)
