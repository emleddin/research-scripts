regions_file="regions.inp"
vmd_filename="test-vmd-selections.vmd"

def readreg(regions_file):
    '''Read the LICHEM regions file and save the atom indices. These indices
    correspond to VMD numbering (starting at 0), not TINKER or PDB indices.
    Parameters
    ----------
    regions_file : str
        The file name of the LICHEM regions file.
    Returns
    -------
    all_QM : list
        List of all the indices for atoms in the QM region.
    all_PB : list
        List of all the indices for the psuedobond atoms.
    all_boundary : list
        List of all the indices for the boundary atoms.
    all_frozen : list
        List of all the indices for the frozen atoms.
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
    save_PB_atoms = False
    all_PB = []
    #
    save_boundary_atoms = False
    all_boundary = []
    #
    save_frozen_atoms = False
    all_frozen = []
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
            if save_PB_atoms == True:
                save_PB_atoms = False
            if save_boundary_atoms == True:
                save_boundary_atoms = False
            if save_frozen_atoms == True:
                save_false_atoms = False
        elif "pseudobond_atoms:" in line.strip().lower():
            save_PB_atoms = True
            junk,total_PB = line.split(":")
            total_PB = total_PB.strip("\n").strip()
            if save_QM_atoms == True:
                save_QM_atoms = False
            if save_boundary_atoms == True:
                save_boundary_atoms = False
            if save_frozen_atoms == True:
                save_false_atoms = False
        elif "boundary_atoms:" in line.strip().lower():
            save_boundary_atoms = True
            junk,total_boundary = line.split(":")
            total_boundary = total_boundary.strip("\n").strip()
            if save_QM_atoms == True:
                save_QM_atoms = False
            if save_PB_atoms == True:
                save_PB_atoms = False
            if save_frozen_atoms == True:
                save_false_atoms = False
        elif "frozen_atoms:" in line.strip().lower():
            save_frozen_atoms = True
            junk,total_frozen = line.split(":")
            total_frozen = total_frozen.strip("\n").strip()
            if save_QM_atoms == True:
                save_QM_atoms = False
            if save_PB_atoms == True:
                save_PB_atoms = False
            if save_boundary_atoms == True:
                save_boundary_atoms = False
        ## Check for keyword that follows
        elif [i for i in LICHEM_reg_keywords if i in line.strip().lower()]:
            save_QM_atoms = False
            save_PB_atoms = False
            save_boundary_atoms = False
            save_false_atoms = False
        ## Save the lines in between with the QM indices
        elif save_QM_atoms:
            all_QM.append(line.strip('\n'))
        ## Save the lines in between with the QM indices
        elif save_PB_atoms:
            all_PB.append(line.strip('\n'))
        elif save_boundary_atoms:
            all_boundary.append(line.strip('\n'))
        elif save_frozen_atoms:
            all_frozen.append(line.strip('\n'))

    f.close()
    #
    ## Clean up the QM list of lists by first splitting at spaces and then
    ## flattening the list of lists
    all_QM = [i.split(' ') for i in all_QM]
    all_QM = [item for sublist in all_QM for item in sublist]
    ## Repeat for PB
    all_PB = [i.split(' ') for i in all_PB]
    all_PB = [item for sublist in all_PB for item in sublist]
    ## Repeat for boundary
    all_boundary = [i.split(' ') for i in all_boundary]
    all_boundary = [item for sublist in all_boundary for item in sublist]
    ## Repeat for frozen
    all_frozen = [i.split(' ') for i in all_frozen]
    all_frozen = [item for sublist in all_frozen for item in sublist]
    ## Check that all_QM has correct number of items
    if int(total_QM) == len(all_QM):
        print("It looks like I read the correct number of QM atoms.")
    else:
        print("Aw shucks. Your regions says there should be {} QM atoms, \
I read {}. You should check your regions file, because either \
I read it in wrong, or you have the wrong number of QM atoms listed \
for reading.".format(int(total_QM), len(all_QM)))
    #
    return all_QM, all_PB, all_boundary, all_frozen

def write_vmd_selections(all_QM, all_PB, all_boundary, all_frozen, vmd_filename):
    """Creates a VMD preferences file defining macros using the atom lists from
    the regions file. These include qm, QM, quantum, pseudobond, pb, pseudo,
    boundary, bound, frozen, f, unfrozen, and uf.
    Returns
    -------
    vmd_filename : vmd
        Text-based file with VMD macro selections.
    """
    with open(vmd_filename, "w+") as f:
        f.write("\n")
        f.write("atomselect macro qm {index ")
        for i in range(len(all_QM)):
            f.write("{} ".format(all_QM[i]))
        f.write("}\n\n")
        f.write("atomselect macro QM {index ")
        for i in range(len(all_QM)):
            f.write("{} ".format(all_QM[i]))
        f.write("}\n\n")
        f.write("atomselect macro quantum {index ")
        for i in range(len(all_QM)):
            f.write("{} ".format(all_QM[i]))
        f.write("}\n\n")
        f.write("atomselect macro psuedobond {index ")
        for i in range(len(all_PB)):
            f.write("{} ".format(all_PB[i]))
        f.write("}\n\n")
        f.write("atomselect macro pb {index ")
        for i in range(len(all_PB)):
            f.write("{} ".format(all_PB[i]))
        f.write("}\n\n")
        f.write("atomselect macro psuedo {index ")
        for i in range(len(all_PB)):
            f.write("{} ".format(all_PB[i]))
        f.write("}\n\n")
        f.write("atomselect macro boundary {index ")
        for i in range(len(all_boundary)):
            f.write("{} ".format(all_boundary[i]))
        f.write("}\n\n")
        f.write("atomselect macro bound {index ")
        for i in range(len(all_boundary)):
            f.write("{} ".format(all_boundary[i]))
        f.write("}\n\n")
        f.write("atomselect macro ba {index ")
        for i in range(len(all_boundary)):
            f.write("{} ".format(all_boundary[i]))
        f.write("}\n\n")
        f.write("atomselect macro frozen {index ")
        for i in range(len(all_frozen)):
            f.write("{} ".format(all_frozen[i]))
        f.write("}\n\n")
        f.write("atomselect macro f {index ")
        for i in range(len(all_frozen)):
            f.write("{} ".format(all_frozen[i]))
        f.write("}\n\n")
        ## Add unfrozen
        f.write("atomselect macro unfrozen {all not frozen}\n\n")
        f.write("atomselect macro uf {all not frozen}\n")
        f.write("\n")
        f.close()

all_QM, all_PB, all_boundary, all_frozen = readreg(regions_file)

write_vmd_selections(all_QM, all_PB, all_boundary, all_frozen, vmd_filename)
