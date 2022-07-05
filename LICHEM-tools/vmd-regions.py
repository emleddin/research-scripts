#!/bin/env python3

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
        List of all the indices for the pseudobond atoms.
    all_boundary : list
        List of all the indices for the boundary atoms.
    all_frozen : list
        List of all the indices for the frozen atoms.
    '''
    f = open(regions_file, 'r')
    reg_lines = f.readlines()
    #
    LICHEM_reg_keywords = [
            "acceptance_ratio:", "beads:",
            "boundary_atoms:", "box_size:", "calculation_type:",
            "electrostatics:", "ensemble:", "eq_steps:",
            "force_constant:", "frozen_atoms:", "frozen_ends:",
            "init_path_chk:", "keep_files:",
            "lrec_cut:", "lrec_exponent:",
            "max_stepsize:", "max_opt_steps:", "max_qm_steps:",
            "mm_opt_cut:", "mm_opt_tolerance:", "mm_type:",
            "neb_atoms:", "nqsm:",
            "opt_stepsize:", "pbc:", "potential_type:", "pressure:",
            "prod_steps:", "pseudobond_atoms:",
            "qm_atoms", "qm_basis:", "qm_charge:", "qm_max_force_tol:",
            "qm_memory:", "qm_method:", "qm_opt_tolerance:",
            "qm_rms_force_tol:", "qm_spin:", "qm_type:", "qm_units:",
            "restrain_mm:", "solv_model:", "spring_constant:",
            "ts_freqs:", "use_ewald:", "use_lrec:",
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
    save_NEB_atoms = False
    all_NEB = []
    #
    for line in reg_lines:
        ## Turns on the line saving for the line following qm_atoms
        ## THANKFULLY this eliminates the needs to remove number of QM atoms
        ## Also, save the number of QM atoms specified in the regions file
        ## the .lower() fixes potential case problems
        if "qm_atoms:" in line.strip().lower():
            save_QM_atoms = True
            junk, total_QM = line.split(":")
            total_QM = total_QM.strip("\n").strip()
            if save_PB_atoms is True:
                save_PB_atoms = False
            if save_boundary_atoms is True:
                save_boundary_atoms = False
            if save_frozen_atoms is True:
                save_frozen_atoms = False
            if save_NEB_atoms is True:
                    save_NEB_atoms = False
        elif "pseudobond_atoms:" in line.strip().lower():
            save_PB_atoms = True
            junk, total_PB = line.split(":")
            total_PB = total_PB.strip("\n").strip()
            if save_QM_atoms is True:
                save_QM_atoms = False
            if save_boundary_atoms is True:
                save_boundary_atoms = False
            if save_frozen_atoms is True:
                save_frozen_atoms = False
            if save_NEB_atoms is True:
                    save_NEB_atoms = False
        elif "boundary_atoms:" in line.strip().lower():
            save_boundary_atoms = True
            junk, total_boundary = line.split(":")
            total_boundary = total_boundary.strip("\n").strip()
            if save_QM_atoms is True:
                save_QM_atoms = False
            if save_PB_atoms is True:
                save_PB_atoms = False
            if save_frozen_atoms is True:
                save_frozen_atoms = False
            if save_NEB_atoms is True:
                    save_NEB_atoms = False
        elif "frozen_atoms:" in line.strip().lower():
            save_frozen_atoms = True
            junk, total_frozen = line.split(":")
            total_frozen = total_frozen.strip("\n").strip()
            if save_QM_atoms is True:
                save_QM_atoms = False
            if save_PB_atoms is True:
                save_PB_atoms = False
            if save_boundary_atoms is True:
                save_boundary_atoms = False
            if save_NEB_atoms is True:
                    save_NEB_atoms = False
        elif "neb_atoms:" in line.strip().lower():
            save_NEB_atoms = True
            junk, total_NEB = line.split(":")
            total_NEB = total_NEB.strip("\n").strip()
            if save_QM_atoms is True:
                save_QM_atoms = False
            if save_PB_atoms is True:
                save_PB_atoms = False
            if save_boundary_atoms is True:
                save_boundary_atoms = False
            if save_frozen_atoms is True:
                save_frozen_atoms = False
        ## Check for keyword that follows
        elif [i for i in LICHEM_reg_keywords if i in line.strip().lower()]:
            save_QM_atoms = False
            save_PB_atoms = False
            save_boundary_atoms = False
            save_frozen_atoms = False
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
    # Remove empty strings in QM
    all_QM = list(filter(None, all_QM))
    ## Repeat for PB
    all_PB = [i.split(' ') for i in all_PB]
    all_PB = [item for sublist in all_PB for item in sublist]
    all_PB = list(filter(None, all_PB))
    ## Repeat for boundary
    all_boundary = [i.split(' ') for i in all_boundary]
    all_boundary = [item for sublist in all_boundary for item in sublist]
    all_boundary = list(filter(None, all_boundary))
    ## Repeat for frozen
    all_frozen = [i.split(' ') for i in all_frozen]
    all_frozen = [item for sublist in all_frozen for item in sublist]
    all_frozen = list(filter(None, all_frozen))
    # Repeat for NEB
    all_NEB = [i.split(' ') for i in all_NEB]
    all_NEB = [item for sublist in all_NEB for item in sublist]
    all_NEB = list(filter(None, all_NEB))
    ## Check that all_QM has correct number of items
    if int(total_QM) == len(all_QM):
        print("It looks like I read the correct number of QM atoms.")
    else:
        print("WARNING: Your input regions says there should be\n"
                f"         {int(total_QM)} QM atoms, but\n"
                f"         {len(all_QM)} were read.\n"
                f"         You should verify the specified QM atoms"
                "in generated output.")
    #
    return all_QM, all_PB, all_boundary, all_frozen, all_NEB

def write_vmd_selections(all_QM, all_PB, all_boundary, all_frozen, all_NEB,
                         vmd_filename):
    """
    Creates a VMD preferences file defining macros.

    Returns
    -------
    vmd_filename : vmd
        Text-based file with VMD macro selections.

    Notes
    -----
    VMD macros include qm, QM, quantum, pseudobond, pb, pseudo, boundary, 
    bound, frozen, f, neb, NEB, mm, MM, unfrozen, and uf.
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
        f.write("atomselect macro pseudobond {index ")
        for i in range(len(all_PB)):
            f.write("{} ".format(all_PB[i]))
        f.write("}\n\n")
        f.write("atomselect macro pb {index ")
        for i in range(len(all_PB)):
            f.write("{} ".format(all_PB[i]))
        f.write("}\n\n")
        f.write("atomselect macro pseudo {index ")
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
        # Add NEB
        f.write("atomselect macro neb {index ")
        for i in range(len(all_NEB)):
            f.write("{} ".format(all_NEB[i]))
        f.write("}\n\n")
        f.write("atomselect macro NEB {index ")
        for i in range(len(all_NEB)):
            f.write("{} ".format(all_NEB[i]))
        f.write("}\n\n")
        ## Add mm
        f.write("atomselect macro mm {all not qm or pb}\n")
        f.write("atomselect macro MM {all not qm or pb}\n")
        f.write("\n")
        f.close()

all_QM, all_PB, all_boundary, all_frozen, all_NEB = readreg(regions_file)

write_vmd_selections(all_QM, all_PB, all_boundary, all_frozen, all_NEB,
                     vmd_filename)
