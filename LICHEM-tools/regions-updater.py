#!/bin/env python3

class LICHEMParameters:
    """
    Parameters from a LICHEM "regions.inp" file.

    Methods
    -------
    use_dfp(DFP_criteria="tight")
        Set parameters for DFP optimization.
    use_qsm(restrain_QSM=True, beads=3)
        Set parameters for restrained or unrestrained QSM.
    write(outfile="updated_regions.inp")
        Write out an updated parameter file.
    """
    def __init__(self, regions_file="regions.inp"):
        '''
        Parse the LICHEM regions file.

        Parameters
        ----------
        regions_file : str
            The file name of the LICHEM regions file.
        '''
        self.regions_file = regions_file
        # Set up option to append value
        def append_val(line):
            """Save keyword option specified on one-line."""
            junk, var = line.split(":")
            var = var.strip()
            return var
        #
        f = open(self.regions_file, 'r')
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
            "neb_atoms:", "nqsm",
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
            # Ignore lines starting with comments
            if line.startswith("#") or line.startswith("!") or \
              line.startswith("//"):
                continue
            # Turns on the line saving for the line following qm_atoms
            # THANKFULLY this eliminates the needs to remove number of
            #  QM atoms
            # Also, save the number of QM atoms specified in the regions file
            #  the .lower() fixes potential case problems
            elif "qm_atoms:" in line.strip().lower():
                save_QM_atoms = True
                junk, total_QM = line.split(":")
                self.total_QM = total_QM.strip("\n").strip()
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
                self.total_PB = total_PB.strip("\n").strip()
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
                self.total_boundary = total_boundary.strip("\n").strip()
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
                self.total_frozen = total_frozen.strip("\n").strip()
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
                self.total_NEB = total_NEB.strip("\n").strip()
                if save_QM_atoms is True:
                    save_QM_atoms = False
                if save_PB_atoms is True:
                    save_PB_atoms = False
                if save_boundary_atoms is True:
                    save_boundary_atoms = False
                if save_frozen_atoms is True:
                    save_frozen_atoms = False
            # Check for keyword that follows
            elif [i for i in LICHEM_reg_keywords if i in line.strip().lower()]:
                save_QM_atoms = False
                save_PB_atoms = False
                save_boundary_atoms = False
                save_frozen_atoms = False
                # Set other keywords
                lsl = line.strip().lower()
                if "acceptance_ratio:" in lsl:
                    self.acceptance_ratio = append_val(line)
                elif "pbc:" in lsl:
                    self.pbc = append_val(line)
                elif "beads:" in lsl:
                    self.beads = append_val(line)
                elif "box_size:" in lsl:
                    self.box_size = append_val(line)
                    if len(self.box_size.split()) > 3:
                        print("\nWARNING: Read box size has more than 3 "
                              "coordinates (XYZ).\n"
                              f"         Check that {self.box_size} "
                              "was intentional.")
                    elif len(self.box_size.split()) < 3:
                        print("\nWARNING: Read box size has less than 3 "
                              "coordinates (XYZ).\n"
                              f"         Check that {self.box_size} "
                              "was intentional.")
                elif "calculation_type:" in lsl:
                    self.calculation_type = append_val(line)
                elif "electrostatics:" in lsl:
                    self.electrostatics = append_val(line)
                    # TODO: Verify AMBER is acceptable
                    if self.electrostatics.lower() in ("charges", "charge", 
                      "point-charge"):
                        print("\nNOTE: Setting LREC_exponent to 2 "
                                "because {self.electrostatics} was given "
                                "for electrostatics.")
                        self.lrec_exponent = 2
                    elif self.electrostatics.lower() == "amoeba":
                        print("\nNOTE: Setting LREC_exponent to 3 "
                                "because {self.electrostatics} was given "
                                "for electrostatics.")
                        self.lrec_exponent = 3
                    # "gem" case?
                elif "ensemble:" in lsl:
                    self.box_size = append_val(line)
                elif "eq_steps:" in lsl:
                    self.eq_steps = append_val(line)
                elif "force_constant:" in lsl:
                    self.force_constant = append_val(line)
                elif "frozen_atoms:" in lsl:
                    self.frozen_atoms = append_val(line)
                elif "frozen_ends:" in lsl:
                    self.frozen_ends = append_val(line)
                elif "init_path_chk:" in lsl:
                    self.init_path_chk = append_val(line)
                elif "keep_files:" in lsl:
                    self.keep_files = append_val(line)
                elif "lrec_cut:" in lsl:
                    self.lrec_cut = append_val(line)
                elif "max_stepsize:" in lsl:
                    self.max_stepsize = append_val(line)
                elif "max_opt_steps:" in lsl:
                    self.max_opt_steps = append_val(line)
                elif "max_qm_steps:" in lsl:
                    self.max_qm_steps = append_val(line)
                elif "mm_opt_cut:" in lsl:
                    self.mm_opt_cut = append_val(line)
                elif "mm_opt_tolerance:" in lsl:
                    self.mm_opt_tolerance = append_val(line)
                elif "mm_type:" in lsl:
                    self.mm_type = append_val(line)
                elif "nqsm:" in lsl:
                    self.nqsm = append_val(line)
                elif "opt_stepsize:" in lsl:
                    self.opt_stepsize = append_val(line)
                elif "pbc:" in lsl:
                    self.pbc = append_val(line)
                elif "potential_type:" in lsl:
                    self.potential_type = append_val(line)
                elif "pressure:" in lsl:
                    self.pressure = append_val(line)
                elif "prod_steps:" in lsl:
                    self.prod_steps = append_val(line)
                elif "qm_basis:" in lsl:
                    self.qm_basis = append_val(line)
                elif "qm_charge:" in lsl:
                    self.qm_charge = append_val(line)
                elif "qm_max_force_tol:" in lsl:
                    self.qm_max_force_tol = append_val(line)
                elif "qm_memory:" in lsl:
                    self.qm_memory = append_val(line)
                elif "qm_method:" in lsl:
                    self.qm_method = append_val(line)
                elif "qm_opt_tolerance:" in lsl:
                    self.qm_opt_tolerance = append_val(line)
                elif "qm_rms_force_tol:" in lsl:
                    self.qm_rms_force_tol = append_val(line)
                elif "qm_spin:" in lsl:
                    self.qm_spin = append_val(line)
                elif "qm_type:" in lsl:
                    self.qm_type = append_val(line)
                elif "qm_units:" in lsl:
                    self.qm_units = append_val(line)
                elif "restrain_mm:" in lsl:
                    self.restrain_mm = append_val(line)
                elif "solv_model:" in lsl:
                    self.solv_model = append_val(line)
                elif "spring_constant:" in lsl:
                    self.spring_constant = append_val(line)
                elif "ts_freqs:" in lsl:
                    self.ts_freqs = append_val(line)
                elif "use_ewald:" in lsl:
                    self.use_ewald = append_val(line)
                elif "use_lrec:" in lsl:
                    self.use_lrec = append_val(line)
                elif "use_mm_cutoff:" in lsl:
                    self.use_mm_cutoff = append_val(line)
                elif "use_solvent:" in lsl:
                    self.use_solvent = append_val(line)
                # elif ":" in lsl:
                #     self. = append_val(line)
                # else:
            # Save the lines in between with the QM indices
            elif save_QM_atoms:
                all_QM.append(line.strip('\n'))
            # Save the lines in between with the QM indices
            elif save_PB_atoms:
                all_PB.append(line.strip('\n'))
            elif save_boundary_atoms:
                all_boundary.append(line.strip('\n'))
            elif save_frozen_atoms:
                all_frozen.append(line.strip('\n'))
            elif save_NEB_atoms:
                all_NEB.append(line.strip('\n'))
        #
        f.close()
        #
        # Clean up the QM list of lists by first splitting at spaces and then
        # flattening the list of lists
        all_QM = [i.split(' ') for i in all_QM]
        all_QM = [item for sublist in all_QM for item in sublist]
        # Remove empty strings in QM
        all_QM = list(filter(None, all_QM))
        # Repeat for PB
        all_PB = [i.split(' ') for i in all_PB]
        all_PB = [item for sublist in all_PB for item in sublist]
        # Remove empty strings in PB
        all_PB = list(filter(None, all_PB))
        # Repeat for boundary
        all_boundary = [i.split(' ') for i in all_boundary]
        all_boundary = [item for sublist in all_boundary for item in sublist]
        # Remove empty strings in BA
        all_boundary = list(filter(None, all_boundary))
        # Repeat for frozen
        all_frozen = [i.split(' ') for i in all_frozen]
        all_frozen = [item for sublist in all_frozen for item in sublist]
        # Remove empty strings in FA
        all_frozen = list(filter(None, all_frozen))
        # Repeat for NEB
        all_NEB = [i.split(' ') for i in all_NEB]
        all_NEB = [item for sublist in all_NEB for item in sublist]
        # Remove empty strings in NEB
        all_NEB = list(filter(None, all_NEB))
        # Check that all_QM has correct number of items
        if int(total_QM) == len(all_QM):
            print("It looks like I read the correct number of QM atoms.")
        else:
            print("WARNING: Your input regions says there should be\n"
                  f"         {int(total_QM)} QM atoms, but\n"
                  f"         {len(all_QM)} were read.\n"
                  f"         You should verify the specified QM atoms"
                  "in generated output.")
        #
        # Update overall params
        self.qm_atoms = all_QM
        self.pseudobond_atoms = all_PB
        self.boundary_atoms = all_boundary
        self.frozen_atoms = all_frozen
        self.neb_atoms = all_NEB
    #
    def use_dfp(self, DFP_criteria="tight"):
        """
        Set parameters for DFP optimizations.

        Parameters
        ----------
        DFP_criteria : str
            Criteria to use: "loose", "medium", or "tight".
        """
        self.DFP_criteria = DFP_criteria
        def check_attr(attribute):
            """Check whether the attribute is assigned and delete."""
            if hasattr(self, attribute):
                delattr(self, attribute)
            return
        # Set for all
        self.calculation_type = "DFP"
        self.opt_stepsize = 0.50
        self.max_stepsize = 0.10
        self.max_opt_steps = 30
        self.max_qm_steps = 15
        # Base on given type
        if DFP_criteria.lower() == "loose":
            self.qm_opt_tolerance = 0.15
            self.qm_rms_force_tol = 0.10
            self.qm_max_force_tol = 0.02
            self.mm_opt_tolerance = 0.20
        elif DFP_criteria.lower() == "medium":
            self.qm_opt_tolerance = 0.05
            self.qm_rms_force_tol = 0.10
            self.qm_max_force_tol = 0.015
            self.mm_opt_tolerance = 0.05
        elif DFP_criteria.lower() in ("tight", "default"):
            self.qm_opt_tolerance = 0.001
            self.qm_rms_force_tol = 0.005
            self.qm_max_force_tol = 0.015
            self.mm_opt_tolerance = 0.01
        else:
            print("\nWARNING: DFP criteria not understood.\n "
                  "        Accepted values include 'loose', 'medium', "
                  "and 'tight'.\n "
                  "        Using 'tight' criteria as default.")
            self.qm_opt_tolerance = 0.001
            self.qm_rms_force_tol = 0.005
            self.qm_max_force_tol = 0.015
            self.mm_opt_tolerance = 0.01
        # Unset the QSM parameters
        check_attr("beads")
        check_attr("nqsm")
        check_attr("ts_freq")
        check_attr("frozen_ends")
        check_attr("restrain_mm")
        check_attr("force_constant")
        return self
    #
    def use_qsm(self, restrain_QSM=True, beads=3):
        """
        Set parameters for QSM optimizations.

        Parameters
        ----------
        restrain_QSM : bool
            Criteria for restrained QSM (True) or unrestrained QSM (False).
        beads : int
            Number of images in QSM path (including reactant and product).
        """
        self.restrain_QSM = restrain_QSM
        self.beads = beads
        self.calculation_type = "QSM"
        if self.beads == 3:
            print("\nWARNING: Using 3 beads for QSM.\n"
                  "         This value should likely be changed. Users "
                  "commonly request 7+ beads\n"
                  "         (1 reactant, 1 product, 5+ images between).")
        def check_attr(attribute):
            """Check whether the attribute is assigned and delete."""
            if hasattr(self, attribute):
                delattr(self, attribute)
            return
        if restrain_QSM:
            self.opt_stepsize = 0.01
            self.max_stepsize = 0.05
            self.nqsm = 0
            self.ts_freq = "yes"
            self.frozen_ends = "yes"
            self.restrain_mm = "yes"
            self.force_constant = 50
            self.qm_opt_tolerance = 0.001
            self.qm_rms_force_tol = 0.01
            self.qm_max_force_tol = 0.02
            self.mm_opt_tolerance = 0.1
            self.max_opt_steps = 40
            self.max_qm_steps = 10
        else:
            self.opt_stepsize = 0.01
            self.max_stepsize = 0.05
            self.nqsm = 0
            self.ts_freq = "yes"
            self.frozen_ends = "yes"
            self.restrain_mm = "no"
            # Delete attribute if it exists
            check_attr("force_constant")
            self.qm_opt_tolerance = 0.001
            self.qm_rms_force_tol = 0.005
            self.qm_max_force_tol = 0.01
            self.mm_opt_tolerance = 0.1
            self.max_opt_steps = 3
            self.max_qm_steps = 10
        # print("\nNOTE: Prepared for QSM. "
        #       "You may want to add a section,\n"
        #       "      NEB_atoms, for those atoms involved in the "
        #       "reaction.\n")
        return self
    #
    # TODO: Add SP, FREQ, OPT, STEEP, NEB, PIMC, and FBNEB configurations
    # TODO: Standardize yes/true vs no/false responses
    # TODO: Verify that given values are typed correctly (float, int, etc.)
    #
    def write(self, outfile="updated_regions.inp"):
        """
        Write out an updated regions.inp file.

        Parameters
        ----------
        outfile : str
            Name of the output file.
        """
        self.outfile = outfile
        print(f"\nWriting output to {outfile}.")
        def check_attr(ofile_name, attribute):
            """Check whether the attribute is assigned."""
            if hasattr(self, attribute):
                f.write(f"{ofile_name}: {getattr(self, attribute)}\n")
            return
        # Write the update regions file
        with open(self.outfile, 'w') as f:
            # If attributes exist, write in specific order
            # Name_in_output_file, Name_of_class_attribute
            check_attr("Potential_type", "potential_type")
            check_attr("QM_type", "qm_type")
            check_attr("QM_method", "qm_method")
            check_attr("QM_basis", "qm_basis")
            check_attr("QM_memory", "qm_memory")
            check_attr("QM_charge", "qm_charge")
            check_attr("QM_spin", "qm_spin")
            check_attr("MM_type", "mm_type")
            check_attr("Electrostatics", "electrostatics")
            check_attr("Calculation_type", "calculation_type")
            check_attr("Opt_stepsize", "opt_stepsize")
            check_attr("Max_stepsize", "max_stepsize")
            check_attr("beads", "beads")
            check_attr("nqsm", "nqsm")
            check_attr("ts_freq", "ts_freq")
            check_attr("frozen_ends", "frozen_ends")
            check_attr("restrain_mm", "restrain_mm")
            check_attr("force_constant", "force_constant")
            check_attr("qm_opt_tolerance", "qm_opt_tolerance")
            check_attr("qm_rms_force_tol", "qm_rms_force_tol")
            check_attr("qm_max_force_tol", "qm_max_force_tol")
            check_attr("mm_opt_tolerance", "mm_opt_tolerance")
            check_attr("max_opt_steps", "max_opt_steps")
            check_attr("max_qm_steps", "max_qm_steps")
            check_attr("PBC", "pbc")
            check_attr("Box_size", "box_size")
            check_attr("Use_LREC", "use_lrec")
            check_attr("LREC_cut", "lrec_cut")
            check_attr("LREC_exponent", "lrec_exponent")
            check_attr("Use_Ewald", "use_ewald")
            check_attr("Keep_files", "keep_files")
            # TODO: Verify that these make sense in this order
            check_attr("acceptance_ratio", "acceptance_ratio")
            check_attr("ensemble", "ensemble")
            check_attr("eq_steps", "eq_steps")
            check_attr("init_path_chk", "init_path_chk")
            check_attr("pressure", "pressure")
            check_attr("prod_steps", "prod_steps")
            check_attr("qm_units", "qm_units")
            check_attr("solv_model", "solv_model")
            check_attr("spring_constant", "spring_constant")
            check_attr("use_mm_cutoff", "use_mm_cutoff")
            check_attr("use_solvent", "use_solvent")
            #
            # check_attr("", "")
            #
            # Finish with the lists
            # Write NEB (only if they exist and doing QSM or NEB)
            if self.calculation_type in ("QSM", "NEB"):
                if hasattr(self, "total_NEB"):
                    f.write(f"NEB_atoms: {self.total_NEB}\n")
                    for i, item in enumerate(self.neb_atoms):
                        # If i is divisible by 10 or last item
                        if (i+1) % 10 == 0 or item == self.neb_atoms[-1]:
                            f.write(f"{item}\n")
                        else:
                            f.write(f"{item} ")
                else:
                    print("NOTE: The calculation_type is "
                          f"{self.calculation_type}, but the "
                          "'NEB_atoms' section is blank.\n"
                          "      That section is required for NEB "
                          "calculations and optional for QSM.")
            # Write QM
            f.write(f"QM_atoms: {self.total_QM}\n")
            for i, item in enumerate(self.qm_atoms):
                # If i is divisible by 10 or last item
                if (i+1) % 10 == 0 or item == self.qm_atoms[-1]:
                    f.write(f"{item}\n")
                else:
                    f.write(f"{item} ")
            # Write PB atoms
            f.write(f"Pseudobond_atoms: {self.total_PB}\n")
            for i, item in enumerate(self.pseudobond_atoms):
                # If i is divisible by 10 or last item
                if (i+1) % 10 == 0 or item == self.pseudobond_atoms[-1]:
                    f.write(f"{item}\n")
                else:
                    f.write(f"{item} ")
            # Write BA
            f.write(f"Boundary_atoms: {self.total_boundary}\n")
            for i, item in enumerate(self.boundary_atoms):
                # If i is divisible by 10 or last item
                if (i+1) % 10 == 0 or item == self.boundary_atoms[-1]:
                    f.write(f"{item}\n")
                else:
                    f.write(f"{item} ")
            # Write frozen
            f.write(f"Frozen_atoms: {self.total_frozen}\n")
            for i, item in enumerate(self.frozen_atoms):
                # If i is divisible by 10 or last item
                if (i+1) % 10 == 0 or item == self.frozen_atoms[-1]:
                    f.write(f"{item}\n")
                else:
                    f.write(f"{item} ")
            f.close()
            return self


# END CLASS
# -------------------------------------------------------------------------- #

# Read regions file into class
params = LICHEMParameters()
# params = LICHEMParameters("qsm-regions.inp")

# Print all attributes
# print(vars(params))

# Update any values here!!!!
# params.qm_spin = 12345

# Write basic
# params.write()

# QSM Examples
# ------------
# params.use_qsm().write("qsm_reg_rest.inp")

# Don't need to label restrain_QSM/beads, I did it to keep myself straight
# You can also do the steps individually if you wish
# params.use_qsm(restrain_QSM=False, beads=15)
# params.write("qsm_reg_unrest.inp")

# DFP Examples
# ------------
# Don't need to label DFP_criteria, I did it to keep myself straight
# params.use_dfp(DFP_criteria="loose").write("dfp_reg_loose.inp")
#
# params.use_dfp(DFP_criteria="default").write("dfp_reg_default.inp")
#
# params.use_dfp(DFP_criteria="abcd").write("dfp_reg_abcd.inp")
