#!/usr/bin/env python 3

# import MDAnalysis as mda
import numpy as np
import pandas as pd
import pathlib
import re
# import time

# Required LICHEM files
basisfile = "BASIS"
confile = "connect.inp"
regfile = "regions.inp"
xyzfile = "xyzfile.xyz"

out_dir = "./output"

# -------------------------- Behind the Curtain -------------------------- #

# Tried to mess with the precision for read_csv, to no avail.
# pd.set_option('display.float_format', '{:,16f}'.format)
# pd.set_option('display.precision', 16)
# pd.set_option('styler.format.precision', 16)
# np.set_printoptions(precision=16)

class GaussianGENFile:
    """
    Information from a Gaussian GEN BASIS file.

    Methods
    -------
    read(read_file="BASIS")
        Read in the BASIS file.
    write(write_dir="./output", write_file="BASIS")
        Write out the BASIS file.
    """
    def __init__(self, basisfile="BASIS", out_dir="./output", outfile="BASIS"):
        self.basisfile = basisfile
        self.out_dir = out_dir
        self.outfile = outfile
        return
    def read(self, read_file=None):
        '''
        Parse the Gaussian basis file used with the `QM_basis: GEN` argument.
        '''
        if read_file is None:
            read_file = self.basisfile
        self.basis_atoms = {}
        #
        def b_rstrip(orig_string, replacement='\s'):
            """
            Strip only the end of a string matching a regex.

            Parameters
            ----------
            orig_string : str
                The complete string.
            replacement : str
                What to remove from the string.

            Returns
            -------
            stripped_content : str
                Cleaned string.
            """
            replace_regex = re.compile(f'^({replacement})*|({replacement})*$')
            stripped_content = replace_regex.sub('', orig_string)
            return stripped_content
        #
        pseudo_pot = False
        skip_next_line = False
        self.pseudo_pot_info = []
        #
        with open(read_file, 'r+') as f:
            lines = f.readlines()
            for line_pointer, line in enumerate(lines):
                # Store next line (if it exists)
                if line == lines[-1]:
                    last_line = True
                else:
                    next_line = lines[line_pointer+1]
                # Add basis to dictionary
                # Catch lines to skip
                if skip_next_line is True:
                    skip_next_line = False
                # Get lines with functional info
                elif line.endswith("  0\n") and pseudo_pot is False:
                    # Can't use regular rstrip in case 0 is
                    #  elsewhere in the string
                    aline = b_rstrip(line, "  0\n")
                    # Normal strip is okay because it only removes newline
                    functional = next_line.strip()
                    # Save atoms as dictionary key, since many will have
                    #  same functional as value
                    for i in aline.split():
                        self.basis_atoms[str(i)] = functional
                # If pseudopotential line
                # TODO: Figure out how these are defined for other people
                #       so it's not just hardcoded as STO-2G
                elif line.endswith("0 STO-2G\n") and pseudo_pot is False:
                    pseudo_pot = True
                    # Can't use regular rstrip in case 0/2 are
                    #  elsewhere in the string
                    aline = b_rstrip(line, "0 STO-2G\n")
                    self.pseudo_pot_info.append("_PB_ATOM_PH_STO-2G")
                    # Save atoms to dictionary with "PB" value
                    for i in aline.split():
                        self.basis_atoms[str(i)] = "PB"
                    # break
                # TODO: Similarly, this part will break if multiple blocks
                elif pseudo_pot is True:
                    if line.startswith("\n"):
                        # If very last 2 lines in the file
                        if last_line is True:
                            self.pseudo_pot_info.append(line)
                        else:
                            self.pseudo_pot_info.append(line)
                            self.pseudo_pot_info.append("_PB_ATOM_PH")
                            # Skip the line of atoms
                            skip_next_line = True
                    else:
                        # Save list of lines
                        self.pseudo_pot_info.append(line.strip())
                # line.startswith("****"):
                else:
                    continue
        # Verify that it looks right
        # print(self.pseudo_pot_info)
        # print(self.basis_atoms)
        return self
    #
    # TODO: Remove the atom indices you think you should
    # def remove():
    #     return self
    #
    # TODO: Write out a new BASIS file
    def write(self, write_dir=None, write_file=None):
        """
        Write the GEN basis file.

        Parameters
        ----------
        write_dir : str
            Output directory.
        write_file : str
            Name of the file to write to.
        """
        if write_dir is None:
            write_dir = self.out_dir
        if write_file is None:
            write_file = self.outfile
        # Specify output path
        self.out_path = write_dir + '/' + write_file
        # Create the output directory (no error if it doesn't exist)
        pathlib.Path(write_dir).mkdir(parents=True, exist_ok=True)
        # Write to file
        # If n_connections is not 0, put a space before connections
        with open(self.out_path, 'w+') as f:
            # PB atoms
            matching_PB = [k for k,v in self.basis_atoms.items() \
                                    if v == 'PB']
            # Use the values to make a "dictionary of keys" for functionals,
            #  then force to a list.
            # This preserves the order from the original dict,
            #  since sets are unordered.
            all_funcs = list(dict.fromkeys(self.basis_atoms.values()))
            # Move "PB" to the very end of the list
            all_funcs.append(all_funcs.pop(all_funcs.index("PB")))
            for value in all_funcs:
                # Get list of atoms with same functional
                if value == "PB":
                    for line in self.pseudo_pot_info:
                        if line == "_PB_ATOM_PH_STO-2G":
                            for item in matching_PB:
                                if item == matching_PB[-1]:
                                    f.write(f"{item}  0 STO-2G\n")
                                else:
                                    f.write(f"{item} ")
                        elif line == "_PB_ATOM_PH":
                            for item in matching_PB:
                                if item == matching_PB[-1]:
                                    f.write(f"{item}  0\n")
                                else:
                                    f.write(f"{item} ")
                        elif line == "\n":
                            f.write(line)
                        else:
                            f.write(f"{line}\n")
                else:
                    matching_func = [k for k,v in self.basis_atoms.items() \
                                    if v == value]
                    # Write to file in increments of 8
                    for i, item in enumerate(matching_func, start=1):
                        # If i is divisible by 8 or last item
                        if (i) % 8 == 0 or item == matching_func[-1]:
                            f.write(f"{item}  0\n")
                            f.write(f"{value}\n****\n")
                        else:
                            f.write(f"{item} ")
        return


class LICHEMConnectFile:
    """
    Information from a LICHEM connect.inp file.

    Methods
    -------
    read(read_file="connect.inp")
        Read in the connectivity file.
    write(write_dir="./output", write_file="connect.inp")
        Write out the connectivity file.
    """
    def __init__(self, confile="connect.inp", out_dir="./output",
                 outfile="connect.inp"):
        self.confile = confile
        self.out_dir = out_dir
        self.outfile = outfile
        return
    def read(self, read_file=None):
        """
        Parse the LICHEM connectivity file.
        """
        if read_file is None:
            read_file = self.confile
        # Extend to 6 possible connections
        col_names = ['LICHEM_ID', 'atom_name', 'tinker_type', 'mass', 'charge',
                    'n_connections', 'con_a', 'con_b', 'con_c', 'con_d',
                    'con_e', 'con_f']
        # Read in the connectivity file, keep NAs blank
        self.con_df = pd.read_csv(
                    read_file,
                    low_memory=False,
                    names = col_names,
                    delim_whitespace=True, na_filter=False)
        # Combine all con columns into a one-column space-separated list
        self.con_df['CT'] = self.con_df[['con_a', 'con_b', 'con_c', 'con_d',
                            'con_e', 'con_f']].agg(' '.join, axis=1)
        # Strip trailing whitespace in CT column
        self.con_df['CT'] = self.con_df['CT'].map(lambda x: x.rstrip())
        #
        return self
    #
    # TODO: Remove the atom indices you think you should
    # def remove():
    #     return self
    #
    def write(self, write_dir=None, write_file=None):
        """
        Write the LICHEM connectivity file.

        Parameters
        ----------
        write_dir : str
            Output directory.
        write_file : str
            Name of the file to write to.
        """
        if write_dir is None:
            write_dir = self.out_dir
        if write_file is None:
            write_file = self.outfile
        # Specify output path
        self.out_path = write_dir + '/' + write_file
        # Create the output directory (no error if it doesn't exist)
        pathlib.Path(write_dir).mkdir(parents=True, exist_ok=True)
        # Write to file
        # If n_connections is not 0, put a space before connections
        with open(self.out_path, 'w+') as f:
            for r in self.con_df.itertuples(index=True, name='Pandas'):
                f.write("{:} {:} {:} {:.2f} {:.4f} {:}{:}\n".format(\
                r.LICHEM_ID, r.atom_name, r.tinker_type, r.mass, r.charge, \
                r.n_connections,
                (' '+r.CT if r.CT != '' else '')))
        return


class LICHEMRegionsFile:
    """
    Parameters from a LICHEM "regions.inp" file.

    Methods
    -------
    read(read_file="regions.inp")
        Parse the regions file and set initial values.
    use_dfp(DFP_criteria="tight")
        Set parameters for DFP optimization.
    use_qsm(restrain_QSM=True, beads=3)
        Set parameters for restrained or unrestrained QSM.
    write(write_dir="./output", write_file="regions.inp")
        Write out an updated parameter file.
    """
    def __init__(self, regfile="regions.inp", out_dir="./output",
                 outfile="regions.inp"):
        self.regions_file = regfile
        self.out_dir = out_dir
        self.outfile = outfile
        return
    def read(self, read_file=None):
        '''
        Parse the LICHEM regions file.

        Parameters
        ----------
        read_file : str
            The file name of the LICHEM regions file.
        '''
        if read_file is None:
            read_file = self.regions_file
        # Set up option to append value
        def append_val(line):
            """Save keyword option specified on one-line."""
            junk, var = line.split(":")
            var = var.strip()
            return var
        #
        f = open(read_file, 'r')
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
                              f"because {self.electrostatics} was given "
                              "for electrostatics.")
                        self.lrec_exponent = 2
                    elif self.electrostatics.lower() == "amoeba":
                        print("\nNOTE: Setting LREC_exponent to 3 "
                              f"because {self.electrostatics} was given "
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
        return self
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
    def write(self, write_dir=None, write_file=None):
        """
        Write the LICHEM regions file.

        Parameters
        ----------
        write_dir : str
            Output directory.
        write_file : str
            Name of the file to write to.
        """
        if write_dir is None:
            write_dir = self.out_dir
        if write_file is None:
            write_file = self.outfile
        # Specify output path
        self.out_path = write_dir + '/' + write_file
        #
        # Create the output directory (no error if it doesn't exist)
        pathlib.Path(write_dir).mkdir(parents=True, exist_ok=True)
        #
        def check_attr(ofile_name, attribute):
            """Check whether the attribute is assigned."""
            if hasattr(self, attribute):
                f.write(f"{ofile_name}: {getattr(self, attribute)}\n")
            return
        # Write the update regions file
        with open(self.out_path, 'w+') as f:
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
                    for i, item in enumerate(self.neb_atoms, start=1):
                        # If i is divisible by 10 or last item
                        if (i) % 10 == 0 or item == self.neb_atoms[-1]:
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
            for i, item in enumerate(self.qm_atoms, start=1):
                # If i is divisible by 10 or last item
                if (i) % 10 == 0 or item == self.qm_atoms[-1]:
                    f.write(f"{item}\n")
                else:
                    f.write(f"{item} ")
            # Write PB atoms
            f.write(f"Pseudobond_atoms: {self.total_PB}\n")
            for i, item in enumerate(self.pseudobond_atoms, start=1):
                # If i is divisible by 10 or last item
                if (i) % 10 == 0 or item == self.pseudobond_atoms[-1]:
                    f.write(f"{item}\n")
                else:
                    f.write(f"{item} ")
            # Write BA
            f.write(f"Boundary_atoms: {self.total_boundary}\n")
            for i, item in enumerate(self.boundary_atoms, start=1):
                # If i is divisible by 10 or last item
                if (i) % 10 == 0 or item == self.boundary_atoms[-1]:
                    f.write(f"{item}\n")
                else:
                    f.write(f"{item} ")
            # Write frozen
            f.write(f"Frozen_atoms: {self.total_frozen}\n")
            for i, item in enumerate(self.frozen_atoms, start=1):
                # If i is divisible by 10 or last item
                if (i) % 10 == 0 or item == self.frozen_atoms[-1]:
                    f.write(f"{item}\n")
                else:
                    f.write(f"{item} ")
            f.close()
            return


class LICHEMXYZFile:
    """
    Information from a LICHEM XYZ file.

    Methods
    -------
    read(read_file="xyzfile.xyz")
        Read in the XYZ file.
    write(write_dir="./output", write_file="xyzfile.xyz")
        Write out the XYZ file.
    """
    # TODO: Only read last frame!
    def __init__(self, xyzfile="xyzfile.xyz", out_dir="./output",
                 outfile="xyzfile.xyz"):
        self.xyzfile = xyzfile
        self.out_dir = out_dir
        self.outfile = outfile
        return
    def read(self, read_file=None):
        """
        Parse the LICHEM XYZ file.
        """
        if read_file is None:
            read_file = self.xyzfile
        # Read firstline for number of atoms
        self.natom = open(read_file, 'r+').readlines()[0].strip()
        # Read in the XYZ file
        col_names = ["element", "x_coord", "y_coord", "z_coord"]
        self.xyz_df = pd.read_csv(
                read_file,
                low_memory=False,
                names = col_names,
                delim_whitespace=True, na_filter=False, skiprows=1)
        # Strip trailing whitespace in element column
        self.xyz_df['element'] = self.xyz_df['element'].map(lambda x: \
                                                            x.rstrip())
        # Save all information to lists
        self.elements = self.xyz_df['element'].to_list()
        self.x_coord = self.xyz_df['x_coord'].to_list()
        self.y_coord = self.xyz_df['y_coord'].to_list()
        self.z_coord = self.xyz_df['z_coord'].to_list()
        # print(self.xyz_df)
        return self
    #
    def write(self, write_dir=None, write_file=None):
        """
        Write the LICHEM XYZ file.

        Parameters
        ----------
        write_dir : str
            Output directory.
        write_file : str
            Name of the file to write to.
        """
        if write_dir is None:
            write_dir = self.out_dir
        if write_file is None:
            write_file = self.outfile
        # Specify output path
        self.out_path = write_dir + '/' + write_file
        # Create the output directory (no error if it doesn't exist)
        pathlib.Path(write_dir).mkdir(parents=True, exist_ok=True)
        #
        with open(self.out_path, 'w+') as f:
            # Write the number of atoms
            f.write(f"{self.natom}\n\n")
            for x in self.xyz_df.values:
                f.write(f"{x[0]:<3} "
                        f"{x[1]:0<6.3f} {x[2]:0<6.3f} {x[3]:0<6.3f}\n")
                # This will trailing zeros (but precision different than orig)
                # f.write(f"{x[0]:<3} "
                #          "{x[1]:0<16f} {x[2]:0<16f} {x[3]:0<16f}\n")
        return


class LICHEMSystem:
    """
    Information for running LICHEM with Gaussian/TINKER.

    Methods
    -------
    remove_atoms(mfr)
        Remove atom(s) from the system.
    write()
        Write all of the files.
    """
    def __init__(self, basis, connectivity, regions, xyz):
        """
        Create a dataframe of all the LICHEM information.

        Parameters:
        basis : GaussianGENFile
        connectivity : LICHEMConnectFile
        regions : LICHEMRegionsFile
        xyz : LICHEMXYZFile
        """
        self.basis = basis
        self.connectivity = connectivity
        self.regions = regions
        self.xyz = xyz
        #
        self.all_atom_df = self.connectivity.con_df
        #
        # Get info from XYZ
        self.natom = self.xyz.natom
        self.all_atom_df["element"] = self.xyz.elements
        self.all_atom_df["x_coord"] = self.xyz.x_coord
        self.all_atom_df["y_coord"] = self.xyz.y_coord
        self.all_atom_df["z_coord"] = self.xyz.z_coord
        #
        # Create a column if index matches an atom in the QM list
        self.all_atom_df['QM'] = self.all_atom_df['LICHEM_ID'].astype(
            str).isin(self.regions.qm_atoms)
        # Create a column if index matches an atom in the QM list
        self.all_atom_df['PB'] = self.all_atom_df['LICHEM_ID'].astype(
            str).isin(self.regions.pseudobond_atoms)
        # Create a column if index matches an atom in the QM list
        self.all_atom_df['Boundary'] = self.all_atom_df['LICHEM_ID'].astype(
            str).isin(self.regions.boundary_atoms)
        # Create a column if index matches an atom in the QM list
        self.all_atom_df['Frozen'] = self.all_atom_df['LICHEM_ID'].astype(
            str).isin(self.regions.frozen_atoms)
        # Create a column if index matches an atom in the QM list
        self.all_atom_df['NEB'] = self.all_atom_df['LICHEM_ID'].astype(
            str).isin(self.regions.neb_atoms)
        #
        # Create a set of all qm_pb_atoms
        qm_pb_atoms = set(self.regions.qm_atoms + \
                          self.regions.pseudobond_atoms)
        # Make them integers, then sort
        qm_pb_atoms = [int(x) for x in qm_pb_atoms]
        qm_pb_atoms = sorted(qm_pb_atoms)
        #
        if len(self.basis.basis_atoms) == len(qm_pb_atoms):
            self.all_atom_df['BASIS'] = self.all_atom_df['LICHEM_ID'].isin(
            qm_pb_atoms)
            self.all_atom_df['Func'] = None
            for i, atom in enumerate(qm_pb_atoms, start=1):
                # BASIS counting starts at 1 and keys are strings
                self.all_atom_df.loc[atom, 'Func'] = \
                    self.basis.basis_atoms[str(i)]
        else:
            print(f"\nWARNING: Information for {self.basis.in_basis_file} "
                  "does not match\n"
                  f"         the {self.regions.in_regions_file} file.\n"
                  f"         {self.basis.in_basis_file} will NOT be written "
                  f"correctly!")
        #
        #
        # Test that it worked
        # print(self.all_atom_df.loc[[0]])
        # print(self.all_atom_df.loc[[1234]])
        return
    #
    def remove_atoms(self, mfr):
        """
        Remove an atom from the overall LICHEM system.

        Parameters
        ----------
        mfr : int or list(int)
            LICHEM_ID of atom(s) marked for removal from the system.
            The ID is based on the original input, not updated values in the
            string. (So to remove 1, 2, 3, mfr = [1, 2, 3] and not [1, 1, 1]).
        """
        # Sorted list of atoms to remove
        self.mfr = sorted(list(mfr))
        #
        # Update connectivity info!
        self.orig_df = self.all_atom_df
        self.updated_df = self.all_atom_df
        #
        # print("Converting....")
        # Insert NaNs, but allows them to be updated...
        self.updated_df['n_connections'] = self.updated_df[
            'n_connections'].apply(
            pd.to_numeric, errors='coerce')
        self.updated_df['con_a'] = self.updated_df['con_a'].apply(
            pd.to_numeric, errors='coerce')
        self.updated_df['con_b'] = self.updated_df['con_b'].apply(
            pd.to_numeric, errors='coerce')
        self.updated_df['con_c'] = self.updated_df['con_c'].apply(
            pd.to_numeric, errors='coerce')
        self.updated_df['con_d'] = self.updated_df['con_d'].apply(
            pd.to_numeric, errors='coerce')
        self.updated_df['con_e'] = self.updated_df['con_e'].apply(
            pd.to_numeric, errors='coerce')
        self.updated_df['con_f'] = self.updated_df['con_f'].apply(
            pd.to_numeric, errors='coerce')
        #
        # print("Done converting....")
        #
        for i, atom in enumerate(self.mfr):
            # Break loop if input is non-integer
            if type(atom) is not int:
                print(f"\nWARNING: Non-integer value {atom} given for "
                      f"removal. Skipping this ID.\n"
                      f"         Maybe you meant to remove LICHEM_ID "
                      f"{int(atom)}?")
                break
            # Since atom to pull changes as you remove rows, subtract by
            #  increment. (i starts at 0, so 1st atom in list is unaffected)
            print(f"\nRemoving atom with LICHEM_ID {atom} from system.")
            atom -= i
            # Figure out rows after the atom that's being changed
            row_list = self.updated_df[atom+1:]
            row_list = row_list["LICHEM_ID"].to_list()
            # Test that rows are right
            # print(f"Rows: {row_list}")
            # Update index column by 1 for rows following atom
            self.updated_df.loc[row_list, 'LICHEM_ID'] = self.updated_df.loc[
                row_list, 'LICHEM_ID'] - 1
            #
            # Dict is a faster way to iterate...
            df_dict = self.updated_df.to_dict('records')
            #
            # for x in row_list:
            for row in df_dict:
                # Update N connections
                if row['con_a'] == atom:
                    new_key = row['n_connections'] - 1
                    row['n_connections'] = new_key
                    # Set to nothing
                    row['con_a'] = np.nan
                if row['con_b'] == atom:
                    new_key = row['n_connections'] - 1
                    row['n_connections'] = new_key
                    row['con_b'] = np.nan
                if row['con_c'] == atom:
                    new_key = row['n_connections'] - 1
                    row['n_connections'] = new_key
                    row['con_c'] = np.nan
                if row['con_d'] == atom:
                    new_key = int(row['n_connections']) - 1
                    row['n_connections'] = new_key
                    row['con_d'] = np.nan
                if row['con_e'] == atom:
                    new_key = int(row['n_connections']) - 1
                    row['n_connections'] = new_key
                    row['con_e'] = np.nan
                if row['con_f'] == atom:
                    new_key = row['n_connections'] - 1
                    row['n_connections'] = new_key
                    row['con_f'] = np.nan
            #
            # Convert back to data frame
            self.updated_df = pd.DataFrame.from_dict(df_dict)
            #
            # Lookup where atom exists, then decrease by 1
            for col in ['con_a', 'con_b', 'con_c', 'con_d', 'con_e',
                        'con_f']:
                self.updated_df[col] = np.where((
                    # If column value greater than atom of interest
                    self.updated_df[col] >= atom)
                    & (self.updated_df[col] != np.nan),
                    # This is the value that is inserted
                    self.updated_df[col] - 1,
                    # This is the column that is changed
                    self.updated_df[col])
            #
            # Drop the atom row
            self.updated_df.drop(self.updated_df[self.updated_df[
                'LICHEM_ID'] == atom].LICHEM_ID, axis='index', inplace=True)
            # Reset the index -- use inplace to avoid NaNs in dropped row
            self.updated_df.reset_index(inplace=True)
            # Make new index the LICHEM_ID column
            self.updated_df['LICHEM_ID'] = \
                self.updated_df.index.values.tolist()
        #
        # End each atom
        #
        # Update NaNs
        for col in ['con_a', 'con_b', 'con_c', 'con_d', 'con_e', 'con_f']:
            # Fill NaNs with -1
            self.updated_df[col] = self.updated_df[col].fillna(-1)
            # Change to int ('Int64' does not help...)
            self.updated_df[col] = self.updated_df[col].astype(int)
            # Convert to string
            self.updated_df[col] = self.updated_df[col].astype(str)
            # Replace the -1 NaNs with ''
            self.updated_df[col] = self.updated_df[col].replace('-1', '')
        #
        # Update CT
        # Combine all con columns into a one-column space-separated list
        self.updated_df['CT'] = self.updated_df[['con_a', 'con_b', 'con_c',
                        'con_d', 'con_e', 'con_f']].agg(' '.join, axis=1)
        #
        # Strip whitespace in CT column since not all con columns have values
        #  and spaces get introduced as the separator
        self.updated_df['CT'] = self.updated_df['CT'].map(lambda x: x.strip())
        #
        # Now that update is complete, reset
        self.all_atom_df = self.updated_df
        #
        # Connectivity Updates
        self.connectivity.con_df = self.all_atom_df[[
            'LICHEM_ID', 'atom_name', 'tinker_type', 'mass', 'charge',
            'n_connections', 'con_a', 'con_b', 'con_c', 'con_d', 'con_e',
            'con_f', 'CT']]
        #
        # BASIS updates
        # Find where BASIS is true
        basis_indicies = list(self.all_atom_df['BASIS'].index[
            self.all_atom_df['BASIS']])
        # Recreate dictionary
        self.basis.basis_atoms = {}
        # BASIS file starts at 1
        for i, row in enumerate(basis_indicies, start=1):
            self.basis.basis_atoms[str(i)] = self.all_atom_df.loc[row, 'Func']
        #
        # Verify that it looks right
        # print(self.basis.basis_atoms)
        #
        # Regions updates
        self.regions.qm_atoms = list(self.all_atom_df['QM'].index[
            self.all_atom_df['QM']])
        self.regions.pseudobond_atoms = list(self.all_atom_df['PB'].index[
            self.all_atom_df['PB']])
        self.regions.boundary_atoms = list(self.all_atom_df['Boundary'].index[
            self.all_atom_df['Boundary']])
        self.regions.frozen_atoms = list(self.all_atom_df['Frozen'].index[
            self.all_atom_df['Frozen']])
        self.regions.neb_atoms = list(self.all_atom_df['NEB'].index[
            self.all_atom_df['NEB']])
        #
        # Test that it worked
        # print(self.all_atom_df.loc[[0]])
        # print(self.all_atom_df.loc[[1234]])
        #
        # Update numbers
        self.regions.total_QM = len(self.regions.qm_atoms)
        self.regions.total_PB = len(self.regions.pseudobond_atoms)
        self.regions.total_boundary = len(self.regions.boundary_atoms)
        self.regions.total_frozen = len(self.regions.frozen_atoms)
        self.regions.total_NEB = len(self.regions.neb_atoms)
        #
        # XYZ updates
        self.natom = len(self.all_atom_df)
        self.xyz.natom = self.natom
        # Only relevant columns
        self.xyz.xyz_df = self.all_atom_df[
            ['element', 'x_coord', 'y_coord', 'z_coord']]
        # Get elements/coords
        self.xyz.elements = self.all_atom_df["element"].to_list()
        self.xyz.x_coord = self.all_atom_df["x_coord"].to_list()
        self.xyz.y_coord = self.all_atom_df["y_coord"].to_list()
        self.xyz.z_coord = self.all_atom_df["z_coord"].to_list()
        #
        return self
    #
    def write(self, out_dir="./output"):
        """
        Write all of the LICHEM system files.
        """
        self.out_dir = out_dir
        # Rewrite the files
        self.basis.write(write_dir=out_dir)
        self.connectivity.write(write_dir=out_dir)
        self.regions.write(write_dir=out_dir)
        self.xyz.write(write_dir=out_dir)
        return

# TODO: Create new Tinker XYZ?
# TODO: Update BASIS_verification.txt?

# ------------------------------------------------------------------------ #

# Examples:
# basis = GaussianGENFile(outfile="BASIS").read()
# basis.write(write_dir="./funky", write_file="WHY")
# connectivity = LICHEMConnectFile()
# c = connectivity.read()
# c.write()
# w = LICHEMSystem(basis, connectivity, regions, xyz).remove_atoms([2, 3])
# w.write(out_dir="./updated-files")

# start = time.time()

# Read and write the BASIS file
basis = GaussianGENFile().read()
# basis.write()

# Read and write the connect file
connectivity = LICHEMConnectFile().read()
# connectivity.write()

# Read and write the regions file
regions = LICHEMRegionsFile().read()
# regions.write()

# Read and write the xyzfile
xyz = LICHEMXYZFile().read()
# xyz.write()

# List of atoms to remove (using LICHEM_ID!!!!)
remove_me = [1234]

# Create a LICHEM system, remove atoms, and then write the new output.
w = LICHEMSystem(basis, connectivity, regions, xyz).remove_atoms(remove_me)
w.write()

# end = time.time()
# total_time = end - start
# print(f"\nTotal time: {total_time}")
