import parmed as pmd
import pandas as pd

## Code to source a single parm file (not a leaprc)
# source_params = "parm99.dat"
# param_dat = pmd.load_file(source_params)

## It looks like anything sourced in the leaprc needs to have an absolute
## path to it, so you might need to modify the leaprc to incorporate the
## absolute path. It is STRONGLY RECOMMENDED that you copy the leaprc file to
## do this, and then reference that copy!!

source_params = "param_files/leaprc.ff14SB.OL15.tip3p"
param_dat = pmd.amber.AmberParameterSet().from_leaprc(source_params)

## Give your FF a name (it will be preceded by AMBER-)
ff_name = "ff14SB"

##################
##  Definitions ##
##################

def get_ATs(param_dat):
    """Create a dictionary for the atom types, where the values are the TINKER
    atom types.
    """
    ## If you use enumerate, then the keys are numeric, which doesn't help.
    atom_types = {}
    i = 1
    for key in param_dat.atom_types.keys():
        atom_types[key] = i
        i += 1
    ## Add 'X' atom type
    print(" If you encounter analyze errors, 999 is the numeric X placeholder.\n",
    "Just add in the missing parameters that would use those.")
    atom_types['X'] = 999
    return atom_types

def get_VDW(param_dat):
    """Create dictionaries for the rmin_14 and epsilon_14 based on the atom
    types.
    """
    at_test = list(param_dat.atom_types.items())
    rmin14_dict = {}
    eps14_dict = {}
    for i in range(len(at_test)):
        rmin14_dict[at_test[i][1].name] = at_test[i][1].rmin_14
        eps14_dict[at_test[i][1].name] = at_test[i][1].epsilon_14
    return rmin14_dict, eps14_dict

def get_bonds(param_dat):
    """Create lists of the bond terms, k values, and req.
    """
    bonds = list(param_dat.bond_types.items())
    bond1a = []
    bond2a = []
    bond_k = []
    bond_req = []
    for i in range(len(bonds)):
        bond1a.append(bonds[i][0][0])
        bond2a.append(bonds[i][0][1])
        bond_k.append(bonds[i][1].k)
        bond_req.append(bonds[i][1].req)
    bond1 = [atom_types[key] for key in bond1a]
    bond2 = [atom_types[key] for key in bond2a]
    return bond1, bond2, bond_k, bond_req

def get_angles(param_dat):
    """Create lists of the angle terms, force constant values (k), and
    equilibrium angles (theteq).
    """
    angles = list(param_dat.angle_types.items())
    angle1a = []
    angle2a = []
    angle3a = []
    ang_k = []
    ang_theteq = []
    for i in range(len(angles)):
        angle1a.append(angles[i][0][0])
        angle2a.append(angles[i][0][1])
        angle3a.append(angles[i][0][2])
        ang_k.append(angles[i][1].k)
        ang_theteq.append(angles[i][1].theteq)
    angle1 = [atom_types[key] for key in angle1a]
    angle2 = [atom_types[key] for key in angle2a]
    angle3 = [atom_types[key] for key in angle3a]
    return angle1, angle2, angle3, ang_k, ang_theteq

def get_dihedrals(param_dat):
    """Create lists of the dihedral terms and extract out the force constant
    (k), periodicity (per), and phase information.
    """
    dihedrals = list(param_dat.dihedral_types.items())
    dihedral1a = []
    dihedral2a = []
    dihedral3a = []
    dihedral4a = []
    dval = []
    for i in range(len(dihedrals)):
        dihedral1a.append(dihedrals[i][0][0])
        dihedral2a.append(dihedrals[i][0][1])
        dihedral3a.append(dihedrals[i][0][2])
        dihedral4a.append(dihedrals[i][0][3])
        #
        di_type = str(dihedrals[i][1])
        di_type = str(di_type)
        di_type.replace('<','').replace('>','').replace('[','').replace(']','').replace('DihedralTypes','').replace('DihedralType;','')
        w_di_type = di_type.split()
        for i in range(len(w_di_type)):
            w_di_type[i] = w_di_type[i].strip(',')
        #
        saved_k = []
        saved_per = []
        saved_phase = []
        for i in range(len(w_di_type)):
            if w_di_type[i].rstrip('=.1234567890') == 'phi_k':
                saved_k.append(w_di_type[i].lstrip('phi_k='))
            if w_di_type[i].rstrip('=.1234567890') == 'per':
                saved_per.append(w_di_type[i].lstrip('per='))
            if w_di_type[i].rstrip('=.1234567890') == 'phase':
                saved_phase.append(w_di_type[i].lstrip('phase='))
        #
        dval.append(list(zip(saved_k, saved_per, saved_phase)))
        #
    dihedral1 = [atom_types[key] for key in dihedral1a]
    dihedral2 = [atom_types[key] for key in dihedral2a]
    dihedral3 = [atom_types[key] for key in dihedral3a]
    dihedral4 = [atom_types[key] for key in dihedral4a]
    return dihedral1, dihedral2, dihedral3, dihedral4, dval

def clean_dihedrals(dval):
    """Create the print string for the force constant (k), periodicity (per),
    and phase.
    """
    di_line = []
    for i in dval:
        i = str(i)
        di_line.append(i.replace(')','').replace('(','').replace(',','').replace('\'','').replace('[','').replace(']',''))
    return di_line

def get_improper_torsions(param_dat):
    """Create lists of the improper terms, force constant values (phi_k), and
    equilibrium angles (theteq).
    """
    imptors = list(param_dat.improper_periodic_types.items())
    imptor1a = []
    imptor2a = []
    imptor3a = []
    imptor4a = []
    imp_diheds_phik = []
    imp_diheds_phase = []
    imp_diheds_per = []
    for i in range(len(imptors)):
        imptor1a.append(imptors[i][0][0])
        imptor2a.append(imptors[i][0][1])
        imptor3a.append(imptors[i][0][2])
        imptor4a.append(imptors[i][0][3])
        imp_diheds_phik.append(imptors[i][1].phi_k)
        imp_diheds_phase.append(imptors[i][1].phase)
        imp_diheds_per.append(imptors[i][1].per)
    imptor1 = [atom_types[key] for key in imptor1a]
    imptor2 = [atom_types[key] for key in imptor2a]
    imptor3 = [atom_types[key] for key in imptor3a]
    imptor4 = [atom_types[key] for key in imptor4a]
    return imptor1, imptor2, imptor3, imptor4, imp_diheds_phik,\
     imp_diheds_phase, imp_diheds_per

# def get_atomlist(param_dat, rmin14_dict, eps14_dict):
def get_atomlist(param_dat):
    """Create a dataframe of all the different residues and update element and
    mass.
    """
    res_types = list(param_dat.residues.items())
    #
    res_df = pd.DataFrame()
    for i in range(len(res_types)):
        test_df = res_types[i][1].to_dataframe()
        res_df = res_df.append(test_df, ignore_index=True, sort=False)
    #
    ## Add element column that uses atomic_num to get element
    res_df['element'] = ""
    ## Update mass based on element
    for i in range(len(res_df)):
        res_df.loc[i,'element'] = list(pmd.periodic_table.AtomicNum.keys())[list(pmd.periodic_table.AtomicNum.values()).index(res_df.loc[i, 'atomic_number'])]
        res_df.loc[i,'mass'] = pmd.periodic_table.Mass[res_df.loc[i,'element']]
        # res_df.loc[i,'rmin_14'] = rmin14_dict[res_df.loc[i,'type']]
        # res_df.loc[i,'epsilon_14'] = eps14_dict[res_df.loc[i,'type']]
    #
    # ## Update index to start at 1, not zero (for mapping)
    # res_df.index = res_df.index + 1
    return res_df

def write_params(param_dat, atom_types, bond1, bond2, bond_k, bond_req,\
 angle1, angle2, angle3, ang_k, ang_theteq,\
 imptor1, imptor2, imptor3, imptor4, imp_diheds_phik, imp_diheds_phase,\
 imp_diheds_per, dihedral1, dihedral2, dihedral3, dihedral4, di_line, res_df,\
 rmin14_dict, eps14_dict, ff_name):
    """Generate the new TINKER parameter file.
    """
    with open("TINKER_params.prm", "w+") as tp_out:
        tp_out.write("\n")
        tp_out.write("      ##############################\n")
        tp_out.write("      ##                          ##\n")
        tp_out.write("      ##  Force Field Definition  ##\n")
        tp_out.write("      ##                          ##\n")
        tp_out.write("      ##############################\n")
        tp_out.write("\n\n")
        tp_out.write("forcefield              AMBER-{}\n\n".format(ff_name))
        tp_out.write("vdwtype                 LENNARD-JONES\n")
        tp_out.write("radiusrule              ARITHMETIC\n")
        tp_out.write("radiustype              R-MIN\n")
        tp_out.write("radiussize              RADIUS\n")
        tp_out.write("epsilonrule             GEOMETRIC\n")
        tp_out.write("vdw-14-scale            {}\n".format(param_dat.default_scnb))
        tp_out.write("chg-14-scale            {}\n".format(param_dat.default_scee))
        tp_out.write("electric                {:.7f}\n".format(pmd.constants.AMBER_ELECTROSTATIC**2))
        tp_out.write("dielectric              1.0\n")
        tp_out.write("\n\n")
        tp_out.write("      #############################\n")
        tp_out.write("      ##                         ##\n")
        tp_out.write("      ##  Literature References  ##\n")
        tp_out.write("      ##                         ##\n")
        tp_out.write("      #############################\n")
        tp_out.write("\n\n")
        tp_out.write("\nThis was generated using parmed from the following parameter sets:\n")
        for name in param_dat.titles:
            tp_out.write("{}.\n".format(name))
        tp_out.write("\n")
        tp_out.write("Current parameter values are available from the Amber site, located\nat http://ambermd.org/\n")
        tp_out.write("\n\n")
        tp_out.write("   ##################################\n")
        tp_out.write("   ##                              ##\n")
        tp_out.write("   ##  Tinker Atom Class Numbers   ##\n")
        tp_out.write("   ##      to Amber Atom Types     ##\n")
        tp_out.write("   ##                              ##\n")
        for atom in atom_types.items():
            tp_out.write("   ##           {:3}  {:<4}          ##\n".format(atom[1], atom[0]))
        tp_out.write("   ##                              ##\n")
        tp_out.write("   ##################################\n")
        tp_out.write("\n\n")
        tp_out.write("      #############################\n")
        tp_out.write("      ##                         ##\n")
        tp_out.write("      ##  Atom Type Definitions  ##\n")
        tp_out.write("      ##                         ##\n")
        tp_out.write("      #############################\n")
        tp_out.write("\n\n")
        for i in range(len(res_df)):
            if res_df.loc[i,'type'] in ('Ag+', 'Cu+'):
                tp_out.write("atom       {:4}  {:2}    {:<2}    \"{:<30} {:3}   {:>7.2f}\n".format((i+1),\
                 atom_types[res_df.loc[i,'type'][:-1]+"2+"], res_df.loc[i,'type'],\
                 (res_df.loc[i, 'resname'] + " " + res_df.loc[i, 'name']+'\"'),\
                 res_df.loc[i, 'atomic_number'], res_df.loc[i, 'mass']))
            else:
                tp_out.write("atom       {:4}  {:2}    {:<2}    \"{:<30} {:3}   {:>7.2f}\n".format((i+1),\
                 atom_types[res_df.loc[i,'type']], res_df.loc[i,'type'],\
                 (res_df.loc[i, 'resname'] + " " + res_df.loc[i, 'name']+'\"'),\
                 res_df.loc[i, 'atomic_number'], res_df.loc[i, 'mass']))
        tp_out.write("\n\n")
        tp_out.write("      ################################\n")
        tp_out.write("      ##                            ##\n")
        tp_out.write("      ##  Van der Waals Parameters  ##\n")
        tp_out.write("      ##                            ##\n")
        tp_out.write("      ################################\n")
        tp_out.write("\n\n")
        for i, (value1, value2) in enumerate(zip(rmin14_dict.items(), eps14_dict.items())):
            ## Get around immutable tuple by setting as a list
            value1 = list(value1)
            value2 = list(value2)
            if value1[1] is None:
                value1[1] = -0.0000
                print("WARNING! {} has no listed rmin_14 value (VDW). Setting to -0.0000.".format(value1[0]))
            if value2[1] is None:
                value2[1] = -0.000
                print("WARNING! {} has no listed epsilon_14 value (VDW). Setting to -0.0000.".format(value2[0]))
            tp_out.write("vdw          {:2}               {:>8.4f}     {:>12.7f}\n".format(i, value1[1], value2[1]))
        tp_out.write("\n\n")
        tp_out.write("      ##################################\n")
        tp_out.write("      ##                              ##\n")
        tp_out.write("      ##  Bond Stretching Parameters  ##\n")
        tp_out.write("      ##                              ##\n")
        tp_out.write("      ##################################\n")
        tp_out.write("\n\n")
        for i in range(len(bond1)):
            tp_out.write("bond        {:3}  {:3}          {:3.2f}     {:1.4f}\n".format(bond1[i],\
             bond2[i], bond_k[i], bond_req[i]))
        tp_out.write("\n\n")
        tp_out.write("      ################################\n")
        tp_out.write("      ##                            ##\n")
        tp_out.write("      ##  Angle Bending Parameters  ##\n")
        tp_out.write("      ##                            ##\n")
        tp_out.write("      ################################\n")
        tp_out.write("\n\n")
        ## You might be able to loop through for any of the 999 choices!
        for i in range(len(angle1)):
            tp_out.write("angle        {:2}   {:2}   {:2}     {:>5.1f}     {:>6.2f}\n".format(angle1[i],\
             angle2[i], angle3[i], ang_k[i], ang_theteq[i], sep=' '))
        tp_out.write("\n\n")
        tp_out.write("      #####################################\n")
        tp_out.write("      ##                                 ##\n")
        tp_out.write("      ##  Improper Torsional Parameters  ##\n")
        tp_out.write("      ##                                 ##\n")
        tp_out.write("      #####################################\n")
        tp_out.write("\n\n")
        for i in range(len(imptor1)):
            tp_out.write("imptors      {:3}  {:3}  {:3}  {:3}           {:6}  {:5}  {:1}\n".format(imptor1[i],\
             imptor2[i], imptor3[i], imptor4[i], imp_diheds_phik[i], imp_diheds_phase[i], imp_diheds_per[i]))
        tp_out.write("\n\n")
        tp_out.write("            ############################\n")
        tp_out.write("            ##                        ##\n")
        tp_out.write("            ##  Torsional Parameters  ##\n")
        tp_out.write("            ##                        ##\n")
        tp_out.write("            ############################\n")
        tp_out.write("\n\n")
        for i in range(len(dihedral1)):
            tp_out.write("torsion      {:3}  {:3}  {:3}  {:3}           {:6}\n".format(dihedral1[i],\
             dihedral2[i], dihedral3[i], dihedral4[i], di_line[i]))
        tp_out.write("\n\n")
        tp_out.write("      ########################################\n")
        tp_out.write("      ##                                    ##\n")
        tp_out.write("      ##  Atomic Partial Charge Parameters  ##\n")
        tp_out.write("      ##                                    ##\n")
        tp_out.write("      ########################################\n")
        tp_out.write("\n\n")
        for i in range(len(res_df)):
            tp_out.write("charge     {:>4}              {:>8.4f}\n".format((i+1), res_df.loc[i,'charge']))
        tp_out.write("\n\n")
        ## You might need biotype, but I don't know if it's actually used.
        tp_out.close()

###########
##  RUN  ##
###########
atom_types = get_ATs(param_dat)

rmin14_dict, eps14_dict = get_VDW(param_dat)

bond1, bond2, bond_k, bond_req = get_bonds(param_dat)

angle1, angle2, angle3, ang_k, ang_theteq = get_angles(param_dat)

dihedral1, dihedral2, dihedral3, dihedral4, dval = get_dihedrals(param_dat)

di_line = clean_dihedrals(dval)

imptor1, imptor2, imptor3, imptor4, imp_diheds_phik, imp_diheds_phase,\
 imp_diheds_per = get_improper_torsions(param_dat)

# res_df = get_atomlist(param_dat, rmin14_dict, eps14_dict)
res_df = get_atomlist(param_dat)

write_params(param_dat, atom_types, bond1, bond2, bond_k, bond_req, \
 angle1, angle2, angle3, ang_k, ang_theteq,\
 imptor1, imptor2, imptor3, imptor4, imp_diheds_phik, imp_diheds_phase,\
 imp_diheds_per, dihedral1, dihedral2, dihedral3, dihedral4, di_line, res_df,\
 rmin14_dict, eps14_dict, ff_name)
