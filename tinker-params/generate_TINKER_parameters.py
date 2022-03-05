## !!! This uses Pandas 1.0 !!! Without it, remove 'ignore_index' in the
## drop_duplicates lines.

import parmed as pmd
import pandas as pd
import numpy as np
import copy
from collections import OrderedDict

## Code to source a single parm file (not a leaprc)
# source_params = "parm99.dat"
# param_dat = pmd.load_file(source_params)

## Prmtop Method
source_params = "my_amber_system.prmtop"
param_dat = pmd.load_file(source_params)

## It looks like anything sourced in the leaprc needs to have an absolute
## path to it, so you might need to modify the leaprc to incorporate the
## absolute path. It is STRONGLY RECOMMENDED that you copy the leaprc file to
## do this, and then reference that copy!!

## Leaprc Method
# source_params = "param_files/leaprc.ff14SB.OL15.tip3p"
# param_dat = pmd.amber.AmberParameterSet().from_leaprc(source_params)

## Leave the X dihedrals as X to fix by hand
## Setting this false will likely result in a MAXPRM issue with TINKER
## Because... well... you'll likely have over a million dihedral angles.
leave_as_X = True

## Give your FF a name (it will be preceded by AMBER-)
ff_name = "ff14SB"

## Give your new parameter file a name
param_file_name = "TINKER_params.prm"

##################
##  Definitions ##
##################

def det_structure(param_dat):
    if type(param_dat) == pmd.amber._amberparm.AmberParm:
        print(" Using an AMBER prmtop file.\n")
        ## Deal with C and N Terminals
        ## If there's a terminal OXT for the protein, use as CTERM
        for residue in param_dat.residues:
            for atom in residue.atoms:
                if atom.name == 'OXT':
                    residue.name = 'C'+residue.name
                else:
                    continue
            ## If there's an H3 in a protein residue, use NTERM
            if residue.name in ('ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'CYX', 'GLN',
             'GLU', 'GLY', 'HID', 'HIE', 'HIP', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
             'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'):
                for atom in residue.atoms:
                    if atom.name == 'H3':
                        residue.name = 'N'+residue.name
        #
        struct_dat = pmd.amber.AmberParameterSet.from_structure(param_dat)
    else:
        struct_dat = None
    return struct_dat

def get_ATs(param_dat, struct_dat):
    """Create a dictionary for the atom types, where the values are the TINKER
    atom types.
    """
    ## AMBER parmtops already have this, but it'll cause issues later...
    # if type(param_dat) == pmd.amber._amberparm.AmberParm:
    #     atom_types = copy.deepcopy(param_dat.LJ_types)
    ## If you use enumerate, then the keys are numeric, which doesn't help.
    atom_types = {}
    i = 1
    if struct_dat is None:
        for key in param_dat.atom_types.keys():
            atom_types[key] = i
            i += 1
    else:
        for key in struct_dat.atom_types.keys():
            atom_types[key] = i
            i += 1
    ## Add 'X' atom type
    print(" If you encounter analyze errors, 999 is the numeric X placeholder.\n",
    "Just add in the missing parameters that would use those.\n")
    atom_types['X'] = 999
    #
    return atom_types

def guess_connectivity(atom_types):
    """Guess the atom connectivity based on previous atom class assignments.
    """
    atom_connect = {
    ## CX is new for ff14SB, splits off from CA
    ## CO and C8 new for ff14SB; CO was C, C8 was CT
    ## 2C and 3C were CT
    ## C2 was CK, CJ was CB but TINKER complained about 3, C7 was CT, C1 was CM
    "CT": 4, "C" : 3, "CA": 3, "CM": 3, "CC": 3, "CV": 3, "CW": 3, "CR": 3,
    "CB": 3, "C*": 3, "CN": 3, "CK": 3, "CQ": 3, "CX": 4, "C8": 4, "CO": 3,
    "C2": 3, "CJ": 4, "C7": 4, "C1": 3,
    "2C": 4, "3C": 4,
    "N" : 3, "NA": 3, "NB": 2, "NC": 2, "N*": 3, "N2": 3, "N3": 4,
    "OW": 2, "OH": 2, "OS": 2, "O" : 1, "O2": 1, "S" : 2, "SH": 2, "P" : 4,
    "H" : 1, "HW": 1, "HO": 1, "HS": 1, "HA": 1, "HC": 1, "H1": 1, "H2": 1,
    "H3": 1, "HP": 1, "H4": 1, "H5": 1,
    ## Add common ions as zero
    "FE"  : 0, "Zn"  : 0, "Li+" : 0, "Na+" : 0, "K+"  : 0, "Rb+" : 0, "Cs+" : 0,
    "F-"  : 0, "Cl-" : 0, "Br-" : 0, "I-"  : 0, "Be2+": 0, "Cu2+": 0, "Ni2+": 0,
    "Pt2+": 0, "Zn2+": 0, "Co2+": 0, "Pd2+": 0, "Ag2+": 0, "Cr2+": 0, "Fe2+": 0,
    "Mg2+": 0, "V2+" : 0, "Mn2+": 0, "Hg2+": 0, "Cd2+": 0, "Yb2+": 0, "Ca2+": 0,
    "Sn2+": 0, "Pb2+": 0, "Eu2+": 0, "Sr2+": 0, "Sm2+": 0, "Ba2+": 0, "Ra2+": 0,
    "Al3+": 0, "Fe3+": 0, "Cr3+": 0, "In3+": 0, "Tl3+": 0, "Y3+" : 0, "La3+": 0,
    "Ce3+": 0, "Pr3+": 0, "Nd3+": 0, "Sm3+": 0, "Eu3+": 0, "Gd3+": 0, "Tb3+": 0,
    "Dy3+": 0, "Er3+": 0, "Tm3+": 0, "Lu3+": 0, "Hf4+": 0, "Zr4+": 0, "Ce4+": 0,
    "U4+" : 0, "Pu4+": 0, "Th4+": 0
    }
    return atom_connect

def get_VDW(param_dat, struct_dat):
    """Create dictionaries for the rmin_14 and epsilon_14 based on the atom
    types.
    """
    rmin14_dict = {}
    eps14_dict = {}
    # if type(param_dat) == pmd.amber._amberparm.AmberParm:
    #     at_test = list(param_dat.LJ_types)
    #     for i in range(len(at_test)):
    #         rmin14_dict[at_test[i]] = param_dat.LJ_radius[param_dat.LJ_types[at_test[i]]-1]
    #         eps14_dict[at_test[i]] = param_dat.LJ_depth[param_dat.LJ_types[at_test[i]]-1]
    if struct_dat is None:
        at_test = list(param_dat.atom_types.items())
    else:
        at_test = list(struct_dat.atom_types.items())
    for i in range(len(at_test)):
        rmin14_dict[at_test[i][1].name] = at_test[i][1].rmin_14
        eps14_dict[at_test[i][1].name] = at_test[i][1].epsilon_14
    #
    return rmin14_dict, eps14_dict

def get_bonds(param_dat, atom_types):
    """Create lists of the bond terms, k values, and req.
    """
    ## If an AMBER prmtop
    if type(param_dat) == pmd.amber._amberparm.AmberParm:
        origin = 1
        ## Create a dictionary of the bonds
        bo_OD = { 'columns' : ('B1', 'B2', 'k', 'req')}
        for i in range(len(param_dat.bonds)):
            bo_OD.update({origin: [ param_dat.bonds[i].atom1.type,\
             param_dat.bonds[i].atom2.type, param_dat.bonds[i].type.k,\
             param_dat.bonds[i].type.req ] })
            origin += 1
        init_df = pd.DataFrame.from_dict(bo_OD, "index", columns=['B1', 'B2',
         'k', 'req'])
        #
        # Remove the row used to set-up the dictionary
        bo_df = init_df.drop(['columns'])
        #
        # Drop duplicate values, keeping one copy
        bo_df.drop_duplicates(keep = 'first', inplace = True, ignore_index=True)
        #
        # Lookup AMBER types in dictionary and create column with TINKER type
        bo_df['B1T'] = bo_df['B1'].map(atom_types)
        bo_df['B2T'] = bo_df['B2'].map(atom_types)
        #
        # Write them to a list for printing
        bond1 = bo_df['B1T'].tolist()
        bond2 = bo_df['B2T'].tolist()
        bond_k = bo_df['k'].tolist()
        bond_req = bo_df['req'].tolist()
    else:
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
    print(" Achievement unlocked: the names bond, atom bond.\n")
    return bond1, bond2, bond_k, bond_req

def get_angles(param_dat, atom_types):
    """Create lists of the angle terms, force constant values (k), and
    equilibrium angles (theteq).
    """
    # ## If an AMBER prmtop
    # ## You can do it this way, but there's point need to...
    # if type(param_dat) == pmd.amber._amberparm.AmberParm:
    #     origin = 1
    #     an_OD = { 'columns' : ('A1', 'A2', 'A3', 'k', 'theteq')}
    #     for i in range(len(param_dat.angles)):
    #         an_OD.update({origin: [ param_dat.angles[i].atom1.type,\
    #          param_dat.angles[i].atom2.type, param_dat.angles[i].atom3.type,
    #          param_dat.angles[i].type.k, param_dat.angles[i].type.theteq ] })
    #         origin += 1
    #     init_df = pd.DataFrame.from_dict(an_OD, "index", columns=['A1', 'A2',
    #      'A3', 'k', 'theteq'])
    #     #
    #     # Remove the row used to set-up the dictionary
    #     an_df = init_df.drop(['columns'])
    #     #
    #     # Drop duplicate values, keeping one copy
    #     an_df.drop_duplicates(keep = 'first', inplace = True, ignore_index=True)
    #     #
    #     # Lookup AMBER types  in dictionary and create column with TINKER type
    #     an_df['A1T'] = an_df['A1'].map(atom_types)
    #     an_df['A2T'] = an_df['A2'].map(atom_types)
    #     an_df['A3T'] = an_df['A3'].map(atom_types)
    #     #
    #     # Write out to a list for printing
    #     angle1 = an_df['A1T'].tolist()
    #     angle2 = an_df['A2T'].tolist()
    #     angle3 = an_df['A3T'].tolist()
    #     ang_k = an_df['k'].tolist()
    #     ang_theteq = an_df['theteq'].tolist()
    if struct_dat is None:
        angles = list(param_dat.angle_types.items())
    else:
        angles = list(struct_dat.angle_types.items())
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
    print(" Achievement unlocked: angles managed.\n")
    return angle1, angle2, angle3, ang_k, ang_theteq

def get_dihedrals(param_dat, struct_dat, atom_types):
    """Create lists of the dihedral terms and extract out the force constant
    (k), periodicity (per), and phase information.
    """
    if struct_dat is None:
        dihedrals = list(param_dat.dihedral_types.items())

    else:
        dihedrals = list(struct_dat.dihedral_types.items())
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
        ## Combine k, periodicity, and phase into a single line
        # dval.append(list(zip(saved_k, saved_per, saved_phase)))
        ## Tinker does k, phase, and periodicity!!!!
        dval.append(list(zip(saved_k, saved_phase, saved_per)))
        #
    dihedral1 = [atom_types[key] for key in dihedral1a]
    dihedral2 = [atom_types[key] for key in dihedral2a]
    dihedral3 = [atom_types[key] for key in dihedral3a]
    dihedral4 = [atom_types[key] for key in dihedral4a]
    return dihedral1, dihedral2, dihedral3, dihedral4, dval

## Should figure out how to remove duplicates (X-A-B-X == X-B-A-X)
def clean_dihedrals(dval):
    """Create the print string for the force constant (k), periodicity (per),
    and phase.
    """
    di_line = []
    for i in dval:
        i = str(i)
        di_line.append(i.replace(')','').replace('(','').replace(',','').replace('\'','').replace('[','').replace(']',''))
    print(" Achievement unlocked: dihedrals dealt with.\n")
    return di_line


def build_X_dihedrals(dihedral1, dihedral2, dihedral3, dihedral4, atom_types, \
 di_line):
    ## Copy atoms types to get rid of ions from torsions
    c_atom_types = copy.deepcopy(atom_types)
    ## Remove know ion types AND X to not have repeats with X
    ## Each of the iteration lines will print the original X-type for posterity
    ions = [
    ## KEEP THESE
     #'HO', 'H1', 'CT', 'OH', 'OS', 'C7', 'CJ', 'H2', 'N*', 'C', 'O', 'NA',
     #'CM', 'H4', 'HC', 'P', 'H'
    ## These might not be important
     #'C', 'CA', 'CB', 'CC', 'CD', 'CI', 'CK', 'CP', 'CM', 'CS', 'CN', 'CQ',
     #'CR', 'CT', 'CV', 'CW', 'C*', 'CX', 'CY', 'CZ', 'C5', 'C4', 'C0', 'H',
     #'HC', 'H1', 'H2', 'H3', 'HA', 'H4', 'H5', 'HO', 'HS', 'HW', 'HP', 'HZ',
     #'N', 'NA', 'NB', 'NC', 'N2', 'N3', 'NT', 'N*', 'NY', 'O', 'O2', 'OW',
     #'OH', 'OS', 'OP', 'P', 'S', 'SH', 'CU',
    ## REMOVE REST
     'F', 'Cl', 'Br', 'I', 'MG',
     'FE', 'Zn', 'EP', 'CO', '2C', '3C', 'C8', 'C7', 'C2', 'C1', 'CJ',
     'Li+', 'Na+', 'K+', 'Rb+', 'Cs+', 'F-', 'Cl-', 'Br-', 'I-', 'Be2+',
     'Cu2+', 'Ni2+', 'Pt2+', 'Zn2+', 'Co2+', 'Pd2+', 'Ag2+', 'Cr2+', 'Fe2+',
     'Mg2+', 'V2+', 'Mn2+', 'Hg2+', 'Cd2+', 'Yb2+', 'Ca2+', 'Sn2+', 'Pb2+',
     'Eu2+', 'Sr2+', 'Sm2+', 'Ba2+', 'Ra2+', 'Al3+', 'Fe3+', 'Cr3+', 'In3+',
     'Tl3+', 'Y3+', 'La3+', 'Ce3+', 'Pr3+', 'Nd3+', 'Sm3+', 'Eu3+', 'Gd3+',
     'Tb3+', 'Dy3+', 'Er3+', 'Tm3+', 'Lu3+', 'Hf4+', 'Zr4+', 'Ce4+', 'U4+',
     'Pu4+', 'Th4+', 'X'
     ]
    for i in ions:
        try:
            c_atom_types.pop(i)
        except KeyError:
            continue
    ## Set up the dictionary to build the table from
    di_OD = { 'columns' : ('tinker', 'DH1', 'DH2', 'DH3', 'DH4', 'KPP')}
    ## use atom_types.values() for the numbers!
    origin = 1
    for i in range(len(dihedral1)):
    # for i in range(10):
        ## if only dihedral1 is X
        if dihedral1[i] == 999 and dihedral2[i] != 999 and dihedral3[i] != 999 and dihedral4[i] != 999:
            for a_key in c_atom_types.values():
                # di_OD.update({origin: ['torsion      ', dihedral1[i], dihedral2[i],\
                # dihedral3[i], dihedral4[i], di_line[i]] })
                di_OD.update({origin: ['torsion      ', a_key, dihedral2[i],\
                dihedral3[i], dihedral4[i], di_line[i]] })
                origin += 1
        ## if only dihedral2 is X
        elif dihedral2[i] == 999 and dihedral1[i] != 999 and dihedral3[i] != 999 and dihedral4[i] != 999:
            for a_key in c_atom_types.values():
                di_OD.update({origin: ['torsion      ', dihedral1[i], a_key,\
                dihedral3[i], dihedral4[i], di_line[i]] })
                origin += 1
        ## if only dihedral3 is X
        elif dihedral3[i] == 999 and dihedral1[i] != 999 and dihedral2[i] != 999 and dihedral4[i] != 999:
            for a_key in c_atom_types.values():
                di_OD.update({origin: ['torsion      ', dihedral1[i], dihedral2[i],\
                a_key, dihedral4[i], di_line[i]] })
                origin += 1
        ## if only dihedral4 is X
        elif dihedral4[i] == 999 and dihedral1[i] != 999 and dihedral2[i] != 999 and dihedral3[i] != 999:
            for a_key in c_atom_types.values():
                di_OD.update({origin: ['torsion      ', dihedral1[i], dihedral2[i],\
                dihedral3[i], a_key, di_line[i]] })
                origin += 1
        ## if dihedral1 and dihedral2 are X
        elif dihedral1[i] == 999 and dihedral2[i] == 999 and dihedral3[i] != 999 and dihedral4[i] != 999:
            for a_key in c_atom_types.values():
                for b_key in c_atom_types.values():
                    di_OD.update({origin: ['torsion      ', a_key, b_key,\
                    dihedral3[i], dihedral4[i], di_line[i]] })
                    origin += 1
                    di_OD.update({origin: ['torsion      ', b_key, a_key,\
                    dihedral3[i], dihedral4[i], di_line[i]] })
                    origin += 1
        ## if dihedral1 and dihedral3 are X
        elif dihedral1[i] == 999 and dihedral3[i] == 999 and dihedral2[i] != 999 and dihedral4[i] != 999:
            for a_key in c_atom_types.values():
                for b_key in c_atom_types.values():
                    di_OD.update({origin: ['torsion      ', a_key, dihedral2[i],\
                    b_key, dihedral4[i], di_line[i]] })
                    origin += 1
                    di_OD.update({origin: ['torsion      ', b_key, dihedral2[i],\
                    a_key, dihedral4[i], di_line[i]] })
                    origin += 1
        ## if dihedral1 and dihedral4 are X
        elif dihedral1[i] == 999 and dihedral4[i] == 999 and dihedral2[i] != 999 and dihedral3[i] != 999:
            for a_key in c_atom_types.values():
                for b_key in c_atom_types.values():
                    di_OD.update({origin: ['torsion      ', a_key, dihedral2[i],\
                    dihedral3[i], b_key, di_line[i]] })
                    origin += 1
                    di_OD.update({origin: ['torsion      ', b_key, dihedral2[i],\
                    dihedral3[i], a_key, di_line[i]] })
                    origin += 1
        ## if dihedral2 and dihedral3 are X
        elif dihedral2[i] == 999 and dihedral3[i] == 999 and dihedral1[i] != 999 and dihedral4[i] != 999:
            for a_key in c_atom_types.values():
                for b_key in c_atom_types.values():
                    di_OD.update({origin: ['torsion      ', dihedral1[i], a_key,\
                    b_key, dihedral4[i], di_line[i]] })
                    origin += 1
                    di_OD.update({origin: ['torsion      ', dihedral1[i], b_key,\
                    a_key, dihedral4[i], di_line[i]] })
                    origin += 1
        ## if dihedral2 and dihedral4 are X
        elif dihedral2[i] == 999 and dihedral4[i] == 999 and dihedral1[i] != 999 and dihedral2[i] != 999:
            for a_key in c_atom_types.values():
                for b_key in c_atom_types.values():
                    di_OD.update({origin: ['torsion      ', dihedral1[i], a_key,\
                    dihedral3[i], b_key, di_line[i]] })
                    origin += 1
                    di_OD.update({origin: ['torsion      ', dihedral1[i], b_key,\
                    dihedral3[i], a_key, di_line[i]] })
                    origin += 1
        ## if dihedral3 and dihedral4 are X
        elif dihedral3[i] == 999 and dihedral4[i] == 999 and dihedral1[i] != 999 and dihedral2[i] != 999:
            for a_key in c_atom_types.values():
                for b_key in c_atom_types.values():
                    di_OD.update({origin: ['torsion      ', dihedral1[i], dihedral2[i],\
                    a_key, b_key, di_line[i]] })
                    origin += 1
                    di_OD.update({origin: ['torsion      ', dihedral1[i], dihedral2[i],\
                    b_key, a_key, di_line[i]] })
                    origin += 1
        ## To date, AMBER doesn't have any dihedrals that have 3 X's (which
        ## seems like a good thing!)
        # if none are X
        else:
            di_OD.update({origin: ['torsion      ', dihedral1[i], dihedral2[i],\
            dihedral3[i], dihedral4[i], di_line[i]] })
            origin += 1
    #
    # Turn that ordered dictionary into a pd.DataFrame
    init_df = pd.DataFrame.from_dict(di_OD, "index", columns=['tinker', 'DH1', 'DH2',
     'DH3', 'DH4', 'KPP'])
    #
    # Remove the row used to set-up the dictionary
    di_df = init_df.drop(['columns'])
    #
    # Drop duplicate values, keeping one copy
    di_df.drop_duplicates(keep = 'first', inplace = True, ignore_index=True)
    #
    # https://stackoverflow.com/questions/49538497/how-to-apply-function-to-slice-of-columns-using-loc
    # https://stackoverflow.com/questions/40474799/remove-reverse-duplicates-from-dataframe
    # https://stackoverflow.com/questions/51603520/pandas-remove-duplicates-that-exist-in-any-order
    #
    ## Set-up drop reverse duplicates
    columns = ['DH1', 'DH2', 'DH3', 'DH4']
    ## Attach extra columns of the dihedral
    dihed_df = pd.concat([di_df, pd.DataFrame(np.sort(di_df[columns], axis=1))], axis=1)
    ## Remove the duplicates, based on what seems like ABCD DCBA/DCBA ABCD
    dihed_df = dihed_df.drop_duplicates(dihed_df.columns.difference(di_df.columns))[di_df.columns]
    #
    dihed_df = dihed_df.reset_index()
    #
    # Remove new index that appears because of tinker????
    dihed_df = dihed_df.drop('index', axis=1)
    #
    print(" Achievement unlocked: dihedrals dealt with.\n")
    return dihed_df

def get_improper_torsions(param_dat, struct_dat):
    """Create lists of the improper terms, force constant values (phi_k), and
    equilibrium angles (theteq).
    """
    if struct_dat is None:
        imptors = list(param_dat.improper_periodic_types.items())
    else:
        imptors = list(struct_dat.improper_periodic_types.items())
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
    print(" Achievement unlocked: improper torsions acquired.\n")
    return imptor1, imptor2, imptor3, imptor4, imp_diheds_phik,\
     imp_diheds_phase, imp_diheds_per

# def get_atomlist(param_dat, rmin14_dict, eps14_dict):
def get_atomlist(param_dat, struct_dat):
    """Create a dataframe of all the different residues and update element and
    mass.
    """
    if struct_dat is None:
        res_types = list(param_dat.residues.items())
    else:
        # res_types = list(struct_dat.residues.items())
        ## Create a dictionary of all the different residue types
        residue_dict = {}
        for residue in param_dat.residues:
            residue_dict[residue.name] = pmd.modeller.residue.ResidueTemplate(name=residue.name).from_residue(residue)
        #
        res_types = list(residue_dict.items())
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
    print(" Achievement unlocked: generated an atom list.\n")
    return res_df

def write_params_noX(param_dat, struct_dat, atom_types,\
 bond1, bond2, bond_k, bond_req, \
 angle1, angle2, angle3, ang_k, ang_theteq,\
 imptor1, imptor2, imptor3, imptor4, imp_diheds_phik, imp_diheds_phase,\
 imp_diheds_per, dihed_df, res_df,\
 rmin14_dict, eps14_dict, ff_name, atom_connect, param_file_name):
    """Generate the new TINKER parameter file.
    """
    with open(param_file_name, "w+") as tp_out:
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
        if struct_dat is None:
            tp_out.write("vdw-14-scale            {}\n".format(param_dat.default_scnb))
            tp_out.write("chg-14-scale            {}\n".format(param_dat.default_scee))
        else:
            tp_out.write("vdw-14-scale            {}\n".format(struct_dat.default_scnb))
            tp_out.write("chg-14-scale            {}\n".format(struct_dat.default_scee))
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
        if struct_dat is None:
            for name in param_dat.titles:
                tp_out.write("{}.\n".format(name))
        else:
            for name in struct_dat.titles:
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
            ## Test_ion is to try and figure out letters of Atom Type
            ## Test_AT is to try the guess atom connectivity from Atom Type
            test_ion = atom_types.get(res_df.loc[i,'type'])
            test_AT = atom_connect.get(res_df.loc[i,'type'])
            if test_ion is None:
                guess_type = "999"
                guess_connect = "0 !! GUESSED TYPE AND CONNECTION"
                tp_out.write("atom       {:4}  {:2}    {:<4}  \"{:<30} {:3}   {:>7.2f}    {}\n".format((i+1),\
                 guess_type, res_df.loc[i,'type'],\
                 (res_df.loc[i, 'resname'] + " " + res_df.loc[i, 'name']+'\"'),\
                 res_df.loc[i, 'atomic_number'], res_df.loc[i, 'mass'],\
                 guess_connect))
            ## Attempt key retrieval
            elif test_AT is None:
                if res_df.loc[i,'type'][0] in ('C', 'N'):
                    guess_connect = "3 !! GUESSED CONNECTION"
                    tp_out.write("atom       {:4}  {:2}    {:<4}  \"{:<30} {:3}   {:>7.2f}    {}\n".format((i+1),\
                     atom_types[res_df.loc[i,'type']], res_df.loc[i,'type'],\
                     (res_df.loc[i, 'resname'] + " " + res_df.loc[i, 'name']+'\"'),\
                     res_df.loc[i, 'atomic_number'], res_df.loc[i, 'mass'],\
                     guess_connect))
                elif res_df.loc[i,'type'][0] in ('O', 'S'):
                    guess_connect = "2 !! GUESSED CONNECTION"
                    tp_out.write("atom       {:4}  {:2}    {:<4}  \"{:<30} {:3}   {:>7.2f}    {}\n".format((i+1),\
                     atom_types[res_df.loc[i,'type']], res_df.loc[i,'type'],\
                     (res_df.loc[i, 'resname'] + " " + res_df.loc[i, 'name']+'\"'),\
                     res_df.loc[i, 'atomic_number'], res_df.loc[i, 'mass'],\
                     guess_connect))
                elif res_df.loc[i,'type'][0] == 'H':
                    guess_connect = "1 !! GUESSED CONNECTION"
                    tp_out.write("atom       {:4}  {:2}    {:<4}  \"{:<30} {:3}   {:>7.2f}    {}\n".format((i+1),\
                     atom_types[res_df.loc[i,'type']], res_df.loc[i,'type'],\
                     (res_df.loc[i, 'resname'] + " " + res_df.loc[i, 'name']+'\"'),\
                     res_df.loc[i, 'atomic_number'], res_df.loc[i, 'mass'],\
                     guess_connect))
                else:
                    guess_connect = "0 !! GUESSED CONNECTION"
                    tp_out.write("atom       {:4}  {:2}    {:<4}  \"{:<30} {:3}   {:>7.2f}    {}\n".format((i+1),\
                     atom_types[res_df.loc[i,'type']], res_df.loc[i,'type'],\
                     (res_df.loc[i, 'resname'] + " " + res_df.loc[i, 'name']+'\"'),\
                     res_df.loc[i, 'atomic_number'], res_df.loc[i, 'mass'],\
                     guess_connect))
            else:
                tp_out.write("atom       {:4}  {:2}    {:<4}  \"{:<30} {:3}   {:>7.2f}    {}\n".format((i+1),\
                 atom_types[res_df.loc[i,'type']], res_df.loc[i,'type'],\
                 (res_df.loc[i, 'resname'] + " " + res_df.loc[i, 'name']+'\"'),\
                 res_df.loc[i, 'atomic_number'], res_df.loc[i, 'mass'],\
                 atom_connect[res_df.loc[i,'type']]))
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
            tp_out.write("vdw          {:2}               {:>8.4f}     {:>12.7f}\n".format((i+1), value1[1], value2[1]))
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
            tp_out.write("imptors      {:3}  {:3}  {:3}  {:3}           {:6.3f}  {:5.1f}  {:1}\n".format(imptor1[i],\
             imptor2[i], imptor3[i], imptor4[i], imp_diheds_phik[i], imp_diheds_phase[i], imp_diheds_per[i]))
        tp_out.write("\n\n")
        tp_out.write("            ############################\n")
        tp_out.write("            ##                        ##\n")
        tp_out.write("            ##  Torsional Parameters  ##\n")
        tp_out.write("            ##                        ##\n")
        tp_out.write("            ############################\n")
        tp_out.write("\n\n")
        ## If you don't want Xs in output
        dihed_df.to_string(tp_out, header=False, index=False, index_names=False)
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
    print(" Always remember to check the new parameters with TINKER analyze.")
    print(" You may be missing important terms, especially from solvent.\n")

def write_params(param_dat, struct_dat, atom_types,\
 bond1, bond2, bond_k, bond_req,\
 angle1, angle2, angle3, ang_k, ang_theteq,\
 imptor1, imptor2, imptor3, imptor4, imp_diheds_phik, imp_diheds_phase,\
 imp_diheds_per, dihedral1, dihedral2, dihedral3, dihedral4, di_line, res_df,\
 rmin14_dict, eps14_dict, ff_name, atom_connect, param_file_name):
    """Generate the new TINKER parameter file.
    """
    with open(param_file_name, "w+") as tp_out:
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
        if struct_dat is None:
            tp_out.write("vdw-14-scale            {}\n".format(param_dat.default_scnb))
            tp_out.write("chg-14-scale            {}\n".format(param_dat.default_scee))
        else:
            tp_out.write("vdw-14-scale            {}\n".format(struct_dat.default_scnb))
            tp_out.write("chg-14-scale            {}\n".format(struct_dat.default_scee))
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
        if struct_dat is None:
            for name in param_dat.titles:
                tp_out.write("{}.\n".format(name))
        else:
            for name in struct_dat.titles:
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
            ## Test_ion is to try and figure out letters of Atom Type
            ## Test_AT is to try the guess atom connectivity from Atom Type
            test_ion = atom_types.get(res_df.loc[i,'type'])
            test_AT = atom_connect.get(res_df.loc[i,'type'])
            if test_ion is None:
                guess_type = "999"
                guess_connect = "0 !! GUESSED TYPE AND CONNECTION"
                tp_out.write("atom       {:4}  {:2}    {:<4}  \"{:<30} {:3}   {:>7.2f}    {}\n".format((i+1),\
                 guess_type, res_df.loc[i,'type'],\
                 (res_df.loc[i, 'resname'] + " " + res_df.loc[i, 'name']+'\"'),\
                 res_df.loc[i, 'atomic_number'], res_df.loc[i, 'mass'],\
                 guess_connect))
            ## Attempt key retrieval
            elif test_AT is None:
                if res_df.loc[i,'type'][0] in ('C', 'N'):
                    guess_connect = "3 !! GUESSED CONNECTION"
                    tp_out.write("atom       {:4}  {:2}    {:<4}  \"{:<30} {:3}   {:>7.2f}    {}\n".format((i+1),\
                     atom_types[res_df.loc[i,'type']], res_df.loc[i,'type'],\
                     (res_df.loc[i, 'resname'] + " " + res_df.loc[i, 'name']+'\"'),\
                     res_df.loc[i, 'atomic_number'], res_df.loc[i, 'mass'],\
                     guess_connect))
                elif res_df.loc[i,'type'][0] in ('O', 'S'):
                    guess_connect = "2 !! GUESSED CONNECTION"
                    tp_out.write("atom       {:4}  {:2}    {:<4}  \"{:<30} {:3}   {:>7.2f}    {}\n".format((i+1),\
                     atom_types[res_df.loc[i,'type']], res_df.loc[i,'type'],\
                     (res_df.loc[i, 'resname'] + " " + res_df.loc[i, 'name']+'\"'),\
                     res_df.loc[i, 'atomic_number'], res_df.loc[i, 'mass'],\
                     guess_connect))
                elif res_df.loc[i,'type'][0] == 'H':
                    guess_connect = "1 !! GUESSED CONNECTION"
                    tp_out.write("atom       {:4}  {:2}    {:<4}  \"{:<30} {:3}   {:>7.2f}    {}\n".format((i+1),\
                     atom_types[res_df.loc[i,'type']], res_df.loc[i,'type'],\
                     (res_df.loc[i, 'resname'] + " " + res_df.loc[i, 'name']+'\"'),\
                     res_df.loc[i, 'atomic_number'], res_df.loc[i, 'mass'],\
                     guess_connect))
                else:
                    guess_connect = "0 !! GUESSED CONNECTION"
                    tp_out.write("atom       {:4}  {:2}    {:<4}  \"{:<30} {:3}   {:>7.2f}    {}\n".format((i+1),\
                     atom_types[res_df.loc[i,'type']], res_df.loc[i,'type'],\
                     (res_df.loc[i, 'resname'] + " " + res_df.loc[i, 'name']+'\"'),\
                     res_df.loc[i, 'atomic_number'], res_df.loc[i, 'mass'],\
                     guess_connect))
            else:
                tp_out.write("atom       {:4}  {:2}    {:<4}  \"{:<30} {:3}   {:>7.2f}    {}\n".format((i+1),\
                 atom_types[res_df.loc[i,'type']], res_df.loc[i,'type'],\
                 (res_df.loc[i, 'resname'] + " " + res_df.loc[i, 'name']+'\"'),\
                 res_df.loc[i, 'atomic_number'], res_df.loc[i, 'mass'],\
                 atom_connect[res_df.loc[i,'type']]))
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
            tp_out.write("vdw          {:2}               {:>8.4f}     {:>12.7f}\n".format((i+1), value1[1], value2[1]))
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
            tp_out.write("imptors      {:3}  {:3}  {:3}  {:3}           {:6.3f}  {:5.1f}  {:1}\n".format(imptor1[i],\
             imptor2[i], imptor3[i], imptor4[i], imp_diheds_phik[i], imp_diheds_phase[i], imp_diheds_per[i]))
        tp_out.write("\n\n")
        tp_out.write("            ############################\n")
        tp_out.write("            ##                        ##\n")
        tp_out.write("            ##  Torsional Parameters  ##\n")
        tp_out.write("            ##                        ##\n")
        tp_out.write("            ############################\n")
        tp_out.write("\n\n")
        ## If you WANT Xs in output
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
    print(" Always remember to check the new parameters with TINKER analyze.")
    print(" You may be missing important terms, especially from solvent.\n")

###########
##  RUN  ##
###########
struct_dat = det_structure(param_dat)

atom_types = get_ATs(param_dat, struct_dat)

atom_connect = guess_connectivity(atom_types)

rmin14_dict, eps14_dict = get_VDW(param_dat, struct_dat)

bond1, bond2, bond_k, bond_req = get_bonds(param_dat, atom_types)

angle1, angle2, angle3, ang_k, ang_theteq = get_angles(param_dat, atom_types)

dihedral1, dihedral2, dihedral3, dihedral4, dval = get_dihedrals(param_dat, \
 struct_dat, atom_types)

di_line = clean_dihedrals(dval)

if leave_as_X == False:
    dihed_df = build_X_dihedrals(dihedral1, dihedral2, dihedral3, dihedral4, atom_types, di_line)

imptor1, imptor2, imptor3, imptor4, imp_diheds_phik, imp_diheds_phase,\
 imp_diheds_per = get_improper_torsions(param_dat, struct_dat)

# res_df = get_atomlist(param_dat, rmin14_dict, eps14_dict)
res_df = get_atomlist(param_dat, struct_dat)

if leave_as_X == True:
    write_params(param_dat, struct_dat, atom_types,\
     bond1, bond2, bond_k, bond_req, \
     angle1, angle2, angle3, ang_k, ang_theteq,\
     imptor1, imptor2, imptor3, imptor4, imp_diheds_phik, imp_diheds_phase,\
     imp_diheds_per, dihedral1, dihedral2, dihedral3, dihedral4, di_line, res_df,\
     rmin14_dict, eps14_dict, ff_name, atom_connect, param_file_name)
else:
    write_params_noX(param_dat, struct_dat, atom_types,\
     bond1, bond2, bond_k, bond_req, \
     angle1, angle2, angle3, ang_k, ang_theteq,\
     imptor1, imptor2, imptor3, imptor4, imp_diheds_phik, imp_diheds_phase,\
     imp_diheds_per, dihed_df, res_df,\
     rmin14_dict, eps14_dict, ff_name, atom_connect, param_file_name)
