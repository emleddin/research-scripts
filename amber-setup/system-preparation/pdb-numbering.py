import MDAnalysis as mda
import pandas as pd

RCSB_PDB = "4NM6.pdb"
tleap_PDB = "TET2_WT_struct_from_LEAP.pdb"

PDBID = "4NM6"
save_list = "TET2_numbering.txt"

def load_PDB(in_pdb):
    """
    Load in the LICHEM XYZ.
    Parameters
    ----------
    in_pdb : str
        The path to/name of the input PDB file.
    Returns
    -------
    system : MDAnalysis.core.universe.Universe
        Information from in_pdb.
    """
    system = mda.Universe(in_pdb, format="PDB", dt=1.0, in_memory=True)
    #
    ## Remove the segment IDs (aka `SYSTEM`) for prettier AtomGroup printing
    for atom in system.atoms:
        atom.segment.segid = ''
    #
    return system

def set_vars(system):
    """
    Set up MDA variables for the system.
    Parameters
    ----------
    system : MDAnalysis.core.universe.Universe
        The MDA Universe containing input information.
    Returns
    -------
    sys_prot : MDAnalysis.core.groups.AtomGroup
        Atom Group of protein residues.
    sys_nuc : MDAnalysis.core.groups.AtomGroup
        Atom Group of nucleic residues.
    sys_prot_res : int
        Number of residues in sys_prot.
    sys_nuc_res : int
        Number of residues in sys_nuc.
    """
    ## Select Atom Groups
    sys_prot = system.select_atoms("protein")
    #sys_nuc = system.select_atoms("nucleic")
    ## Add non-standard nucleic acid to nucleic group!!
    sys_nuc = system.select_atoms("nucleic or resname 5CM or resname 5MC")
    ## Get counters
    sys_prot_res = len(sys_prot.residues)
    sys_nuc_res = len(sys_nuc.residues)
    return sys_prot, sys_prot_res, sys_nuc, sys_nuc_res

def add_res(res_list, res1, res2, b1_count, a1_count):
    """
    Add a residue to the residue list.
    Parameters
    ----------
    res_list : list
        A list containing the Resnames and ResIDs between the two lists.
    res1 : MDAnalysis.core.groups.Residue
        Specific residue in the RCSB structure.
    res2 : MDAnalysis.core.groups.Residue
        Specific residue in the LEAP structure.
    b1_count : int
        Counter for the RCSB residues.
    a1_count : int
        Counter for the LEAP residues.
    Returns
    -------
    Updated parameters.
    """
    #print(res1.resname, res1.resid, res2.resname, res2.resid)
    res_list.append((res1.resname, res1.resid, res2.resname, res2.resid))
    b1_count += 1
    a1_count += 1
    return res_list, res1, res2, b1_count, a1_count

## Cases for the if statements
def b1_res_missing(res2, a1_count):
    #print("MISSING", "NA", res2.resname, res2.resid)
    res_list.append(("MISSING", "NA", res2.resname, res2.resid))
    ## Continue counting a1 since it's present
    a1_count += 1
    return a1_count

def a1_res_missing(res1, b1_count):
    #print(res1.resname, res1.resid, "MISSING", "NA")
    res_list.append((res1.resname, res1.resid, "MISSING", "NA"))
    ## Continue counting b1 since it's present
    b1_count += 1
    return b1_count

def break_func(res1, res2, b1_count, a1_count):
    print("I broke somewhere along the way...")
    print("Maybe you used a non-wild-type structure???")
    print(res1.resname, res1.resid, res2.resname, res2.resid)
    print(f"b1_count: {b1_count} a1_count: {a1_count}")
    print("Exiting...")

#-----------------------------------------------------------------------------#
## Letter Code Dictionary
lc_d = {'ALA': 'A', 'ARG': 'R', 'ASH': 'D', 'ASP': 'D', 'ASN': 'N',
        'CYM': 'C', 'CYS': 'C', 'CYX': 'C', 'GLH': 'E', 'GLU': 'E',
        'GLN': 'Q', 'GLY': 'G', 'HID': 'H', 'HIE': 'H', 'HIP': 'H',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYN': 'K', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T',
        'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'DA':  'dA', 'DC':  'dC', 'DG':  'dG', 'DT':  'dT',
        'DA5': 'dA', 'DC5': 'dC', 'DG5': 'dG', 'DT5': 'dT',
        'DA3': 'dA', 'DC3': 'dC', 'DG3': 'dG', 'DT3': 'dT',
        'RA':  'rA', 'RC':  'rC', 'RG':  'rG', 'RU':  'rU',
        'RA5': 'rA', 'RC5': 'rC', 'RG5': 'rG', 'RU5': 'rU',
        'RA3': 'rA', 'RC3': 'rC', 'RG3': 'rG', 'RU3': 'rU',
        'A':   'rA', 'C':   'rC', 'G':   'rG', 'U':   'rU',
        'A5':  'rA', 'RC':  'rC', 'RG':  'rG', 'RU':  'rU',
        'A3':  'rA', 'RC':  'rC', 'RG':  'rG', 'RU':  'rU',
        'ACE': 'X', 'NME': 'X', 'MISSING': 'X'
        }

#-----------------------------------------------------------------------------#
## Load in structures
biol_struct = load_PDB(RCSB_PDB)
amber_struct = load_PDB(tleap_PDB)

## Set up variables for loop
b1_prot, b1_prot_res, b1_nuc, b1_nuc_res = set_vars(biol_struct)
a1_prot, a1_prot_res, a1_nuc, a1_nuc_res = set_vars(amber_struct)


## Set empty list for resid and resname information
res_list = []

print("Read in the structures correctly!\n")
print("Checking protein...")

## Figure out the protein portion
b1_count = 0
a1_count = 0
while b1_count < b1_prot_res or a1_count < a1_prot_res:
    if b1_count == b1_prot_res:
        res2 = a1_prot.residues[a1_count]
        a1_count = b1_res_missing(res2, a1_count)
        continue
    elif a1_count == a1_prot_res:
        res1 = b1_prot.residues[b1_count]
        b1_count = a1_res_missing(res1, b1_count)
        continue
    else:
        res1 = b1_prot.residues[b1_count]
        res2 = a1_prot.residues[a1_count]
    #
    if res1.resname == res2.resname:
        res_list, res1, res2, b1_count, a1_count = add_res(res_list, \
         res1, res2, b1_count, a1_count)
    ## RCSB vs LEAP naming
    elif res2.resname in ('HID', 'HIE', 'HIP'):
        if res1.resname == "HIS":
            res_list, res1, res2, b1_count, a1_count = add_res(res_list, \
             res1, res2, b1_count, a1_count)
        elif b1_prot.residues[b1_count].resid != (b1_prot.residues[b1_count].resid+1):
            a1_count = b1_res_missing(res2, a1_count)
        elif a1_prot.residues[a1_count].resid != (a1_prot.residues[a1_count].resid+1):
            b1_count = a1_res_missing(res1, b1_count)
        else:
            break_func(res1, res2, b1_count, a1_count)
            break
    elif res2.resname in ('ASH'):
        if res1.resname == "ASP":
            res_list, res1, res2, b1_count, a1_count = add_res(res_list, \
             res1, res2, b1_count, a1_count)
        elif b1_prot.residues[b1_count].resid != (b1_prot.residues[b1_count].resid+1):
            a1_count = b1_res_missing(res2, a1_count)
        elif a1_prot.residues[a1_count].resid != (a1_prot.residues[a1_count].resid+1):
            b1_count = a1_res_missing(res1, b1_count)
        else:
            break_func(res1, res2, b1_count, a1_count)
            break
    elif res2.resname in ('GLH'):
        if res1.resname == "GLU":
            res_list, res1, res2, b1_count, a1_count = add_res(res_list, \
             res1, res2, b1_count, a1_count)
        elif b1_prot.residues[b1_count].resid != (b1_prot.residues[b1_count].resid+1):
            a1_count = b1_res_missing(res2, a1_count)
        elif a1_prot.residues[a1_count].resid != (a1_prot.residues[a1_count].resid+1):
            b1_count = a1_res_missing(res1, b1_count)
        else:
            break_func(res1, res2, b1_count, a1_count)
            break
    elif res2.resname in ('CYM', 'CYX'):
        if res1.resname == "CYS":
            res_list, res1, res2, b1_count, a1_count = add_res(res_list, \
             res1, res2, b1_count, a1_count)
        elif b1_prot.residues[b1_count].resid != (b1_prot.residues[b1_count].resid+1):
            a1_count = b1_res_missing(res2, a1_count)
        elif a1_prot.residues[a1_count].resid != (a1_prot.residues[a1_count].resid+1):
            b1_count = a1_res_missing(res1, b1_count)
        else:
            break_func(res1, res2, b1_count, a1_count)
            break
    elif res2.resname in ('LYN'):
        if res1.resname == "LYS":
            res_list, res1, res2, b1_count, a1_count = add_res(res_list, \
             res1, res2, b1_count, a1_count)
        elif b1_prot.residues[b1_count].resid != (b1_prot.residues[b1_count].resid+1):
            a1_count = b1_res_missing(res2, a1_count)
        elif a1_prot.residues[a1_count].resid != (a1_prot.residues[a1_count].resid+1):
            b1_count = a1_res_missing(res1, b1_count)
        else:
            break_func(res1, res2, b1_count, a1_count)
            break
    elif res2.resname in ('ACE', 'NME'):
        a1_count = b1_res_missing(res2, a1_count)
    ## Deal with missing residues
    else:
        if b1_prot.residues[b1_count].resid != (b1_prot.residues[b1_count].resid+1):
            a1_count = b1_res_missing(res2, a1_count)
        elif a1_prot.residues[a1_count].resid != (a1_prot.residues[a1_count].resid+1):
            b1_count = a1_res_missing(res1, b1_count)
        else:
            break_func(res1, res2, b1_count, a1_count)
            break

print("FYI: I may have had a bit of trouble if you have a linker with")
print("repeating residues.\n")
print("Checking nucleics...")

## Figure out the nucleic acid portion
b1_count = 0
a1_count = 0
while b1_count < b1_nuc_res or a1_count < a1_nuc_res:
    if b1_count >= b1_nuc_res:
        res2 = a1_nuc.residues[a1_count]
        a1_count = b1_res_missing(res2, a1_count)
        continue
    elif a1_count >= a1_nuc_res:
        res1 = b1_nuc.residues[a1_count]
        b1_count = a1_res_missing(res1, b1_count)
        continue
    else:
        res1 = b1_nuc.residues[b1_count]
        res2 = a1_nuc.residues[a1_count]
    #
    if res1.resname == res2.resname:
        res_list, res1, res2, b1_count, a1_count = add_res(res_list, \
         res1, res2, b1_count, a1_count)
    ## RCSB vs LEAP naming -- Last base in strand
    elif res2.resname in ('DC3', 'DG3', 'DA3', 'DT3', 'RC3', 'RG5',
     'RA3', 'RU3', 'C3', 'G3', 'A3', 'U3'):
        if res1.resname == res2.resname[:-1]:
            res_list, res1, res2, b1_count, a1_count = add_res(res_list, \
             res1, res2, b1_count, a1_count)
        ## If the b1 strand moves to second while a1 continues on first
        elif b1_nuc.residues[b1_count].resid > (b1_nuc.residues[b1_count].resid+1):
            a1_count = b1_res_missing(res2, a1_count)
        ## LEAP numbers nucleics sequentially, so this shouldn't happen
        # elif a1_nuc.residues[a1_count].resid > (a1_nuc.residues[a1_count].resid+1):
        #     b1_count = a1_res_missing(res1, b1_count)
        else:
            break_func(res1, res2, b1_count, a1_count)
            break
    ## RCSB vs LEAP naming -- First base in strand
    elif res2.resname in ('DC5', 'DG5', 'DA5', 'DT5', 'RC5', 'RG5', 'RA5',
     'RU5', 'C5', 'G5', 'A5', 'U5'):
        if res1.resname == res2.resname[:-1]:
            res_list, res1, res2, b1_count, a1_count = add_res(res_list, \
             res1, res2, b1_count, a1_count)
        ## If the a1 strand moves to second while b1 continues on first
        elif b1_nuc.residues[b1_count].resid < (b1_nuc.residues[b1_count].resid+1):
            b1_count = a1_res_missing(res1, b1_count)
        ## LEAP numbers nucleics sequentially, so this shouldn't happen
        # elif a1_nuc.residues[a1_count].resid < (a1_nuc.residues[a1_count].resid+1):
        #     a1_count = b1_res_missing(res2, a1_count)
        else:
            break_func(res1, res2, b1_count, a1_count)
            break
    ## Deal with missing residues
    else:
        if b1_nuc.residues[b1_count].resid != (b1_nuc.residues[b1_count].resid+1):
            a1_count = b1_res_missing(res1, a1_count)
        elif a1_nuc.residues[a1_count].resid != (a1_nuc.residues[a1_count].resid+1):
            b1_count = a1_res_missing(res1, b1_count)
        else:
            break_func(res1, res2, b1_count, a1_count)
            break

print("Huzzah! Writing out key!")

## Create a df with res_list
res_list_df = pd.DataFrame(res_list, columns =['RCSB_ResName', 'RCSB_ResID',
 'LEAP_ResName', 'LEAP_ResID'])

## Attach 1 Letter Codes
res_list_df['LC1'] = res_list_df["LEAP_ResName"].map(lc_d)
## Replace `nan` in LC1 column with `X`
na_values = {"LC1": 'X'}
res_list_df = res_list_df.fillna(value=na_values)

with open(save_list, "w+") as f:
    f.write(f"Mapping RCSB/Uniprot Numbering to Modified {PDBID}\n\n")
    f.write("{:<12} {:<12} {:<12} {:<4} {:<12}\n".format(\
    "LEAP_ResName", "LEAP_ResID", "RCSB_ResName", "1LC", "RCSB_ResID"))
    for r in res_list_df.itertuples(index=True, name='Pandas'):
        f.write("{:<12} {:<12} {:<12} {:<4} {:<12}\n".format(\
        r.LEAP_ResName, r.LEAP_ResID, r.RCSB_ResName, r.LC1, r.RCSB_ResID))
    f.close()

## This won't print headers/columns correctly, but it's not "bad"
# ## Write title info
# with open(save_list, "w+") as f:
#     f.write("Mapping RCSB/Uniprot numbering to Model\n\n")
# f.close()
# ## Append res_list_df to existing file with title
# res_list_df.to_csv(save_list, sep='\t', index=False, encoding='utf8',
# header=True, mode='a')

#-----------------------------------------------------------------------------#
# ## BIOL  DF
# res1_df =  pd.DataFrame(columns = ['ResName', 'ResNum'])
# for res in biol_struct.residues:
#     res1_df = res1_df.append({
#     "ResName" : res.resname,
#     "ResNum" : res.resnum
#     }, ignore_index=True)
#
# ## AMBER DF
# res2_df =  pd.DataFrame(columns = ['ResName', 'ResNum'])
# for res in amber_struct.residues:
#     res2_df = res2_df.append({
#     "ResName" : res.resname,
#     "ResNum" : res.resnum
#     }, ignore_index=True)
