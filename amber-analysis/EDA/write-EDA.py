import sys

# python write-EDA.py alloc system replicate sys_tag n_res n_atom n_prot_at
# tot_residues nas_traj file_sep short_tag
# Ex: python write-EDA.py gac.cpu 5hmC r1 thing 1430 7848 12000 234973 blah_1-30.mdcrd -
# script_name = sys.argv[0]

# alloc: Queue allocation
alloc = "my_cpu_alloc"
# system: System specifier for job -N
system = "WT"
# replicate: which replicate to specify (ex: 1, 2, 3)
replicate = "1"
# traj_tag: how trajectories are named
# ex: use `WT_protein_system_wat` for `WT_protein_system_wat.prmtop`
traj_tag = "WT_protein_system_wat"
# n_res: Number of protein (aka non-solvent) residues
n_res = 455
# n_atom: Total number of atoms
n_atom = 59659
# n_prot_at: Number of protein (aka non-solvent) atoms
n_prot_at = "7275"
# tot_residues: Total number of residues
tot_residues = "17929"
# nas_traj: Name of the traj file (no autoimage, no stripping), before X-X.mdcrd
# Ex: nas_traj = "WT_protein_system_wat_imaged_" for "WT_protein_system_wat_imaged_1-100.mdcrd
nas_traj = "WT_protein_system_wat_imaged_"
# nas_start and nas_end: Start and end ranges, assuming that files
# are chunked in 10 or less
# Make sure if you want a nice, even number of ns
#  to set it up for 26-100 (75 total)
# Ex: if range is 26-100; the first file would cover 26-29
nas_start = 26
nas_end = 275
# fs: File separator
fs = "_"

## Set up the output file names
sh_file = ("EDA"+str(fs)+"script.sh")
inp_file = "EDA_new.inp"
ans_file = "ans.txt"

def write_eda_bash(outfile, queue, rep, sys, ans):
    """Creates the PBS script for submitting EDA jobs.

    Parameters
    ----------
    outfile : str
        The name of the output file.
    queue : str
        The name of the queue to submit the job through.
    rep : str
        The replicate the job is for.
    tag : str
        The tag used in the file names. This should be comprised of
        the project ID and system.
    ans : str
        The name of the file containing the answers to the program prompts.

    Returns
    -------
    outfile : txt
        The written file containing PBS script for submitting the EDA job.
    """
    f = open(outfile, "w+")
    f.write("#!/bin/bash\n")
    f.write(f"#PBS -q {queue}\n")
    f.write("#PBS -l nodes=1:ppn=10,mem=120GB\n")
    f.write("#PBS -j oe\n")
    f.write("#PBS -r n\n")
    f.write("#PBS -o EDA.error\n")
    f.write(f"#PBS -N EDA.{rep}.{sys}\n\n")
    f.write("# Load in the Intel compiler\n")
    f.write("module load intel/17.0\n\n")
    f.write("# Access the folder where the files are\n")
    f.write("cd $PBS_O_WORKDIR\n")
    f.write("cat $PBS_NODEFILE  > $PWD/PBS_NODEFILE\n\n")
    f.write("# Compile the EDA program\n")
    f.write("ifort Residue_E_Decomp_openmp.f90 -o Residue_E_Decomp_openmp.x -qopenmp\n\n")
    f.write("date\n\n")
    f.write("# Run the program; read in the prompt answers\n")
    f.write("# [Line 1: Name of input; Line 2: Name of prmtop]\n")
    f.write(f"./Residue_E_Decomp_openmp.x < {ans}\n\n")
    f.write("# Acquire the process ID for the program execution\n")
    f.write("proc_PID=$!\n\n")
    f.write("# Wait until the program execution is over before ending the script\n")
    f.write("wait $proc_PID\n\n")
    f.write("date\n\n")
    f.write("echo 'All done!'\n\n")
    f.close()


def write_eda_ans(outfile, inpfile, tag):
    """Write a file with the answers for running EDA non-interactively.

    Parameters
    ----------
    outfile : str
        The name of the output file.
    inpfile : str
        The name of the input file containing trajectory and atom-specific
        information.
    tag : str
        The tag used in the prmtop file name. This should be comprised of
        the project ID and system.

    Returns
    -------
    outfile : txt
        The file with the answers.
    """
    f = open(outfile, "w+")
    f.write(f"{inpfile}\n")
    f.write(f"{tag}.prmtop")


def write_eda_inp(outfile, n_residues, tot_atom, prot_at, tot_res, traj,
    traj_s, traj_e):
    """Creates the system-specific input file used for EDA.

    Parameters
    ----------
    outfile : str
        The name of the output file.
    n_residues : str
        The number of residues to calculate interactions for.
        This number will only work for 1-n_residues (ex: use 450 for residues
        1 to 450).
    tot_atom : str
        The total number of atoms in the prmtop/trajectory.
    prot_at: str
        The is the number of atoms for the quantity given in n_residues.
        This is the number of protein atoms to calculate interactions for based
        on the prmtop/trajectory.
    tot_res : str
        The total number of residues in the prmtop/trajectory.
    traj : str
        Starting string for the trajectory file names.
    traj_s : int
        The start number for the trajectory file range,.
    traj_e : int
        The end number for the trajectory file range.

    Returns
    -------
    outfile : txt
        The written file containing PBS script for submitting the cpptraj
        strip job.
    """
    f = open(outfile, "w+")
    f.write(f"{n_residues} !number of protein residues\n")
    # Figure out number of files before printing
    count = traj_s
    while count%10!=0:
        count += 1
    first_pass = count - traj_s
    # Determine if first_pass is used
    if first_pass != 0:
        nfiles = 1
    else:
        nfiles = 0
    # start, end+1 for correct range
    pvals = range(count, traj_e+1)
    for i in range(0, len(pvals), 10):
        nfiles += 1
    f.write(f"{nfiles} !number of files\n")
    f.write(f"{tot_atom} !total number of atoms\n")
    f.write(f"{prot_at} !number of protein atoms\n")
    f.write(f"{tot_res} !number of total residues\n")
    f.write("2000 !max number of types\n")
#    count = traj_s
#    while count%10!=0:
#        count += 1
#    first_pass = count - traj_s
    # Print first group, if currently non-divisible by 10
    # Ex: if range is 5-100; this covers 5-9
    if count != traj_s:
        f.write(f"{traj}{traj_s}-{count-1}.mdcrd\n")
    # start, end+1 for correct range
    pvals = range(count, traj_e+1)
    for i in range(0, len(pvals), 10):
        c = pvals[i:i + 10]
        c = list(c)
        f.write(f"{traj}{c[0]}-{c[-1]}.mdcrd\n")
    # f.write(f"{traj}")


write_eda_bash(sh_file, alloc, replicate, system, ans_file)
write_eda_ans(ans_file, inp_file, traj_tag)
# write_eda_inp(inp_file, n_res, n_atom, n_prot_at, tot_residues, nas_traj)
write_eda_inp(inp_file, n_res, n_atom, n_prot_at, tot_residues, nas_traj,
    nas_start, nas_end)
