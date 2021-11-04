import sys
import os

# alloc: Queue allocation
alloc = "my_cpu_alloc"
# replicate: which replicate to specify (ex: 1, 2, 3)
replicate = "1"
# full_path: full path to the production files
full_path = "/absolute/path/to/production/files"
# traj_tag: how trajectories are named
# ex: use `WT_protein_system_wat` for `WT_protein_system_wat.prmtop`
traj_tag = "WT_protein_system_wat"
# file_ext: file extension on production MD (ex: mdcrd, nc)
file_ext = "nc"
# output_tag: tag for labeling output data files
output_tag = "WT_protein_system"
# start_range: first production MD file to read in for analysis
start_range = 26
# end_range: final production MD file to read in for analysis
end_range = 275
# mask: atom mask for cpptraj. Ex: `1-455`
mask = "1-455"
# fs: file separator in generated file names
fs = "_"
# work_dir: here to put the generated cpptraj files
work_dir = "/absolute/path/to/output/analysis"
# system: how to mark SNP/WT systems
system = "WT"
# Number of amino acids (for setting up things like EDA)
n_aa = "455"
# use_offset: have an offset for cpptraj_EDA? Options: True/False
use_offset = False
# offset: offset value (how many frames to keep. 10 = 1/10) Any integer.
# offset = "10"

#--------------------------------------------------

# Remove the trailing slash on a path, if present
full_path = full_path.rstrip('/')

sh_strip = ('cpptraj'+str(fs)+'strip.sh')
out_strip = ('cpptraj'+str(fs)+'strip.in')

sh_analy = ('cpptraj'+str(fs)+'analysis.sh')
out_analy = ('cpptraj'+str(fs)+'analysis.in')

sh_EDA = ('cpptraj'+str(fs)+'EDA.sh')
out_EDA = ('cpptraj'+str(fs)+'EDA.in')


def write_strip_bash(outfile, queue, rep, f_path, cpp_strip, tag, sys):
    """Creates the PBS script for submitting cpptraj strip jobs.

    Parameters
    ----------
    outfile : str
        The name of the output file.
    queue : str
        The name of the queue to submit the job through.
    rep : str
        The replicate the job is for.
    f_path: str
        The full path to the trajectory files, for finding the prmtop.
    cpp_strip : str
        The name of the cpptraj strip file that will be submitted.
    tag : str
        The tag used for the trajectory and input files.
    sys : str
        The system to make the analysis for.

    Returns
    -------
    outfile : txt
        The written file containing PBS script for submitting the cpptraj
        strip job.
    """
    f = open(outfile, "w+")
    f.write("#!/bin/bash\n")
    f.write(f"#PBS -q {queue}\n")
    # f.write("#PBS -l nodes=1:ppn=1,mem=20GB\n")
    f.write("#PBS -l nodes=1:ppn=20,mem=100GB\n")
    f.write("#PBS -j oe\n")
    f.write("#PBS -r n\n")
    f.write("#PBS -o err.error\n")
    f.write(f"#PBS -N S.{rep}.{sys}\n\n")
    f.write(f"prmfile={f_path}/{tag}.prmtop\n")
    f.write(f"cppfile={cpp_strip}\n\n")
    f.write("cd $PBS_O_WORKDIR\n")
    f.write("cat $PBS_NODEFILE  > $PWD/PBS_NODEFILE\n\n")
    # f.write("module load amber/19-serial\n\n")
    # f.write("$AMBERHOME/bin/cpptraj -p $prmfile -i $cppfile\n\n")
    f.write("module load amber/19-mvapich2\n\n")
    f.write("mpirun -np 20 -hostfile $PWD/PBS_NODEFILE \\\n")
    f.write(" $AMBERHOME/bin/cpptraj.MPI -p $prmfile -i $cppfile\n\n")
    f.close()


def write_strip_traj(outfile, file_sep, f_path, tag, f_ext, f_start, f_end,
                     cwd, sys, rep):
    """Creates the cpptraj script for writing an ASCII trajectory for EDA and a
    stripped NetCDF trajectory for analysis.

    Parameters
    ----------
    outfile : str
        The name of the output cpptraj strip file.
    file_sep: str
        Determines the separator to use for the file name.
    f_path : str
        The full path to the trajectory files (parm_path).
    tag : str
        The tag used for the trajectory and input files.
    f_ext : str
        Current file extension used for trajectory files (ex. mdcrd or nc).
    f_start : int
        The number associated with the first trajectory file to read in.
    f_end : int
        The number associated with the final trajectory file to read in.
    cwd : str
        The path to the current working directory for the overall analysis.
    sys : str
        The system to make the analysis for.
    rep : str
        The replicate of a given `sys`.

    Returns
    -------
    outfile : txt
        The input file for cpptraj stripping.
    """
    f = open(outfile, "w+")
    # start, end+1 for correct range
    for i in range(f_start, f_end + 1):
        f.write(f"trajin {f_path}/{tag}{file_sep}md{i}.{f_ext}\n")
    f.write("\n")
    #f.write(f"trajout {cwd}/EDA/{tag}{file_sep}{f_start}-{f_end}.mdcrd crd\n\n")
    #f.write(f"trajout {cwd}/analysis/EDA/{sys}/{rep}/{tag}{file_sep}{f_start}-{f_end}.mdcrd crd\n\n")
    #f.write("go\n\n")
    f.write("autoimage\n")
    f.write("strip :WAT,K+ outprefix strip nobox\n\n")
    f.write(f"trajout {tag}{file_sep}imaged{file_sep}{f_start}-{f_end}.nc cdf\n\n")
    f.close()


def write_analy_bash(outfile, queue, rep, cpp_analy, tag, sys):
    """Creates the PBS script for submitting cpptraj analysis jobs.

    Parameters
    ----------
    outfile : str
        The name of the output file.
    queue : str
        The name of the queue to submit the job through.
    rep : str
        The replicate the job is for.
    cpp_analy : str
        The name of the cpptraj analysis file that will be submitted.
    tag : str
        The tag used for the trajectory and input files.
    sys : str
        The system to make the analysis for.

    Returns
    -------
    outfile : txt
        The written file containing PBS script for submitting the cpptraj
        analysis job.
    """
    f = open(outfile, "w+")
    f.write("#!/bin/bash\n")
    f.write(f"#PBS -q {queue}\n")
    f.write("#PBS -l nodes=1:ppn=20,mem=20GB\n")
    f.write("#PBS -j oe\n")
    f.write("#PBS -r n\n")
    f.write("#PBS -o err.error\n")
    f.write(f"#PBS -N A.{rep}.{sys}\n\n")
    f.write(f"prmfile=strip.{tag}.prmtop\n")
    f.write(f"cppfile={cpp_analy}\n\n")
    f.write("cd $PBS_O_WORKDIR\n")
    f.write("cat $PBS_NODEFILE  > $PWD/PBS_NODEFILE\n\n")
    f.write("module load amber/19-mvapich2\n\n")
    f.write("mpirun -np 20 -hostfile $PWD/PBS_NODEFILE \\\n")
    f.write(" $AMBERHOME/bin/cpptraj.MPI -p $prmfile -i $cppfile\n\n")
    f.close()


def write_analy_traj(outfile, fs, f_path, sys_tag, f_start, f_end,
                     num_aa, res_mask, out_tag):
    """Creates the input file used for analysis with cpptraj.

    Parameters
    ----------
    outfile : str
        The name of the output cpptraj strip file.
    fs: str
        Determines the separator to use for the file name.
    f_path : str
        The full path to the trajectory files (parm_path).
    out_tag : str
        Specifies the inital title information for the output files.
    f_start : int
        The number associated with the first trajectory file to read in.
    f_end : int
        The number associated with the final trajectory file to read in.
    num_aa : int
        Number of amino acids in the system for secondary structure analysis.
    res_mask : str
        The residue mask for specific analyses.
    sys_tag : str
        The tag used for the trajectory and input files.

    Returns
    -------
    outfile : txt
        The input file for cpptraj analysis.
    """
    f = open(outfile, "w+")
    f.write("# Read in the crystal (pre-minimization) structure\n")
    f.write("# You need to specify a prmtop with it because you're reading in\n")
    f.write("# the stripped trajectory\n")
    f.write(f"parm {f_path}/{sys_tag}.prmtop [ref]\n")
    f.write(f"reference {f_path}/{sys_tag}.inpcrd parm [ref]\n")
    f.write("\n")
    f.write("# Read in the stripped trajectory\n")
    f.write(f"trajin {sys_tag}{fs}imaged{fs}{f_start}-{f_end}.nc\n")
    f.write("\n")
    f.write("autoimage\n\n")
    f.write(f"rms reference out test{fs}rms.dat :{res_mask} byres\n\n")
    f.write("# Get for correlation matrix (evecs = eigenvectors)\n")
    f.write(f"matrix out {out_tag}{fs}corr{fs}mat.dat name corr_mat byres :{res_mask} correl\n\n")
    f.write("# Get for normal modes\n")
    f.write(f"matrix out {out_tag}{fs}covar{fs}mat.dat name norm_mode :{res_mask}@CA,P,C4',C2 \\\n")
    f.write(" covar\n")
    f.write(f"diagmatrix norm_mode out {out_tag}{fs}evecs.out vecs 100 reduce \\\n")
    f.write(f" nmwiz nmwizvecs 100 nmwizfile {out_tag}{fs}100.nmd \\\n")
    f.write(f" nmwizmask :{res_mask}@CA,P,C4',C2\n\n")
    f.write(f"hbond out {out_tag}{fs}hbond.dat dist 3.0 \\\n")
    f.write(f" avgout {out_tag}{fs}hbond{fs}avg.dat\n\n")
    f.write(f"rms reference out {out_tag}{fs}total{fs}bb{fs}rms.dat \\\n")
    f.write(f" :{res_mask}@CA,P,O3',O5',C3',C4',C5'\n")
    f.write(f"rmsd :{res_mask} reference perres perresavg range {res_mask} \\\n")
    f.write(f" perresout {out_tag}{fs}rmsd{fs}byres.dat\n\n")
    f.write(f"atomicfluct :{res_mask} out {out_tag}{fs}rmsf{fs}byres.dat byres\n")
    f.write(f"secstruct :1-{num_aa} out {out_tag}{fs}secstruct.gnu\n\n")
#    # Add pucker
#    for i in range(292,301+1):
#        f.write(f"pucker pucker{i} :{i}@C1' :{i}@C2' :{i}@C3' :{i}@C4' :{i}@O4' \\\n")
#        f.write(f" out {out_tag}{fs}pucker{fs}{i}.dat type pucker range360\n")
#    f.write("\n")
    # End add
    f.close()

def write_EDA_bash(outfile, queue, rep, f_path, cpp_EDA, tag, sys):
    """Creates the PBS script for submitting cpptraj strip jobs.

    Parameters
    ----------
    outfile : str
        The name of the output file.
    queue : str
        The name of the queue to submit the job through.
    rep : str
        The replicate the job is for.
    f_path: str
        The full path to the trajectory files, for finding the prmtop.
    cpp_EDA : str
        The name of the cpptraj EDA file that will be submitted.
    tag : str
        The tag used for the trajectory and input files.
    sys : str
        The system to make the analysis for.

    Returns
    -------
    outfile : txt
        The written file containing PBS script for submitting the cpptraj
        strip job.
    """
    f = open(outfile, "w+")
    f.write("#!/bin/bash\n")
    f.write(f"#PBS -q {queue}\n")
    f.write("#PBS -l nodes=1:ppn=1,mem=20GB\n")
    #f.write("#PBS -l nodes=1:ppn=20,mem=100GB\n")
    f.write("#PBS -j oe\n")
    f.write("#PBS -r n\n")
    f.write("#PBS -o err.error\n")
    f.write(f"#PBS -N EDAcpp.{rep}.{sys}\n\n")
    f.write(f"prmfile={f_path}/{tag}.prmtop\n")
    f.write(f"cppfile={cpp_EDA}\n\n")
    f.write("cd $PBS_O_WORKDIR\n")
    f.write("cat $PBS_NODEFILE  > $PWD/PBS_NODEFILE\n\n")
    f.write("module load amber/19-serial\n\n")
    f.write("$AMBERHOME/bin/cpptraj -p $prmfile -i $cppfile\n\n")
    #f.write("module load amber/19-mvapich2\n\n")
    #f.write("mpirun -np 20 -hostfile $PWD/PBS_NODEFILE \\\n")
    #f.write(" $AMBERHOME/bin/cpptraj.MPI -p $prmfile -i $cppfile\n\n")
    f.close()


def write_EDA_traj(outfile, file_sep, f_path, tag, f_ext, f_start, f_end,
                     cwd, sys, rep, use_off, off):
    """Creates the cpptraj script for writing ASCII trajectories for EDA.
    Trajectories are written in groups of 10 (or smaller) from the original
    files.

    Parameters
    ----------
    outfile : str
        The name of the output cpptraj EDA file.
    file_sep: str
        Determines the separator to use for the file name.
    f_path : str
        The full path to the trajectory files (parm_path).
    tag : str
        The tag used for the trajectory and input files.
    f_ext : str
        Current file extension used for trajectory files (ex. mdcrd or nc).
    f_start : int
        The number associated with the first trajectory file to read in.
    f_end : int
        The number associated with the final trajectory file to read in.
    cwd : str
        The path to the current working directory for the overall analysis.
    sys : str
        The system to make the analysis for.
    rep : str
        The replicate of a given `sys`.

    Returns
    -------
    outfile : txt
        The input file for cpptraj EDA.
    """
    f = open(outfile, "w+")
    #
    # Get divisible by 10
    count = f_start
    while count%10!=0:
        count += 1
    first_pass = count - f_start
    # Print first group, if currently non-divisible by 10
    # Ex: if range is 5-100; this covers 5-9
    if count != f_start:
        for i in range(f_start, f_start + first_pass):
            f.write(f"trajin {f_path}/{tag}{file_sep}md{i}.{f_ext}\n")
        f.write("\n")
        if use_off == True:
            f.write(f"trajout {cwd}/EDA/{tag}{file_sep}{f_start}-{count-1}.mdcrd crd offset {off}\n")
        else:
            f.write(f"trajout {cwd}/EDA/{tag}{file_sep}{f_start}-{count-1}.mdcrd crd\n")
        f.write("go\n")
        f.write("clear trajin\n\n")
    # start, end+1 for correct range
    pvals = range(count, f_end+1)
    for i in range(0, len(pvals), 10):
        c = pvals[i:i + 10]
        c = list(c)
        for i in c:
            f.write(f"trajin {f_path}/{tag}{file_sep}md{i}.{f_ext}\n")
        f.write("\n")
        if use_off == True:
            f.write(f"trajout {cwd}/EDA/{tag}{file_sep}{c[0]}-{c[-1]}.mdcrd crd offset {off}\n")
        else:
            f.write(f"trajout {cwd}/EDA/{tag}{file_sep}{c[0]}-{c[-1]}.mdcrd crd\n")
        f.write("go\n")
        f.write("clear trajin\n\n")
    f.close()


# Write the stripped files
write_strip_bash(sh_strip, alloc, replicate, full_path, out_strip, traj_tag,
                 system)
write_strip_traj(out_strip, fs, full_path, traj_tag, file_ext, start_range,
                 end_range, work_dir, system, replicate)

# Write the analysis files
write_analy_bash(sh_analy, alloc, replicate, out_analy, traj_tag, system)
write_analy_traj(out_analy, fs, full_path, traj_tag, start_range, end_range,
                 n_aa, mask, output_tag)

# Make EDA directory if it doesn't exist
os.makedirs(work_dir+'/EDA', exist_ok=True)

if use_offset == True:
    if offset is None:
        raise OffsetError('No offset value given when requested. Exiting...')

# Write the EDA files
write_EDA_bash(sh_EDA, alloc, replicate, full_path, out_EDA, traj_tag, system)
write_EDA_traj(out_EDA, fs, full_path, traj_tag, file_ext, start_range,
                 end_range, work_dir, system, replicate, use_offset, offset)
