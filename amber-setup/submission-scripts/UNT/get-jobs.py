#!/usr/bin/env python3
import subprocess
import textwrap
import argparse

## Set up argument parsing
parser = argparse.ArgumentParser(
    description="Show jobs in queue. By default it prints the full queue."+
    "User-specific information can be requested.")
parser.add_argument('-u', '--user', type=str,
    nargs='?', const=subprocess.getoutput("echo $USER"),
    help="runs at user-level (uses `$USER` if none specified)")

## Parse any given arguments
args = parser.parse_args()

## Set for users (True) or site-wide (False)
if args.user:
    user_level = True
    user = args.user
else:
    user_level = False

# Get a list of all the nodes for cleaning
all_nodes = subprocess.getoutput("pbsnodes -l all")

## Set up flag for grep commands
if user_level == True:
    ## Set to given user
    if user != None:
        extra_flag=f"-u {user}"
else:
    extra_flag=""

## job_id = subprocess.getoutput("qstat -f -u ${USER} | grep 'Job Id:'")
## Escape f-string with double braces if directly including system variables
## Ex: job_id = subprocess.getoutput(f"qstat -f ${{USER}} | grep {grep_job_id}")
grep_job_id = "Job\ Id:"
job_id = subprocess.getoutput(f"qstat -f {extra_flag} | grep {grep_job_id}")

grep_job_owner = "Job_Owner"
job_owner = subprocess.getoutput(f"qstat -f {extra_flag} | grep {grep_job_owner}")

grep_queue = "queue"
queue = subprocess.getoutput(f"qstat -f {extra_flag} | grep {grep_queue}")

grep_job_name = "Job_Name"
job_name = subprocess.getoutput(f"qstat -f {extra_flag} | grep {grep_job_name}")

grep_job_state = "job_state"
job_state = subprocess.getoutput(f"qstat -f {extra_flag} | grep {grep_job_state}")

grep_wall_time = "resources_used.walltime"
wall_time = subprocess.getoutput(f"qstat -f {extra_flag} | grep {grep_wall_time}")

grep_res_list_nodes = "Resource_List.nodes"
res_list_nodes = subprocess.getoutput(f"qstat -f {extra_flag} | grep {grep_res_list_nodes}")

## Check for next match (exec_port) to include mutli-line output
grep_exec_host = "exec_host"
exec_host = subprocess.getoutput(f"qstat -f {extra_flag} | awk '/exec_host/{{flag=1}} /exec_port/{{flag=0}} flag'")

job_dir = subprocess.getoutput(f"qstat -f {extra_flag} | awk '/Variable_List/{{flag=1}} /etime/{{flag=0}} flag' | grep 'PBS_O_WORKDIR'")

def grep_job_id_clean(cmd_output):
    """
    Clean up the `Job ID:` output.
    """
    ## Replace instances of grep locator with nothing
    clean_output = cmd_output.replace("Job Id:", '')
    clean_output = clean_output.replace('.cruntch3.chem.unt.edu', '')
    ## Break into list at newline characters
    clean_output = clean_output.splitlines()
    ## Remove whitespace at ends
    clean_output = [j.strip() for j in clean_output]
    return clean_output

def grep_job_owner_clean(cmd_output, grep_locator):
    """
    Clean up the `Job_Owner` output.
    """
    ## Replace instances of grep locator with nothing
    clean_output = cmd_output.replace(grep_locator+" =", '')
    clean_output = clean_output.replace('@cruntch3.chem.unt.edu', '')
    clean_output = clean_output.replace('@login-15-21.local', '')
    ## Break into list at newline characters
    clean_output = clean_output.splitlines()
    ## Remove whitespace at ends
    clean_output = [j.strip() for j in clean_output]
    return clean_output

def grep_nodes_clean(cmd_output, grep_locator):
    """
    Clean up the `Resource_List.nodes` output.
    """
    ## Replace instances of grep locator with nothing
    clean_output = cmd_output.replace(grep_locator+" =", '')
    clean_output = clean_output.replace("ppn=","")
    ## Remove GPU info (details)
    for x in all_nodes:
        clean_output = clean_output.replace(x,"")
    ## Break into list at newline characters
    clean_output = clean_output.splitlines()
    ## Remove whitespace at ends
    clean_output = [j.strip() for j in clean_output]
    return clean_output

def get_my_nodes(cmd_output):
    """
    Clean up the `pbsnodes` output.
    """
    ## Remove instances of second-line output
    clean_nodes = cmd_output.replace("job-exclusive", "")
    clean_nodes = clean_nodes.replace("free", "")
    clean_nodes = clean_nodes.replace("offline", "")
    clean_nodes = clean_nodes.splitlines()
    clean_nodes = [j.strip() for j in clean_nodes]
    return clean_nodes

def grep_exec_clean(cmd_output, grep_locator):
    """
    Clean up the `exec_host` output.
    """
    ## Break into list at command (since it spans new lines)
    ## This removes the search string
    clean_output = cmd_output.strip().split(grep_locator+" =")
    ## Remove empty strings
    clean_output = list(filter(None, clean_output))
    ## Replace interspersed newlines with nothing
    clean_output = [j.replace('\n', '') for j in clean_output]
    ## Replace any whitespace in long strings
    clean_output = [j.replace('\t', '') for j in clean_output]
    ## Replaced + signs to be readable
    clean_output = [j.replace('+', ', ') for j in clean_output]
    ## Remove whitespace at ends
    clean_output = [j.strip() for j in clean_output]
    return clean_output

def grep_job_dir_clean(cmd_output):
    """
    Clean up the `job_dir` output.
    """
    ## Replace instances of grep locator with nothing
    clean_output = cmd_output.replace("PBS_O_", '')
    ## Break into list at newline characters
    clean_output = clean_output.splitlines()
    ## Remove whitespace at ends
    clean_output = [j.strip() for j in clean_output]
    clean_output = [j.rstrip(',') for j in clean_output]
    return clean_output

def grep_clean(cmd_output, grep_locator):
    """
    Clean up the general output information. It assumes a section is formatted
    like `thing = important_values`.
    """
    ## Replace instances of grep locator with nothing
    clean_output = cmd_output.replace(grep_locator+" =", '')
    ## Break into list at newline characters
    clean_output = clean_output.splitlines()
    ## Remove whitespace at ends
    clean_output = [j.strip() for j in clean_output]
    return clean_output

## Clean up the output
all_nodes = get_my_nodes(all_nodes)
job_id_clean = grep_job_id_clean(job_id)
job_owner_clean = grep_job_owner_clean(job_owner, grep_job_owner)
queue_clean = grep_clean(queue, grep_queue)
job_name_clean = grep_clean(job_name, grep_job_name)
job_state_clean = grep_clean(job_state, grep_job_state)
wall_time_clean = grep_clean(wall_time, grep_wall_time)
res_list_nodes_clean = grep_nodes_clean(res_list_nodes, grep_res_list_nodes)
exec_host_clean = grep_exec_clean(exec_host, grep_exec_host)
job_dir_clean = grep_job_dir_clean(job_dir)

## Set Up Header (in center)
if user_level == True:
    print("{0:^80}".format(f"Jobs for {user} on Cruntch3\n"))
else:
    print("{0:^80}".format(f"Jobs for all users on Cruntch3\n"))
    print("{0:^80}".format(f"{len(job_id_clean)} Total Jobs\n"))

## For some reason, it won't work unless first is forced to be a
## string with `!s`
print('{!s:<7} {:<8} {:<8} {:<18} {:<5} {:<10} {:<5} {:<11}'.format(
     'JobID', 'User', 'Queue', 'JobName', 'Stat', 'WallTime', 'Node',
     'NodeDet'))

print("==================================================================="+
      "=============")

## Print job info
for jidc, joc, qc, jnc, jsc, wc, rlnc, ehc, jdc in zip(
    job_id_clean, job_owner_clean, queue_clean,
    job_name_clean,
    job_state_clean, wall_time_clean, res_list_nodes_clean,
    exec_host_clean, job_dir_clean):
    print("{:<7} {:<8} {:<8} {:<18} {:<5} {:<10} {:<5} {:<11}".format(
    jidc, joc, qc, jnc, jsc, wc, rlnc, ehc[:11]
    ))

    ## If CPU nodes span multiple lines, print them
    if len(ehc) > 11:
        ehc_out = textwrap.TextWrapper(width=78,
           break_long_words=True).wrap(ehc[11:])
        for i in ehc_out:
            print("  {:<78}".format(i))

    ## Print the Job Directory
    print("  {:<78}".format(jdc))

## Print explanations
print("\nStatuses: R=Running, Q=Queued, H=Held, E=Exiting")
print("Nodes are presented as number of nodes : processors per node.")
