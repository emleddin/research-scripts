# Submission Scripts

## Summit
Scripts for running on [Summit](https://docs.olcf.ornl.gov/systems/summit_user_guide.html#).
Summit uses an LSF scheduler.

## UNT
Scripts for running on [Cruntch3](http://chemistry.unt.edu/~cruntch/) and
[Talon3](https://hpc.unt.edu/).
Cruntch3 uses a PBS (Torque) scheduler, and Talon3 uses a SLURM scheduler.

## XSEDE
Scripts for running on [Comet](https://portal.xsede.org/sdsc-comet).
Comet uses a SLURM scheduler.

## Scheduler Commands

### LSF

To submit a job:
```bash
$ bsub name-of-script.sh
```

To kill a job:
```bash
$ bkill JOBID
```

To check on jobs:
```bash
$ bjobs -u $USER
```

### PBS/Torque

To submit a job:
```bash
$ qsub name-of-script.sh
```

To kill a job:
```bash
$ qdel JOBID
```

To check on jobs:
```bash
$ qstat -u $USER
```

### SLURM

To submit a job:
```bash
$ sbatch name-of-script.sh
```

To kill a job:
```bash
$ scancel JOBID
```

To check on jobs:
```bash
$ squeue -u $USER
```
If that doesn't work (certain XSEDE resources), you can use `showq`.
