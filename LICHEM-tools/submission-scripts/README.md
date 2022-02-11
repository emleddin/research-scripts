# submissions-scripts

Various queue scripts for running LICHEM.
Modulefile examples for LICHEM are provided
[elsewhere](https://github.com/emleddin/research-scripts/tree/main/modulefiles)
in this repository.

## Cruntch3

[Cruntch3](https://chemistry.unt.edu/~cruntch/) is UNT Chemistry's cluster.
It uses the PBS queueing system.

### `run-serial-cruntch3.sh`

A PBS script for running SP/DFP/opt LICHEM jobs.

### `run-QSM-cruntch3.sh`

A PBS script for running QSM LICHEM jobs.

## Stampede2

[Stampede2](https://portal.tacc.utexas.edu/user-guides/stampede2) is a
TACC resource. It uses the SLURM queueing system.
To use LICHEM with Gaussian, you will need to sign the TACC
[license agreement](https://portal.tacc.utexas.edu/software/gaussian).

Stampede2-specific modules are provided in the
[lmod](https://github.com/emleddin/research-scripts/tree/main/modulefiles/lmod/lichem)
section of this repository.

### `run-serial-stampede2.sh`

A SLURM script for running SP/DFP/opt LICHEM jobs.

### `run-QSM-stampede2.sh`

A SLURM script for running QSM LICHEM jobs.

### `run-stampede2.sh`

A SLURM script with comments for all LICHEM job types.
