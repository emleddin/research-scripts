# modulefiles

Drafts of modulefiles for different programs are included here.
You will have to modify the paths the specify where things are installed on the
system you are using!!!

## env-mod

This directory contains modules written in Tcl for the
[Environment Modules](http://modules.sourceforge.net/) program.

A good starting point for writing your own modulefiles is the 
[University of Sheffield HPC documentation](https://docs.hpc.shef.ac.uk/en/latest/sharc/software/apps/index.html).

### General Commands for a Modulefile

- `#%Module`: specifies it's a modulefile. You can use a number to denote a
minimum version of the *Environment Modules* program to use, like `#%Module1.0`
- `#`: the comment character
- `set`: used to set-up variables
- `proc`: creates a process
- `global`: sets previous variable as a global variable for use
- `puts stderr`: prints to standard error (terminal)
- `module-whatis`: super short description of the module
- `module load`: tells it to load another module
- `prepend_path`: adds to the beginning of the defined path.
This can apply to `LD_LIBRARY_PATH` and more!
- `setenv`: sets an environment variable (e.g., `AMBERHOME`)
- `set-alias`: sets up an alias (helpful for weirdly named executables)

## lmod

This directory contains modules written in Lua for the
[Lmod](https://lmod.readthedocs.io/en/latest/) program,
which most TACC resources use.

### General Commands for a Modulefile

- `--`: the comment character
- `set`: used to set-up variables
- `local`: sets variables to use in the modulefile
- `help`: prints a defined help message
- `whatis`: used for the description of the module
- `load("")`: tells it to load another module
- `prepend_path`: adds to the beginning of the defined path.
This can apply to `LD_LIBRARY_PATH` and more!
- `setenv`: sets an environment variable (e.g., `AMBERHOME`)

## Setting Up Modules

Modulefiles work by naming a top-level directory with all of your modules,
individual subdirectories named for each program, and finally files for each
version (or true module) of that program.

Environment Module modulefiles typically do not have an extension, but Lmod
files have the `.lua` extension.

Environment Modules Example:
```
$ tree moduleifles/ -a
modulefiles/
├── cmake
│   └── 3.18.2
├── gromacs
│   ├── 2018.8
│   ├── 2019.6-plumed-2.6.1
│   ├── 2020.3
│   └── 4.5.5
├── lichem
│   ├── 2018.07.23
│   ├── 2019.06.05
│   ├── 2020.12.10
│   └── .version
└── tinker
    ├── 7.1
    ├── 8.4
    ├── hp-1.1
    └── .version

4 directories, 13 files
```

Lmod Example:
```
$ tree modulefiles/
modulefiles/
├── amber
│   ├── 18-GEM.lua
│   └── .modulerc.lua
├── lichem
│   ├── 2020-tinker7.1.lua
│   ├── 2020-tinker8.7.lua
│   ├── 2020-tinkerHP-1.2.lua
│   ├── .modulerc.lua
│   └── MPI-tinker7.lua
└── tinker
    ├── 7.1.lua
    ├── 8.7.lua
    ├── HP-1.2.lua
    └── .modulerc.lua

3 directories, 11 files
```

## Recognizing your Modulefiles

You can allow both Environment Modules and Lmod to recognize your modules
by using `module use.`
This command allows you to specify a top-level directory for your
module files.
The top-level directory can have whatever name you choose, since it's directly
specified.
This command can easily be added to your `~/.bashrc` file to source your
modules on login.

```
module use /path/to/top-level/modulefiles/directory
```

You can also use the `use.own` module included in the Environment Modules
implementation to recognize modules.
Here, you must put your modules in `$HOME/privatemodules`, because that is the
path that is expected.
Whenever you want to use your modules, you will first use `module load use.own`
to recognize them.
Often, using this method will load those modules last, as opposed to your
personal modules taking priority.
This can cause issues if you are not specifying versions when
using `module load`.

## Versioning

Environment Modules use a `.version` file in the program-level directory,
and Lmod uses a `.modulerc.lua` file in the program-level directory.

**`.version`**

```tcl
#%Module

set Modulesversion "1.9.3"
```

**`.modulerc.lua`**

```lua
module_version("1.9.3", "default")
```

## Working with Conda Environments

If scripting a task that uses a Conda environment, you'll have to directly
source that Conda environment.
```
module load anaconda/3
source activate $HOME/.conda/envs/name-of-env
```
