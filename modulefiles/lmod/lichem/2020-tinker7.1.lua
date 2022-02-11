-- This is an LMOD module file for LICHEM on Stampede2
-- It calls Gaussian 16 and Tinker 7

local help_message=[[
This version of LICHEM was built with:
x86_64
gcc/7.1.0
gnuparallel/git20180620

and uses Tinker 7.1

More information about LICHEM can be found at
https://github.com/CisnerosResearch/LICHEM
]]

-- help(help_msg)
help(help_message, "\n")

whatis("Name: LICHEM")
whatis("Version: 12-10-2020, commit f6b429a")
whatis("URL: https://github.com/CisnerosResearch/LICHEM")

-- Load original compilers
load("gcc", "gnuparallel")

-- Load Gaussian
-- Note: You MUST have a signed TACC agreement for this to work!
-- See https://portal.tacc.utexas.edu/software/gaussian
load("gaussian/16rA.03")

-- Load Tinker
load("tinker/7.1")

-- Create environment variables.
-- FIXME!
-- Ex: These next two commands would set lichem_root to $HOME/bin/LICHEM/bin
-- local home = os.getenv("HOME")
-- local lichem_root  = pathJoin(home, "bin/LICHEM/bin/")
-- Ex: These two would set lichem_root to $WORK/programs/LICHEM/bin
-- local work = os.getenv("WORK")
-- local lichem_root  = pathJoin(work, "programs/LICHEM/bin/")
-- Ex: This version sets lichem_root directly, without the evironment variable
local lichem_root    = "/path/to/lichem/executable"

-- Set library root
local lib_root     = "/usr/lib64"

-- Set up overall environment when module is loaded!
prepend_path( "PATH",   lichem_root)
prepend_path( "PATH",   lib_root)
setenv( "GAUSS_SCRDIR",  "/tmp")
