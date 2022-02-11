-- This is an LMOD module file for AMBER18 on Stampede2

local help_message=[[
This version of AMBER was built with GEM.

More information about AMBER can be found at:
https://ambermd.org/index.php
]]

-- help(help_msg)
help(help_message, "\n")

whatis("Name: AMBER")
whatis("Version: 2018")
whatis("Keywords: Biology, Chemistry, Molecular Dynamics, Application")
whatis("URL: https://ambermd.org/index.php")
whatis("Description: molecular dynamics simulation package")

-- Load things needed for compilation
load("fftw3/3.3.6", "netcdf/4.6.2", "intel/18.0.2", "impi/18.0.2")

-- Create environment variables.
-- FIXME!
-- Ex: These next two commands would set amber_dir to $HOME/bin/amber18
-- local home = os.getenv("HOME")
-- local lichem_root  = pathJoin(home, "bin/amber18/")
-- Ex: These two would set amber_dir to $WORK/programs/amber18
-- local work = os.getenv("WORK")
-- local lichem_root  = pathJoin(work, "programs/amber18/")
-- Ex: This version sets a full path to amber_dir
local amber_dir    = "/path/to/top-level/amber18"

-- Set up overall environment when module is loaded!
setenv("AMBERHOME",   amber_dir)
prepend_path("PATH",  pathJoin(amber_dir, "bin"))
prepend_path("LD_LIBRARY_PATH", pathJoin(amber_dir, "lib"))
append_path("PYTHONPATH", pathJoin(amber_dir, "lib/python/site-packages"))
