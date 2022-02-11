-- This is an LMOD module file for TINKER on Stampede2

local help_message=[[
Use this module through:
   module load tinker/7.1

For more information on Tinker, visit:
https://dasher.wustl.edu/tinker/

]]

help(help_message,"\n")

whatis("Name: Tinker")
whatis("Version: 7.1")
whatis("URL: https://dasher.wustl.edu/tinker/")

-- Create environment variables.
-- FIXME!
-- Ex: These next two commands would set tinker_root to $HOME/bin/tinker/bin
-- local home = os.getenv("HOME")
-- local tinker_root  = pathJoin(home, "bin/tinker/bin")
-- Ex: These two would set lichem_root to $WORK/programs/tinker/bin
-- local work = os.getenv("WORK")
-- local tinker_root  = pathJoin(work, "programs/tinker/bin/")
-- Ex: This version sets tinker_root directly, without the evironment variable
local tinker_root    = "/path/to/tinker/executables"

prepend_path( "PATH",   tinker_root)
