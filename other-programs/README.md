# other-programs

## `.gnuplot`
A configuration file for `gnuplot` with a `colorsequence` based off of `podo`.
This should be saved as `~/.gnuplot`.

## `.vmdrc`
My configuration file for `VMD` with things such as using black for labels.
This should be saved as `~/.vmdrc`.

## `midasrc`
A configuration file for `UCSF Chimera` (which, in my experience, is applied
after opening Chimera's command line). This should be stored as
`~/.chimera/midasrc`.

## `run-caver.sh`
This is a script for running the command-line [CAVER](http://www.caver.cz/)
from multiple starting point residues in the same protein.

## `traj-pov.tcl`
A VMD script for saving transparent images of a trajectory (since MovieMaker
doesn't allow the render commands to be modified on all devices).
It requires [POV-Ray 3.x](http://www.povray.org/download/).
Specify the step size (`incrementby`) and number of frames to save (`fullframe`;
starting with 1).
The script will suggest commands for converting files with multiple extensions
and making a gif (if intended; requires
[ImageMagick](https://imagemagick.org/script/download.php)]).

To run the script, have it saved in the directory you call `vmd` from and run
`source traj-pov.tcl` in the VMD console.
