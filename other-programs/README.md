# other-programs

## `.bash_aliases`
Aliases and commands you like having on a system.
To use this file directly, add these lines to your
`~/.bashrc` or `~/.bash_profile`.

```bash
if [ -f ~/.bash_aliases ]; then
    . ~/.bash_aliases
fi
```

## `.gnuplot`
A configuration file for `gnuplot` with a `colorsequence` based off of `podo`.
This should be saved as `~/.gnuplot`.

## `.vmdrc`
My configuration file for `VMD` with things such as using black for labels.
This should be saved as `~/.vmdrc`.

## `labeled-gif.sh`
A bash script for taking individual saved PNG or TGA files, labeling them with
the frame number in the top-right corner, and saving as a continuously looping
gif.
An example usage of this would be QM/MM reactions where you need to modify
visible bonds on a per-frame basis.
The script requires [ImageMagick](https://imagemagick.org/script/download.php).
Each file is assumed to be named similarly, for example:
```
$ ls
 my_vmd_image0.png   my_vmd_image1.png    my_vmd_image2.png    my_vmd_image3.png
 my_vmd_image4.png   my_vmd_image5.png    my_vmd_image6.png    my_vmd_image7.png
 my_vmd_image8.png   my_vmd_image9.png   my_vmd_image10.png   my_vmd_image11.png
my_vmd_image12.png
```

## `midasrc`
A configuration file for `UCSF Chimera` (which, in my experience, is applied
after opening Chimera's command line). This should be stored as
`~/.chimera/midasrc`.

## `run-caver.sh`
This is a script for running the command-line [CAVER](http://www.caver.cz/)
from multiple starting point residues in the same protein.

## `traj-optix.tcl`
A VMD script for saving transparent images of a trajectory (since MovieMaker
doesn't allow the render commands to be modified on all devices).
It requires
[VMD 1.9.4](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD)
and a GPU, since it uses the internal Tachyon OptiX rendering with a
transparent background.
Specify the step size (`incrementby`) and number of frames to save (`fullframe`;
starting with 1).
The script will suggest commands for converting files with multiple extensions
and making a gif (if intended; requires
[ImageMagick](https://imagemagick.org/script/download.php)]).

To run the script, have it saved in the directory you call `vmd` from and run
`source traj-optix.tcl` in the VMD console.

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
