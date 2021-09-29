#!/usr/bin/env tclsh
## Transparent 360 Trajectory Gif in VMD
## EML 29 Sept 2021
##
## Modify the `SPECIFY USER VARIABLES` section below!
## Run with `source opt-pov.tcl` in VMD command line

## Tell the program to expect variables
variable basename
variable basefilename
variable mybasefilename
variable imgsavedir

### !!! SPECIFY USER VARIABLES !!!

## Base for output file name (i.e., use `test` to later create `test.gif`)
set basename "my_awesome_trajectory"

## Save directory -- either use a string or `[pwd]` for current working dir
## Note: if not using the `[pwd]` or `/tmp`, you'll need a `povray.conf` file
#set imgsavedir "/full/path/to/save/dir"
set imgsavedir [pwd]
## Number of frames to increment by (integer!)
set incrementby "1"

## Total number of frames to play through (integer!)
set fullframe "13"

### !!! END USER VARIABLES !!!

#--------------------------------------#
#-- Advanced User Modification Below --#
#--------------------------------------#

## Set the file path using basename and imgsavedir
set mybasefilename [format "%s/$basename" $imgsavedir]

## Use a white background and turn off depthcueing for transparency to render!
color Display Background white
display depthcue off

## Go to initial frame
animate goto 0

set myfullframe [expr {$fullframe + 1}]
for {set i 1} {$i < $myfullframe} {incr i $incrementby} {
    puts {}
    puts "Starting to process Frame $i"
    animate goto $i
    set basefilename [format "%s/$basename.%04d.pov" $imgsavedir $i ]
    ## Set the render command
    render POV3 $basefilename povray +W%w +H%h -I%s -O%s.png -D +X +A +FN +UA
}

## Other render options
## --------------------
##    +W: width `+W%w` uses the current window width
##    +H: height. `+H%h` uses the current window height
##    +I: the name of the POV file to use. You may need the `.pov`
##        extension in $basefilename.
##    +O: the name of the final image. If you specify the `.pov` above, it'll
##        be something like `.pov.png`.
##    -D: switches off the graphics display while being rendered. Using `+D`
##        instead means the display will be on!
##    +X:  make rendering killable
##    +FT: write targa file (with `+O%s.tga`)
##    +FN: write PNG file (with `+O%s.png`)
##    +FJ: write JPG file (with `+O%s.jpg`)
##    +FB: write Bitmap file (with `+O%s.bmp`)
##    +FP: write PPM format (with `+O%s.ppm`)
##    +UA: use background transparency. To use this, save as filetype that 
##         supports transparency, like PNG.
##         The VMD background must be white for this option!!!
## Note: If you need super-high-res images, use `+W%w0 +H%h0` (it WILL be slow
##       but very beautiful)
## Note: XXXX is the $basename
##
##    TGA -- not transparent
##    render POV3 XXXX povray +W%w +H%h -I%s -O%s.tga -D +X +A +FT
##
## JPG -- not transparent
##   render POV3 XXXX povray +W%w +H%h -I%s -O%s.jpg -D +X +A +FJ
##
## PNG Format -- yes transparent
##    render POV3 XXXX povray +W%w +H%h -I%s -O%s.png -D +X +A +FN +UA
## --------------------

## Print to VMD Terminal
puts {}
puts "Frames rendered!"
puts {}
puts "You may want to amend the file extensions from .pov.png to .png:"
puts "for filename in *pov.png; do mv \$filename \${filename//.pov.png/.png}; done"
puts {}
puts "To actually generate the gif as you want it, try:"
puts "convert -delay 50 -layers optimize -dispose previous -loop 0 XXX.*.png XXX.gif"
puts "Tip: You may want to mess with the delay!"
puts {}
puts "Finished!"

