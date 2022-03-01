#!/usr/bin/env tclsh
## Transparent Trajectory Gif in VMD (Normal Modes? :D)
## EML 01 March 2022
##
## Modify the `SPECIFY USER VARIABLES` section below!
## Run with `source traj-optix.tcl` in VMD command line

## Tell the program to expect variables
variable basename
variable basefilename
variable mybasefilename
variable imgsavedir

### !!! SPECIFY USER VARIABLES !!!

## Base for output file name (i.e., use `test` to later create `test.gif`)
set basename "my_awesome_trajectory"

## Save directory -- either use a string or `[pwd]` for current working dir
#set imgsavedir "/full/path/to/save/dir"
set imgsavedir [pwd]
## Number of frames to increment by (integer!)
set incrementby "1"

## Total number of frames to play through (integer!)
set fullframe "50"

### !!! END USER VARIABLES !!!

#--------------------------------------#
#-- Advanced User Modification Below --#
#--------------------------------------#

## Set the file path using basename and imgsavedir
set mybasefilename [format "%s/$basename" $imgsavedir]

## Use a white background and turn off depthcueing for transparency to render!
color Display Background white
display depthcue off

## Make Tachyon use a transparent background (1.9.4+)
set env(VMDOPTIXWRITEALPHA) 1

## Go to initial frame
animate goto 0

## Play the trajectory forward
set myfullframe [expr {$fullframe + 1}]
for {set i 1} {$i < $myfullframe} {incr i $incrementby} {
    puts {}
    puts "Starting to process Frame $i"
    animate goto $i
    set basefilename [format "%s/$basename.%04d" $imgsavedir $i ]
    ## Set the render command
    render TachyonLOptiXInternal $basefilename.png
}

## Rock the trajectory back
set myfullframerev [expr {$fullframe + 1}]
for {set i 1} {$i < $myfullframerev} {incr i $incrementby} {
    puts {}
    ## Set the current frame in terminal output
    set currframe [expr {$fullframe + $i}]
    ## Move to the specific frame for animation
    set aniframe [expr {$fullframe - $i}]
    puts "Starting to process Frame $currframe"
    animate goto $aniframe
    set basefilename [format "%s/$basename.%04d" $imgsavedir $currframe ]
    ## Set the render command
    render TachyonLOptiXInternal $basefilename.png
}

## Print to VMD Terminal
puts {}
puts "Frames rendered!"
puts {}
puts "To actually generate the gif as you want it, try:"
puts "convert -delay 4 -layers optimize -dispose previous -loop 0 XXX.*.png XXX.gif"
puts "Tip: You may want to mess with the delay!"
puts {}
puts "Finished!"
