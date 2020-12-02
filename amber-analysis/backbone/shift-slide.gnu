set encoding iso_8859_1
set term pngcairo enhanced color font "Arial,24" size 750,525;

# Colors
# lc rgb "#0000FA" (blue)
# lc rgb "#FA7D00" (orange)
# lc rgb "#00FA00" (green)
# lc rgb "#00FAFA" (cyan)
# lc rgb "#FA007D" (pink)
# lc rgb "#FA8072" (salmon)
# lc rgb "#778899" (slate gray)
# lc rgb "#C1FFC1" (sea-green)
# lc rgb "#9400D3" (dark-violet)

set xlabel "Shift Degrees ({\260})"
set ylabel "Slide Degrees ({\260})"

set key top left font "Arial,20" width -1 height 1

set style fill transparent solid 1.0 noborder
set style circle radius 0.02


## Use 00 in HEX for first plotted and CC or 33 in hex for second plotted
## These are orange and blue
set output "WT-MUTA-shift-slide-C3A4.png";
plot "/abs/path/to/WT-1/BPstep-433-434.dat" u ($4):($5) w circles t "WT C3:G3', A4:T4'" lt rgb "#000000FA", \
"/abs/path/to/WT-2/BPstep-433-434.dat" u ($4):($5) w circles notitle lt rgb "#000000FA", \
"/abs/path/to/WT-3/BPstep-433-434.dat" u ($4):($5) w circles notitle lt rgb "#000000FA", \
"/abs/path/to/MUT-A-1/BPstep-433-434.dat" u ($4):($5) w circles t "MUT-A C3:G3', A4:T4'" lt rgb "#CCFA7D00" fs solid 0.15, \
"/abs/path/to/MUT-A-2/BPstep-433-434.dat" u ($4):($5) w circles notitle lt rgb "#CCFA7D00" fs solid 0.15, \
"/abs/path/to/MUT-A-3/BPstep-433-434.dat" u ($4):($5) w circles notitle lt rgb "#CCFA7D00" fs solid 0.15;

set style fill transparent solid 1.0 noborder
set style circle radius 0.02

set output "WT-MUTA-shift-slide-C3A4-flip.png";
plot "/abs/path/to/MUT-A-1/BPstep-433-434.dat" u ($4):($5) w circles t "MUT-A C3:G3', A4:T4'" lt rgb "#00FA7D00", \
"/abs/path/to/MUT-A-2/BPstep-433-434.dat" u ($4):($5) w circles notitle lt rgb "#00FA7D00", \
"/abs/path/to/MUT-A-3/BPstep-433-434.dat" u ($4):($5) w circles notitle lt rgb "#00FA7D00", \
"/abs/path/to/WT-1/BPstep-433-434.dat" u ($4):($5) w circles t "WT C3:G3', A4:T4'" lt rgb "#330000FA" fs solid 0.15, \
"/abs/path/to/WT-2/BPstep-433-434.dat" u ($4):($5) w circles notitle lt rgb "#330000FA" fs solid 0.15, \
"/abs/path/to/WT-3/BPstep-433-434.dat" u ($4):($5) w circles notitle lt rgb "#330000FA" fs solid 0.15;
