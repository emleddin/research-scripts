set encoding iso_8859_1
set term pngcairo enhanced color font "Arial,30" size 1500,1050;

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

set xlabel "{/Symbol c} Degrees ({\260})"
set ylabel "{/Symbol d} Degrees ({\260})"

set key top left font "Arial,24" width 1 height 1

set style fill transparent solid 1.0 noborder
set style circle radius 0.5

set xrange[0:360]
set yrange[0:360]

## Use 00 in HEX for first plotted and 33 in hex for second plotted
## These are orange and blue
set output "delta-chi-WT-MUTA-resX.png";
plot "/abs/path/to/WT-1/WT_protein_system_backbone-X-dihedral.dat" u ($8+180):($5+180) w circles t "WT" lt rgb "#00FA7D00", \
"/abs/path/to/WT-2/WT_protein_system_backbone-X-dihedral.dat" u ($8+180):($5+180) w circles notitle lt rgb "#00FA7D00", \
"/abs/path/to/WT-3/WT_protein_system_backbone-X-dihedral.dat" u ($8+180):($5+180) w circles notitle lt rgb "#00FA7D00", \
"/abs/path/to/MUT-A-1/MUT_A_system_backbone-X-dihedral.dat" u ($8+180):($5+180) w circles t "MUT-A" lt rgb "#330000FA" fs solid 0.15, \
"/abs/path/to/MUT-A-2/MUT_A_system_backbone-X-dihedral.dat" u ($8+180):($5+180) w circles notitle lt rgb "#330000FA" fs solid 0.15, \
"/abs/path/to/MUT-A-3/MUT_A_system_backbone-X-dihedral.dat" u ($8+180):($5+180) w circles notitle lt rgb "#330000FA" fs solid 0.15;
