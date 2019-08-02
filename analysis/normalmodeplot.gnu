set encoding iso_8859_1
set term postscript eps enhanced color font "Helvetica, 20";

unset ytics
set xlabel "Residue Number"
set ylabel "PCA Square Fluctuations"
set xrange [1:500]
set xtics border nomirror out
set x2tics border nomirror out rotate by 25 ("DNA" 400)
set boxwidth 0.25
set style fill solid

set output "WT_protein_system_modes.eps";
plot "WT-system/WT-system_mode1.dat" u 1:2 w boxes t "Mode 1" lw 3, \
"WT-system/WT-system_mode2.dat" u 1:2 w boxes t "Mode 2" lw 3, \
"WT-system/WT-system_mode3.dat" u 1:2 w boxes t "Mode 3" lw 3;

set output "MUT-A-system_modes.eps";
plot "MUT-A-system/MUT-A-system_mode1.dat" u 1:2 w boxes t "Mode 1" lw 3, \
"MUT-A-system/MUT-A-system_mode2.dat" u 1:2 w boxes t "Mode 2" lw 3, \
"MUT-A-system/MUT-A-system_mode3.dat" u 1:2 w boxes t "Mode 3" lw 3;

set output "MUT-B-system_modes.eps";
plot "MUT-B-system/MUT-B-system_mode1.dat" u 1:2 w boxes t "Mode 1" lw 3, \
"MUT-B-system/MUT-B-system_mode2.dat" u 1:2 w boxes t "Mode 2" lw 3, \
"MUT-B-system/MUT-B-system_mode3.dat" u 1:2 w boxes t "Mode 3" lw 3;

set output "MUT-C-system_modes.eps";
plot "MUT-C-system/MUT-C-system_mode1.dat" u 1:2 w boxes t "Mode 1" lw 3, \
"MUT-C-system/MUT-C-system_mode2.dat" u 1:2 w boxes t "Mode 2" lw 3, \
"MUT-C-system/MUT-C-system_mode3.dat" u 1:2 w boxes t "Mode 3" lw 3;
