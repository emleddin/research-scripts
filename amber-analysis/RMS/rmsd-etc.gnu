set encoding iso_8859_1
set term postscript enhanced color font "Arial,24";

set xlabel "Time (ns)"
set ylabel "RMSD ({\305})"
set key left bottom Left reverse
#set yrange [1.5:4]

set output "rmsds.eps";
plot "WT-system/WT_protein_system_total_bb_rms.dat" u ($1/500):($2) w lines s bezier t "Wild Type" lw 4, \
"MUT-A-system/MUT_A_system_total_bb_rms.dat" u ($1/500):($2) w lines s bezier t "Mutation A" lw 4, \
"MUT-B-system/MUT_B_system_RNA_total_bb_rms.dat" u ($1/500):($2) w lines s bezier t "Mutation B" lw 4, \
"MUT-C-system/MUT_C_system_total_bb_rms.dat" u ($1/500):($2) w lines s bezier t "Mutation C" lw 4;

set xlabel "Time (ns)"
set ylabel "Number of hydrogen bonds"
set output "hbonds.eps";
plot "WT-system/WT_protein_system_hbond.dat" u ($1/500):($2) w lines s bezier t "Wild Type" lw 4, \
"MUT-A-system/MUT_A_system_hbond.dat" u ($1/500):($2) w lines s bezier t "Mutation A" lw 4, \
"MUT-B-system/MUT_B_system_hbond.dat" u ($1/500):($2) w lines s bezier t "Mutation B" lw 4, \
"MUT-C-system/MUT_C_system_hbond.dat" u ($1/500):($2) w lines s bezier t "Mutation C" lw 4;

set xlabel "Residue number"
set ylabel "RMSF ({\305})"
set key top left Left reverse
set xrange [0:455]
#set yrange [0:7]
set output "rmsf.eps";
plot "WT-system/WT_protein_system_rmsf_byres.dat" u 1:2 w lines t "Wild Type" lw 4, \
"MUT-A-system/MUT_A_system_rmsf_byres.dat" u 1:2 w lines t "Mutation A" lw 4, \
"MUT-B-system/MUT_B_system_rmsf_byres.dat" u 1:2 w lines t "Mutation B" lw 4, \
"MUT-C-system/MUT_C_system_rmsf_byres.dat" u 1:2 w lines t "Mutation C" lw 4;
