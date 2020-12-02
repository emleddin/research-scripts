set encoding iso_8859_1
set term postscript enhanced color font "Arial,24";

set xlabel "Time (ns)"
set ylabel "RMSD ({\305})"
set key right bottom Left reverse maxrows 2
#set yrange [1.5:4]

# Other Berendsen inputs are $1/500, these Langevin are $1/100
set output "rmsds-avg.eps";
plot "WT_protein_system_total_bb_rms_avgd.dat" u ($1/100):($2) w lines s bezier t "Wild Type" ls 1 lw 4, \
"WT_protein_system_total_bb_rms_avgd.dat" u ($1/100):($2+$3):($2-$3) w filledcurve fs solid 0.15 t "Std. Dev." ls 1, \
"MUT_A_system_total_bb_rms_avgd.dat" u ($1/100):($2) w lines s bezier t "Mutation A" ls 2 lw 4, \
"MUT_A_system_total_bb_rms_avgd.dat" u ($1/100):($2+$3):($2-$3) w filledcurve fs solid 0.15 t "Std. Dev." ls 2, \
"WT_protein_system_total_bb_rms_avgd.dat" u ($1/100):($2) w lines s bezier notitle ls 1 lw 4, \
"MUT_A_system_total_bb_rms_avgd.dat" u ($1/100):($2) w lines s bezier notitle ls 2 lw 4;

set xlabel "Time (ns)"
set ylabel "Number of hydrogen bonds"
#set key left bottom Left reverse

# Other Berendsen inputs are $1/500, these Langevin are $1/100
set output "hbonds-avg.eps";
plot "WT_protein_system_hbond_3trial_avgd.dat" u ($1/100):($2) w lines s bezier t "Wild Type" ls 1 lw 4, \
"WT_protein_system_hbond_3trial_avgd.dat" u ($1/100):($2+$3):($2-$3) w filledcurve fs solid 0.15 t "Std. Dev." ls 1, \
"MUT_A_system_hbond_3trial_avgd.dat" u ($1/100):($2) w lines s bezier t "Mutation A" ls 2 lw 4, \
"WT_protein_system_hbond_3trial_avgd.dat" u ($1/100):($2+$3):($2-$3) w filledcurve fs solid 0.15 t "Std. Dev." ls 2, \
"MUT_A_system_hbond_3trial_avgd.dat" u ($1/100):($2) w lines s bezier notitle ls 1 lw 4, \
"MUT_A_system_hbond_3trial_avgd.dat" u ($1/100):($2) w lines s bezier notitle ls 2 lw 4;

set xlabel "Residue number"
set ylabel "RMSF ({\305})"
set key left top Left reverse maxrows 2
#set yrange [0:7]
set xrange [0:455]

set output "rmsf-avg.eps";
plot "WT_protein_system_rmsf_byres_avgd.dat" u ($1):($2) w lines t "Wild Type" ls 1 lw 4, \
"WT_protein_system_rmsf_byres_avgd.dat" u ($1):($2+$3):($2-$3) w filledcurve fs solid 0.15 t "Std. Dev." ls 1, \
"MUT_A_system_hbond_3trial_avgd.dat" u ($1):($2) w lines t "Mutation A" ls 2 lw 4, \
"MUT_A_system_hbond_3trial_avgd.dat" u ($1):($2+$3):($2-$3) w filledcurve fs solid 0.15 t "Std. Dev." ls 2, \
"WT_protein_system_rmsf_byres_avgd.dat" u ($1):($2) w lines notitle ls 1 lw 4, \
"MUT_A_system_hbond_3trial_avgd.dat" u ($1):($2) w lines notitle ls 2 lw 4;
