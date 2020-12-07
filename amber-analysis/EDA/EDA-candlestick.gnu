## This script works with data generated from `rmagic-EDA-avg-diffs.r`

set encoding iso_8859_1
set term postscript enhanced color font "Arial,24";

set xlabel "Residue number"
set ylabel "Total Interaction Energy (kcal/mol)"
set xrange [0:450]
set key top left Left reverse width 2 height 1

set xtics ("100" 0, "150" 50, "200" 100, "250" 150, "300" 200, "350" 250, "450" 300, \
     "" 325, "" 335, "1000" 350, "1050" 400) border nomirror out;
set x2tics border offset 0,-0.25 nomirror out norotate left ("BR" 325, "" 335)

## For points use "w points" instead of "w boxes" (boxplots)
## For lines use "w lines"
## Candlesticks are box and whisker plots

## Plot once, plot standard deviations, then plot again so the data can
## actually be seen and the legend is ordered properly
set output "WT-MUTA_EDA_diff.eps";
plot "WT-MUTA_total_interaction_res220_avg.dat" u 1:2 w points t "WT - SysA" lw 4 pt 0 lc "black", \
"WT-MUTA_total_interaction_res220_avg.dat" u 1:2:($2-$3):($2+$3):2 w candlesticks fs solid 0.15 t "Avg. Std. Dev." ls 1,\
"WT-MUTA_total_interaction_res220_avg.dat" u 1:2 w points notitle ls 1 pt 0 lc "black" lw 4;
