#!/bin/bash
## Written for GNU sed (MacOS needs any -i <> to be -i '' <>)

## Name designations, without file extensions
## You can use paths
gnu_files="sys1_secstruct.gnu
sys2_secstruct.gnu
sys3_secstruct.gnu"

out_pngs="sys1_secstruct.png
sys2_secstruct.png
sys3_secstruct.png"

## Loop through gnu_files and out_pngs, modifying the gnu_files appropriately
clean_gnuplot_files()
{
set $out_pngs
for gnu_file in $gnu_files; do

## Print out which fileset of all of them that you're on
echo Gnuplot: "$gnu_file" PNG: "$1"

## new_lines contains the updated header information for the gnuplot files.
## Update `xrange`, `yrange`, and the `splot (\$1/100)`variables  with what
## you need for your system.
## This is set up for custom ytics, but you can toggle the commenting
## if you don't have an infuriatingly numbered system.
## This uses a muted color scheme from https://personal.sron.nl/~pault/
## but you can change it (`set palette defined`).
## Use double quotes around everything so that variables get evaluated,
## with the caveat that all other double quotes must be escaped
## (looking at you, ytics).
## The $1 refers to $out_pngs being set in fun().
## You need \\\\ for the proper number of escapes to print a single \ with ex
new_lines="set terminal pngcairo size 2560,1920 font \"Helvetica,48\";
set size 0.96,1
set encoding iso_8859_1
set pm3d map corners2color c1
set xtics nomirror out
#set ytics nomirror
set ytics (\"100\" 0, \"150\" 50, \"200\" 100, \"250\" 150, \"300\" 200, \\\\
\"350\" 250, \"400\" 300, \"BR\" 325, \"\" 335, \"1000\" 350, \"1050\" 400) \\\\
border nomirror out
set cbrange [  -0.500:   7.500]
set cbtics    0.000,7.000,1.0
set palette maxcolors 8
set palette defined (0 \"#DDDDDD\",1 \"#AA4499\",2 \"#882255\", 3 \"#CC6677\",\\\\
4 \"#DDCC77\",5 \"#999933\", 6 \"#117733\",7 \"#44AA99\")
set cbtics(\"None\"    0.000,\"Para\"    1.000,\"Anti\"    2.000,\"3-10\"    \\\\
3.000,\"Alpha\"    4.000,\"Pi\"    5.000,\"Turn\"    6.000,\"Bend\"    7.000)
set xlabel \"Time (ns)\"
set ylabel \"Residue\"
set yrange [   0.000: 430.000]
set xrange [   0.000: 200.000]
set output \"${1}\"
splot \"-\" u (\$1/100):2:3 with pm3d notitle"

## ex is a command-line version of vi -- the << eof tells it to wait until eof
## :1,13d deletes the first 13 lines (bad header)
## :%s line sets up gnuplot for scripts (pause -1 assumes interactive gnuplot)
## 1 insert will insert before the 1st line (inserts everything until . given)
## $new_lines is what gets inserted (the {%?} evaluates it)
ex ${gnu_file} << eof
:1,13d
:%s/pause -1/set output/g
1 insert
${new_lines%?}
.
:wq
eof

## Increment $out_pngs and exit the function
shift
done
}

## Set up a function to go through the gnuplot scripts
run_gnuplot()
{
  for gnu_file in $gnu_files; do
    ## Print a status report
    echo Processing ${gnu_file} now
    gnuplot ${gnu_file}
  done
}

## Run the function
clean_gnuplot_files

## Run the gnuplot scripts too, while we're at it
run_gnuplot
