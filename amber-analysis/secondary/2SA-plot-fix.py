#!/usr/bin/env python3

import sys

# gnu_in:   the gnuplot file saved by cpptraj
# gnu_out:  the filename to save the fixed gnuplot file as
# pic:      name of saved gnuplot figure
# prot:     number of protein residues
# time:     length of simulation in ns
# divider:  dividend for frames v time (typically 100 Lang, 500 Berend)

## Use: python 2SA-plot-fix.py gnu_in gnu_out pic prot time divider
## Use: python 2SA-plot-fix.py test.gnu test_fix.gnu test.png 300 200 100
# Script is sys.argv[0]
gnu_read = sys.argv[1]
gnu_out = sys.argv[2]
pic = sys.argv[3]
prot = sys.argv[4]
time = sys.argv[5]
divider = sys.argv[6]

# Read the original gnuplot file line-by-line
gnu_in = open(gnu_read, 'r')
lines = gnu_in.readlines()
gnu_in.close()

# Change the final line from `pause -1` to `set output`
lines[-1] = "set output"

# Write out a new file with new header info
with open(gnu_out, 'w+') as f:
    f.write("set terminal pngcairo size 2560,1920 font \"Arial,48\";\n")
    f.write("set size 0.96,1\n")
    f.write("set encoding iso_8859_1\n")
    f.write("set pm3d map corners2color c1\n")
    f.write("set xtics nomirror out\n")
    f.write("#set ytics nomirror\n")
    # Set the explicit y-tics for your system if they're special
    f.write("set ytics (\"100\" 0, \"150\" 50, \"200\" 100, \"250\" 150, \"300\" 200, \\\n")
    f.write("\"350\" 250, \"400\" 300, \"BR\" 325, \"\" 335, \"1000\" 350, \"1050\" 400) \\\n")
    f.write("border nomirror out\n")
    f.write("set cbrange [  -0.500:   7.500]\n")
    f.write("set cbtics    0.000,7.000,1.0\n")
    f.write("set palette maxcolors 8\n")
    # Use the Paul Tol muted color scheme https://personal.sron.nl/~pault/
    f.write("set palette defined (0 \"#DDDDDD\",1 \"#AA4499\",2 \"#882255\", 3 \"#CC6677\",\\\n")
    f.write("4 \"#DDCC77\",5 \"#999933\", 6 \"#117733\",7 \"#44AA99\")\n")
    f.write("set cbtics(\"None\"    0.000,\"Para\"    1.000,\"Anti\"    2.000,\"3-10\"    \\\n")
    f.write("3.000,\"Alpha\"    4.000,\"Pi\"    5.000,\"Turn\"    6.000,\"Bend\"    7.000)\n")
    f.write("set xlabel \"Time (ns)\"\n")
    f.write("set ylabel \"Residue\"\n")
    f.write(f"set yrange [   0.000: {prot}]\n")
    f.write(f"set xrange [   0.000: {time}]\n")
    f.write(f'set output "{pic}"\n')
    f.write(f"splot \"-\" u ($1/{divider}):2:3 with pm3d notitle\n")
    # Rewrite the rest of the file
    # Skipping the first 13 lines with the original header info
    i = 0
    for line in lines:
        if i >= 13:
            f.write(line)
            i +=1
        else:
            i+=1
