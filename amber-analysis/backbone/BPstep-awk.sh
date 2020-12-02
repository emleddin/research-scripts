#!/bin/bash

# BP Step
awk 'NR == 1 || NR % 4 == 2' BPstep.master-10.dat > BPstep-432-433.dat
awk 'NR == 1 || NR % 4 == 3' BPstep.master-10.dat > BPstep-433-434.dat
awk 'NR == 1 || NR % 4 == 0' BPstep.master-10.dat > BPstep-434-435.dat


# BP Helix
awk 'NR == 1 || NR % 4 == 2' Helix.master-10.dat > Helix-432-433.dat
awk 'NR == 1 || NR % 4 == 3' Helix.master-10.dat > Helix-433-434.dat
awk 'NR == 1 || NR % 4 == 0' Helix.master-10.dat > Helix-434-435.dat

# Matches
awk 'NR == 1 || NR % 5 == 2' BP.master-10.dat > BP-432-433.dat
awk 'NR == 1 || NR % 5 == 3' BP.master-10.dat > BP-433-434.dat
awk 'NR == 1 || NR % 5 == 4' BP.master-10.dat > BP-434-435.dat
