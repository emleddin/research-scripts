#!/bin/bash

## Set program location (can be in .bashrc instead)
export CAVER="/home/$USER/bin/caver_3.0/caver"

## Areas of red conservation
##
## 601 is metal center
## 603 is cofactor

listVar=(601 603)

## Backup your original config file
cp config.txt config-start.txt

for RESNUM in "${listVar[@]}"
do
        ## Use double quotes so the shell can substitute variables
        ## Replace the full line with "starting_point_residue $RESNUM"
        sed "22c\starting_point_residue $RESNUM" config-start.txt > config.txt

        ## Run the Caver program
        java -Xmx10000m -cp $CAVER/lib -jar $CAVER/caver.jar -home $CAVER \
         -pdb ./ -conf ./config.txt -out ./caver_output_res$RESNUM

        ## Copy the config file used to the output folder
        cp config.txt ./caver_output_res$RESNUM/

done
