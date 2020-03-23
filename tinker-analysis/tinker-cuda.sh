#!/bin/bash
#PBS -q my_gpu_alloc
#PBS -l nodes=n01-02-03
#PBS -j oe
#PBS -r n
#PBS -o err.error
#PBS -N WT_protein_system

xyz=WT_protein_system_prod.xyz
key=WT_protein_system_prod.key

export CUDA_VISIBLE_DEVICES=0

cd $PBS_O_WORKDIR

module load tinker/8.7.1/gcc550-cuda8

arc=`echo "$xyz" | cut -d'.' -f1`

## natom + 2 = linenumbers
## read trims whitespace/carriage returns of first line
## natom cut processes number when title present
read -r var < ${xyz}
natom=`echo "$var" | cut -d' ' -f1`
linenumber=$((natom+2))

## md = MD steps
## ts = timestep in fs
## int = time interval for printing in ps
## sim = simulation ensemble (2=NVT)
## temp = desired temperature in K
## GPU = print GPU info
md=250000
ts=2.0
int=0.5
sim=2
temp=300
GPU=N

e=1
f=2

## Copy the first one for looping
if [$e -eq 1]
then
    cp $xyz ${xyz}_${e}
fi

## 500 ps = 0.5 ns = 1 file; 20 files is 10 ns
while [ $e -lt 21 ]; do

dynamic_omm ${xyz}_${e} -k $key $md $ts $int $sim $temp $GPU > prod$e.out

## Resave the arc file (otherwise overwritten)
mv ${arc}.arc ${arc}_${e}.arc

## Save the next input XYZ
tail -n ${linenumber} ${arc}_${e}.arc > ${xyz}_${f}

e=$[$e+1]
f=$[$f+1]

done
