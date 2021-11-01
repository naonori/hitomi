#! /bin/sh
#PBS -r y
#PBS -m n
#PBS -q mid
#PBS -N GalaxyData

cd $PBS_O_WORKDIR

python calc_3pcf.py

