#!/bin/sh
#PBS -N 3PCF
#PBS -l nodes=1 
#PBS -j oe
#PBS -q bulk-b 

cd $PBS_O_WORKDIR

./a.out default_param.ini > log

