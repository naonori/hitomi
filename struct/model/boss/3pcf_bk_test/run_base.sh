#! /bin/sh
#PBS -r y
#PBS -m n
#PBS -q mid
#PBS -N Bk3PCF

cd $PBS_O_WORKDIR

python calc_3pcf_bk_test.py 0 2 2 > LOG/log022_decom

