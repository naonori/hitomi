#! /bin/sh
#PBS -r y
#PBS -m n
#PBS -q mid
#PBS -N GalaxyData

cd $PBS_O_WORKDIR

python MontePython.py run -o chains_test -p input/base_galaxy_zbin3.param --conf=my.conf --superupdate 20 -N 200000



#python analysis.py

