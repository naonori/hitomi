#! /bin/sh
#PBS -r y
#PBS -m n
#PBS -q mid
#PBS -o Log.out
#PBS -e Log.err
#PBS -N GalaxyData

cd $PBS_O_WORKDIR

#python fits2dat_galaxy.py NS zbin Weight
#python fits2dat_galaxy_random.py NS zbin Weight
#python fits2dat_mock.py South 3 0
python fits2dat_mock_random.py South 3 0

