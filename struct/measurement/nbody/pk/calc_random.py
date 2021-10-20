########################
## Purpose of this code
########################
## The code generates a uniform random particle distribution to be used when performing the BAO reconstruction. The number of random particles is 100 times the number of data particles such as haloes.
###

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
import random

realization = 1
random_dir = "/mwork0/sugiymnn/WORK/data/nbody/%04d/"  % realization + "rockstar/random"

try:
    os.mkdir(random_dir)
except:
    print("")

for zbin in [0,1,2,3,4]:
    for (Mmin, Mmax) in [(12.5, 13.0), (13.0, 13.5), (13.5, 14.0)]:

        print("zbin =", zbin)
        print("Mmin, Mmax = ", Mmin, Mmax)

        fr = open("results_test_Rockstar_Mmin%s_Mmax%s_zbin%d_RSDFalse/0001/simulation_parameters" % (Mmin, Mmax, zbin))
        AA = fr.readlines()
        fr.close()
        
        print(AA[1])
        print(AA[5])
        
        n_tot = int(AA[1].split()[2])
        boxsize = float(AA[5].split()[3])
        
        n_tot_random = 100 * n_tot
        
        x = np.zeros(n_tot_random)
        y = np.zeros(n_tot_random)
        z = np.zeros(n_tot_random)
        for p in range(n_tot_random):
            x[p] = random.uniform(0, boxsize)
            y[p] = random.uniform(0, boxsize)
            z[p] = random.uniform(0, boxsize)
        
        X = np.array([x, y, z]).T
        np.savetxt("%s/random_Mmin%s_Mmax%s_zbin%d.dat" % (random_dir, Mmin, Mmax, zbin), X, fmt="%.7f")
        
        
