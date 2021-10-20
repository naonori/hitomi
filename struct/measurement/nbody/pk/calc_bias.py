########################
## Purpose of this code
########################
## The code roughly estimates the value of the linear bias by calculating the ratio of the material power spectrum to the halo power spectrum at k=0.03 h/Mpc. Then, we use the value of the linear bias estimated here as an input parameter to perform the reconstruction. 
## We do not need to know the exact value of the linear bias because we can maintain a one-to-one correspondence between the observed results and the theoretical model by using the same value of the bias obtained here as the input parameter to compute the theoretical model.
##

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

b = []
p = []
for zbin in [0,1,2,3,4]:
    k, pk_m = np.loadtxt("results_test_Gadget_zbin%d_RSDFalse/0001/pk0.dat" % (zbin), unpack=True, usecols=(0,1))
    print("k = %3.3f [h/Mpc]" % k[2])
    for (Mmin, Mmax) in [(12.5, 13.0), (13.0, 13.5), (13.5, 14.0)]:
        k, pk_h = np.loadtxt("results_test_Rockstar_Mmin%s_Mmax%s_zbin%d_RSDFalse/0001/pk0.dat" % (Mmin, Mmax, zbin), unpack=True, usecols=(0,1))
        b1 = np.sqrt(pk_h[2]/pk_m[2])
        p.append("zbin%d_Mmin%2.1f_Mmax%2.1f" % (zbin, Mmin, Mmax))
        b.append(b1)

fw = open("bias.dat", "w")
for i in range(len(b)):
    fw.write("%s  %1.3f\n" % (p[i], b[i]))
fw.close()

