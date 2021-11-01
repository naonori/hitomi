########################
# Purpose of this code #
########################
# The code uses a number of 2PCFs and 3PCFs measured from the Patchy mock catalogs and computes their covariance matrices. In addition, this code adds the statistical errors to the 2PCFs and 3PCFs measured from galaxy and mock data, and outputs them. All 2PCFs and 3PCFs are normalized by their respective monopole window functions, xi0_window.dat and zeta000_window.dat.
#####################################################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import numpy as np
from distutils.util import strtobool

WORK = "mwork0/sugiymnn/WORK"
MEASUREMENT = "measurement/boss"
MOCK_2PCF = "2pcf_mock"
MOCK_3PCF = "3pcf_mock"
GALAXY_2PCF = "2pcf"
GALAXY_3PCF = "3pcf"

# Galaxy sample #
NS = "North"
zbin = 1
#NS = sys.argv[1] # North, South
#zbin = int(sys.argv[2]) # 1, 3

print(NS, "zbin = %d" % zbin)

# Reconstruction #
recon = False # True or False
#recon = strtobool(sys.argv[3]) # True or False
print("reconstruction = %s" % recon)

# Directories containing the 2PCF and 3PCF data. #
if recon == False:
    DATA_2PCF = "/%s/%s/%s/mock_%s_zbin%d" % (WORK, MEASUREMENT, MOCK_2PCF, NS, zbin)
    DATA_2PCF_WINDOW = "/%s/%s/%s/mock_%s_zbin%d_Window" % (WORK, MEASUREMENT, MOCK_2PCF, NS, zbin)
    GALAXY_DATA_2PCF = "/%s/%s/%s/galaxy_%s_zbin%d" % (WORK, MEASUREMENT, GALAXY_2PCF, NS, zbin)
    GALAXY_DATA_2PCF_W1 = "/%s/%s/%s/galaxy_%s_zbin%d_Weight1_OnlySys" % (WORK, MEASUREMENT, GALAXY_2PCF, NS, zbin)
    GALAXY_DATA_2PCF_W2 = "/%s/%s/%s/galaxy_%s_zbin%d_Weight2_NoSys" % (WORK, MEASUREMENT, GALAXY_2PCF, NS, zbin)
    GALAXY_DATA_2PCF_W3 = "/%s/%s/%s/galaxy_%s_zbin%d_Weight3_NoWeight" % (WORK, MEASUREMENT, GALAXY_2PCF, NS, zbin)
    GALAXY_DATA_2PCF_WINDOW = "/%s/%s/%s/galaxy_%s_zbin%d_Window" % (WORK, MEASUREMENT, GALAXY_2PCF, NS, zbin)
elif recon == True:
    DATA_2PCF = "/%s/%s/%s/mock_%s_zbin%d_recon_R15" % (WORK, MEASUREMENT, MOCK_2PCF, NS, zbin)
    DATA_2PCF_WINDOW = "/%s/%s/%s/mock_%s_zbin%d_recon_R15_Window" % (WORK, MEASUREMENT, MOCK_2PCF, NS, zbin)
    GALAXY_DATA_2PCF = "/%s/%s/%s/galaxy_%s_zbin%d_recon_R15" % (WORK, MEASUREMENT, GALAXY_2PCF, NS, zbin)
    GALAXY_DATA_2PCF_W1 = "/%s/%s/%s/galaxy_%s_zbin%d_recon_R15_Weight1_OnlySys" % (WORK,MEASUREMENT,GALAXY_2PCF,NS,zbin)
    GALAXY_DATA_2PCF_W2 = "/%s/%s/%s/galaxy_%s_zbin%d_recon_R15_Weight2_NoSys" % (WORK,MEASUREMENT,GALAXY_2PCF,NS,zbin)
    GALAXY_DATA_2PCF_W3 = "/%s/%s/%s/galaxy_%s_zbin%d_recon_R15_Weight3_NoWeight"%(WORK,MEASUREMENT,GALAXY_2PCF,NS,zbin)
    GALAXY_DATA_2PCF_WINDOW = "/%s/%s/%s/galaxy_%s_zbin%d_recon_R15_Window" % (WORK, MEASUREMENT, GALAXY_2PCF, NS, zbin)

DATA_3PCF = "/%s/%s/%s/mock_%s_zbin%d" % (WORK, MEASUREMENT, MOCK_3PCF, NS, zbin)
DATA_3PCF_WINDOW = "/%s/%s/%s/mock_%s_zbin%d_Window" % (WORK, MEASUREMENT, MOCK_3PCF, NS, zbin)
GALAXY_DATA_3PCF = "/%s/%s/%s/galaxy_%s_zbin%d" % (WORK, MEASUREMENT, GALAXY_3PCF, NS, zbin)
GALAXY_DATA_3PCF_W1 = "/%s/%s/%s/galaxy_%s_zbin%d_Weight1_OnlySys" % (WORK, MEASUREMENT, GALAXY_3PCF, NS, zbin)
GALAXY_DATA_3PCF_W2 = "/%s/%s/%s/galaxy_%s_zbin%d_Weight2_NoSys" % (WORK, MEASUREMENT, GALAXY_3PCF, NS, zbin)
GALAXY_DATA_3PCF_W3 = "/%s/%s/%s/galaxy_%s_zbin%d_Weight3_NoWeight" % (WORK, MEASUREMENT, GALAXY_3PCF, NS, zbin)
GALAXY_DATA_3PCF_WINDOW = "/%s/%s/%s/galaxy_%s_zbin%d_Window" % (WORK, MEASUREMENT, GALAXY_3PCF, NS, zbin)

# Number of bins for each of the 2PCF and 3PCF.
num_rbin_2pcf = 25
num_rbin_3pcf = 13

# Total number of mock catalogs #
Rmax = 2048

# Reading 2PCF #
xi_R_dict = {}
for ELL in [0,2,4]:
    xi_R = np.zeros((num_rbin_2pcf, Rmax)) 
    # IMPORTANT: reading the window function #
    if recon == False:
        w0 = np.loadtxt("%s/xi0_window.dat" % (DATA_2PCF_WINDOW), usecols=(1,), unpack=True)
    elif recon == True:
        w0 = np.loadtxt("%s/xi0_recon_window.dat" % (DATA_2PCF_WINDOW), usecols=(1,), unpack=True)
    
    for realization in range(1, Rmax+1):
        print("realization = ", realization)
        # Normalized by w0 #
        if recon == False:
            xi_R[:, realization-1] = np.loadtxt("%s/%04d/xi%d.dat" % (DATA_2PCF, realization, ELL), usecols=(1,), unpack=True)[:] / w0[:]
        if recon == True:
            xi_R[:, realization-1] = np.loadtxt("%s/%04d/xi%d_recon.dat" % (DATA_2PCF, realization, ELL), usecols=(1,), unpack=True)[:] / w0[:]

    xi_R_dict.update({"xi%d" % (ELL): xi_R})

# Reading 3PCF #
zeta_R_dict = {}
for (ell1, ell2, ELL) in [(0,0,0),(1,1,0),(2,0,2),(1,1,2)]:
    zeta_R = np.zeros((num_rbin_3pcf**2, Rmax)) 
    for realization in range(1, Rmax+1):
        print("realization = ", realization)
        zeta_temp = np.zeros((num_rbin_3pcf,num_rbin_3pcf))
        for i in range(num_rbin_3pcf):
            # IMPORTANT: reading the window function #
            w000 = np.loadtxt("%s/zeta000_%02d_window.dat" % (DATA_3PCF_WINDOW, i), usecols=(2,), unpack=True)
            # Normalized by w000 #
            zeta_temp[i,:] = np.loadtxt("%s/%04d/zeta%d%d%d_%02d.dat" % (DATA_3PCF, realization, ell1, ell2, ELL, i), usecols=(2,), unpack=True)[:] / w000[:]
        
        # Rewrite the "zeta" represented by the two-dimensional array of r1 and r2 into a one-dimensional array. #
        zeta_R[:, realization-1] = zeta_temp.reshape(num_rbin_3pcf**2)[:]
    
    zeta_R_dict.update({"zeta%d%d%d" % (ell1, ell2, ELL): zeta_R})
 

# Reading rbin for 2pcf #
rbin_2pcf = np.zeros(num_rbin_2pcf)
if recon == False:
    rbin_2pcf = np.loadtxt("%s/%04d/xi%d.dat" % (DATA_2PCF, 1, 0), usecols=(0,), unpack=True)
elif recon == True:
    rbin_2pcf = np.loadtxt("%s/%04d/xi%d_recon.dat" % (DATA_2PCF, 1, 0), usecols=(0,), unpack=True)

# Reading rbin for 3pcf #
rbin_3pcf = np.zeros((num_rbin_3pcf, num_rbin_3pcf))
for i in range(num_rbin_3pcf):
    rbin_3pcf[i,:] = np.loadtxt("%s/%04d/zeta%d%d%d_%02d.dat" % (DATA_3PCF, 1, 0, 0, 0, i), usecols=(1,), unpack=True)

rbin_3pcf_ith = (rbin_3pcf.T).reshape(num_rbin_3pcf**2)
rbin_3pcf_jth = (rbin_3pcf).reshape(num_rbin_3pcf**2)

# Make an output directory #
if recon == False:
    OUTPUT = "data_and_cov_%s_zbin%d" % (NS, zbin)
elif recon == True:
    OUTPUT = "data_and_cov_%s_zbin%d_recon_R15" % (NS, zbin)

try:
    os.mkdir(OUTPUT)
except:
    print("")


import pickle
# Save the 2PCF measured from each of the mock catalogs #
with open("%s/xi_R.p" % OUTPUT, "wb") as pp:
    pickle.dump(xi_R_dict, pp)

# Save the 3PCF measured from each of the mock catalogs #
with open("%s/zeta_R.p" % OUTPUT, "wb") as pp:
    pickle.dump(zeta_R_dict, pp)

# Save the mean of the 2PCF measured from the mock data #
for ELL in [0,2,4]:
    mean_xi = np.mean(xi_R_dict["xi%d" % (ELL)], axis=1)
    std_xi  = np.std(xi_R_dict["xi%d" % (ELL)], axis=1, ddof=1)
    np.savetxt("%s/xi%d_mean.dat" % (OUTPUT, ELL), np.array([rbin_2pcf, mean_xi, std_xi]).T, fmt="%.5f \t %.7e \t %.7e")

# Save the 2PCF measured from the BOSS DR12 data #
for ELL in [0,2,4]:
    std_xi  = np.std(xi_R_dict["xi%d" % (ELL)], axis=1, ddof=1)
    # IMPORTANT: reading the window function #
    if recon == False:
        w0 = np.loadtxt("%s/xi0_window.dat" % (GALAXY_DATA_2PCF_WINDOW), usecols=(1,), unpack=True)
        xi = np.loadtxt("%s/xi%d.dat" % (GALAXY_DATA_2PCF, ELL), usecols=(1,), unpack=True) / w0
    elif recon == True:
        w0 = np.loadtxt("%s/xi0_recon_window.dat" % (GALAXY_DATA_2PCF_WINDOW), usecols=(1,), unpack=True)
        xi = np.loadtxt("%s/xi%d_recon.dat" % (GALAXY_DATA_2PCF, ELL), usecols=(1,), unpack=True) / w0
    # Normalized by w0 #
    np.savetxt("%s/xi%d_galaxy.dat" % (OUTPUT, ELL), np.array([rbin_2pcf, xi, std_xi]).T, fmt="%.5f \t %.7e \t %.7e")

for ELL in [0,2,4]:
    std_xi  = np.std(xi_R_dict["xi%d" % (ELL)], axis=1, ddof=1)
    # IMPORTANT: reading the window function #
    if recon == False:
        w0 = np.loadtxt("%s/xi0_window.dat" % (GALAXY_DATA_2PCF_WINDOW), usecols=(1,), unpack=True)
        xi = np.loadtxt("%s/xi%d.dat" % (GALAXY_DATA_2PCF_W1, ELL), usecols=(1,), unpack=True) / w0
    elif recon == True:
        w0 = np.loadtxt("%s/xi0_recon_window.dat" % (GALAXY_DATA_2PCF_WINDOW), usecols=(1,), unpack=True)
        xi = np.loadtxt("%s/xi%d_recon.dat" % (GALAXY_DATA_2PCF_W1, ELL), usecols=(1,), unpack=True) / w0
    # Normalized by w0 #
    np.savetxt("%s/xi%d_galaxy_W1.dat" % (OUTPUT, ELL), np.array([rbin_2pcf, xi, std_xi]).T, fmt="%.5f \t %.7e \t %.7e")


for ELL in [0,2,4]:
    std_xi  = np.std(xi_R_dict["xi%d" % (ELL)], axis=1, ddof=1)
    # IMPORTANT: reading the window function #
    if recon == False:
        w0 = np.loadtxt("%s/xi0_window.dat" % (GALAXY_DATA_2PCF_WINDOW), usecols=(1,), unpack=True)
        xi = np.loadtxt("%s/xi%d.dat" % (GALAXY_DATA_2PCF_W2, ELL), usecols=(1,), unpack=True) / w0
    elif recon == True:
        w0 = np.loadtxt("%s/xi0_recon_window.dat" % (GALAXY_DATA_2PCF_WINDOW), usecols=(1,), unpack=True)
        xi = np.loadtxt("%s/xi%d_recon.dat" % (GALAXY_DATA_2PCF_W2, ELL), usecols=(1,), unpack=True) / w0
    # Normalized by w0 #
    np.savetxt("%s/xi%d_galaxy_W2.dat" % (OUTPUT, ELL), np.array([rbin_2pcf, xi, std_xi]).T, fmt="%.5f \t %.7e \t %.7e")


for ELL in [0,2,4]:
    std_xi  = np.std(xi_R_dict["xi%d" % (ELL)], axis=1, ddof=1)
    # IMPORTANT: reading the window function #
    if recon == False:
        w0 = np.loadtxt("%s/xi0_window.dat" % (GALAXY_DATA_2PCF_WINDOW), usecols=(1,), unpack=True)
        xi = np.loadtxt("%s/xi%d.dat" % (GALAXY_DATA_2PCF_W3, ELL), usecols=(1,), unpack=True) / w0
    elif recon == True:
        w0 = np.loadtxt("%s/xi0_recon_window.dat" % (GALAXY_DATA_2PCF_WINDOW), usecols=(1,), unpack=True)
        xi = np.loadtxt("%s/xi%d_recon.dat" % (GALAXY_DATA_2PCF_W3, ELL), usecols=(1,), unpack=True) / w0
    # Normalized by w0 #
    np.savetxt("%s/xi%d_galaxy_W3.dat" % (OUTPUT, ELL), np.array([rbin_2pcf, xi, std_xi]).T, fmt="%.5f \t %.7e \t %.7e")


# Save the mean of the 3PCF measured from the mock data #
for (ell1, ell2, ELL) in [(0,0,0),(1,1,0),(2,0,2),(1,1,2)]:
    mean_zeta = np.mean(zeta_R_dict["zeta%d%d%d" % (ell1, ell2, ELL)], axis=1)
    std_zeta  = np.std(zeta_R_dict["zeta%d%d%d" % (ell1, ell2, ELL)], axis=1, ddof=1)
    np.savetxt("%s/zeta%d%d%d_mean.dat" % (OUTPUT, ell1, ell2, ELL),\
                np.array([rbin_3pcf_ith, rbin_3pcf_jth, mean_zeta, std_zeta]).T,\
                fmt="%.5f \t %.5f \t %.7e \t %.7e")

# Save the 3PCF measured from the BOSS DR12 data #
for (ell1, ell2, ELL) in [(0,0,0),(1,1,0),(2,0,2),(1,1,2)]:
    std_zeta  = np.std(zeta_R_dict["zeta%d%d%d" % (ell1, ell2, ELL)], axis=1, ddof=1)
    zeta = np.zeros((num_rbin_3pcf,num_rbin_3pcf))
    for i in range(num_rbin_3pcf):
        # IMPORTANT: reading the window function #
        w000 = np.loadtxt("%s/zeta000_%02d_window.dat" % (GALAXY_DATA_3PCF_WINDOW, i), usecols=(2,), unpack=True)
        # Normalized by w000 #
        zeta[i,:] = np.loadtxt("%s/zeta%d%d%d_%02d.dat" % (GALAXY_DATA_3PCF, ell1, ell2, ELL, i), usecols=(2,), unpack=True)[:] / w000[:]
        
    # Rewrite the "zeta" represented by the two-dimensional array of r1 and r2 into a one-dimensional array. #
    zeta = zeta.reshape(num_rbin_3pcf**2)
    
    np.savetxt("%s/zeta%d%d%d_galaxy.dat" % (OUTPUT, ell1, ell2, ELL),\
                np.array([rbin_3pcf_ith, rbin_3pcf_jth, zeta, std_zeta]).T,\
                fmt="%.5f \t %.5f \t %.7e \t %.7e")

###

for (ell1, ell2, ELL) in [(0,0,0),(1,1,0),(2,0,2),(1,1,2)]:
    std_zeta  = np.std(zeta_R_dict["zeta%d%d%d" % (ell1, ell2, ELL)], axis=1, ddof=1)
    zeta = np.zeros((num_rbin_3pcf,num_rbin_3pcf))
    for i in range(num_rbin_3pcf):
        # IMPORTANT: reading the window function #
        w000 = np.loadtxt("%s/zeta000_%02d_window.dat" % (GALAXY_DATA_3PCF_WINDOW, i), usecols=(2,), unpack=True)
        # Normalized by w000 #
        zeta[i,:] = np.loadtxt("%s/zeta%d%d%d_%02d.dat" % (GALAXY_DATA_3PCF_W1, ell1, ell2, ELL, i), usecols=(2,), unpack=True)[:] / w000[:]
        
    # Rewrite the "zeta" represented by the two-dimensional array of r1 and r2 into a one-dimensional array. #
    zeta = zeta.reshape(num_rbin_3pcf**2)
    
    np.savetxt("%s/zeta%d%d%d_galaxy_W1.dat" % (OUTPUT, ell1, ell2, ELL),\
                np.array([rbin_3pcf_ith, rbin_3pcf_jth, zeta, std_zeta]).T,\
                fmt="%.5f \t %.5f \t %.7e \t %.7e")


for (ell1, ell2, ELL) in [(0,0,0),(1,1,0),(2,0,2),(1,1,2)]:
    std_zeta  = np.std(zeta_R_dict["zeta%d%d%d" % (ell1, ell2, ELL)], axis=1, ddof=1)
    zeta = np.zeros((num_rbin_3pcf,num_rbin_3pcf))
    for i in range(num_rbin_3pcf):
        # IMPORTANT: reading the window function #
        w000 = np.loadtxt("%s/zeta000_%02d_window.dat" % (GALAXY_DATA_3PCF_WINDOW, i), usecols=(2,), unpack=True)
        # Normalized by w000 #
        zeta[i,:] = np.loadtxt("%s/zeta%d%d%d_%02d.dat" % (GALAXY_DATA_3PCF_W2, ell1, ell2, ELL, i), usecols=(2,), unpack=True)[:] / w000[:]
        
    # Rewrite the "zeta" represented by the two-dimensional array of r1 and r2 into a one-dimensional array. #
    zeta = zeta.reshape(num_rbin_3pcf**2)
    
    np.savetxt("%s/zeta%d%d%d_galaxy_W2.dat" % (OUTPUT, ell1, ell2, ELL),\
                np.array([rbin_3pcf_ith, rbin_3pcf_jth, zeta, std_zeta]).T,\
                fmt="%.5f \t %.5f \t %.7e \t %.7e")

for (ell1, ell2, ELL) in [(0,0,0),(1,1,0),(2,0,2),(1,1,2)]:
    std_zeta  = np.std(zeta_R_dict["zeta%d%d%d" % (ell1, ell2, ELL)], axis=1, ddof=1)
    zeta = np.zeros((num_rbin_3pcf,num_rbin_3pcf))
    for i in range(num_rbin_3pcf):
        # IMPORTANT: reading the window function #
        w000 = np.loadtxt("%s/zeta000_%02d_window.dat" % (GALAXY_DATA_3PCF_WINDOW, i), usecols=(2,), unpack=True)
        # Normalized by w000 #
        zeta[i,:] = np.loadtxt("%s/zeta%d%d%d_%02d.dat" % (GALAXY_DATA_3PCF_W3, ell1, ell2, ELL, i), usecols=(2,), unpack=True)[:] / w000[:]
        
    # Rewrite the "zeta" represented by the two-dimensional array of r1 and r2 into a one-dimensional array. #
    zeta = zeta.reshape(num_rbin_3pcf**2)
    
    np.savetxt("%s/zeta%d%d%d_galaxy_W3.dat" % (OUTPUT, ell1, ell2, ELL),\
                np.array([rbin_3pcf_ith, rbin_3pcf_jth, zeta, std_zeta]).T,\
                fmt="%.5f \t %.5f \t %.7e \t %.7e")

##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################

# Save the auto-covariane matrix of the 2PCF #
for ELL in [0,2,4]:
    for ELL_d in [0,2,4]:
        cov = np.cov(xi_R_dict["xi%d" % (ELL)], xi_R_dict["xi%d" % (ELL_d)])[0:num_rbin_2pcf, num_rbin_2pcf:2*num_rbin_2pcf]
        np.savetxt("%s/cov_xi%d_xi%d.dat" % (OUTPUT, ELL, ELL_d), cov)
 
# Save the cross-covariane matrix between the 2pcf and the 3pcf #
for ELL in [0,2,4]:
    for (ell1_d, ell2_d, ELL_d) in [(0,0,0),(1,1,0),(2,0,2),(1,1,2)]:
        cov = np.cov(xi_R_dict["xi%d" % (ELL)], zeta_R_dict["zeta%d%d%d" % (ell1_d, ell2_d, ELL_d)])[0:num_rbin_2pcf, num_rbin_2pcf:num_rbin_2pcf+num_rbin_3pcf**2]
        np.savetxt("%s/cov_xi%d_zeta%d%d%d.dat" % (OUTPUT, ELL, ell1_d, ell2_d, ELL_d), cov)
 
# Save the cross-covariane matrix between the 2pcf and the 3pcf #
for (ell1, ell2, ELL) in [(0,0,0),(1,1,0),(2,0,2),(1,1,2)]:
    for ELL_d in [0,2,4]:
        cov = np.cov(zeta_R_dict["zeta%d%d%d" % (ell1, ell2, ELL)], xi_R_dict["xi%d" % (ELL_d)])[0:num_rbin_3pcf**2, num_rbin_3pcf**2:num_rbin_3pcf**2+num_rbin_2pcf]
        np.savetxt("%s/cov_zeta%d%d%d_xi%d.dat" % (OUTPUT, ell1, ell2, ELL, ELL_d), cov)
 
# Save the auto-covariane matrix of the 3pcf #
for (ell1, ell2, ELL) in [(0,0,0),(1,1,0),(2,0,2),(1,1,2)]:
    for (ell1_d, ell2_d, ELL_d) in [(0,0,0),(1,1,0),(2,0,2),(1,1,2)]:
        cov = np.cov(zeta_R_dict["zeta%d%d%d" % (ell1, ell2, ELL)],\
                     zeta_R_dict["zeta%d%d%d" % (ell1_d, ell2_d, ELL_d)])[0:num_rbin_3pcf**2, num_rbin_3pcf**2:2*num_rbin_3pcf**2]
        np.savetxt("%s/cov_zeta%d%d%d_zeta%d%d%d.dat" % (OUTPUT, ell1, ell2, ELL, ell1_d, ell2_d, ELL_d), cov)


