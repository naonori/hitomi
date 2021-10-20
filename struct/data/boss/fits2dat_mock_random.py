########################
# Purpose of this code #
########################

# The purpose of this code is to calculate the comoving distance to mock randam particles from their redshifts and to rewrite the mock random distribution as 3D Cartesian coordinates.
 
# In order to calculate the comoving distances, we need to assume certain fiducial cosmological parameters. 
# The values of these parameters are NOT the same as those used for the BOSS galaxy data.

# As input parameters specific to this code, we need to set the "NS" parameter to select the NGC and SGC galaxy samples, 
# and the "zbin" parameter to determine the redshift bins. 

######################################################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
from astropy.io import fits

from classy import Class

## Choose NGC or SGC.
#NS = "North" # "North" or "South"
## Determine the three redshift bins #
#zbin = 1 # 1, 2, or 3
## Control the number of mock realizations to be loaded.
## NR = 0

args = sys.argv
NS = str(args[1])
zbin = int(args[2])
NR = int(args[3])

print("NS = ", NS)
print("zbin = ", zbin)

if zbin == 1:
    (zmin, zmax) = (0.2, 0.5)
elif zbin == 2:
    (zmin, zmax) = (0.4, 0.6)
elif zbin == 3:
    (zmin, zmax) = (0.5, 0.75)
else:
    print("ERROR")
    sys.exit()

if NS == "North":
    SG = "NGC"
elif NS == "South":
    SG = "SGC"
else:
    print("ERROR")
    sys.exit()

#####################################################################################################################
# Patchy-mock parameters #
# Note that these parameters are NOT the same as those used for the BOSS galaxy data.
#####################################################################################################################
h        = 0.6777
Omega_b  = 0.0480
Omega_m  = 0.3071
params_cosmo = {
    'h': h,
    'omega_b': Omega_b * h**2,
    'omega_cdm': (Omega_m - Omega_b) * h**2,
    'n_s': 0.9645,
    'ln10^{10}A_s': 3.094,
    'tau_reio': 0.0826026,
}

# load "class" #
cosmo = Class()
cosmo.set(params_cosmo)
cosmo.compute()

# Function for calculating the comoving distances to galaxies
def chi(z_bin):
    chi_bin = np.zeros(len(z_bin))
    if type(z_bin) == type(np.array([])):
        for i in range(len(z_bin)):
            chi_bin[i] = cosmo.angular_distance(z_bin[i]) * params_cosmo["h"] * (1.0 + z_bin[i])
    return chi_bin

input_fname = "Patchy-Mocks-Randoms-DR12%s-COMPSAM_V6C_x100" % SG
output_dir  = "Patchy-Mocks-DR12%s-COMPSAM_V6C_ZBIN%d" % (SG, zbin)
try:
	os.mkdir(output_dir)
except :
    print(output_dir, " already exists. ")
 
RA, DEC, Z, NBAR, W_f, W_c = np.loadtxt("%s.dat" % (input_fname), usecols=(0,1,2,3,5,6),unpack=True)
W = W_f * W_c / (1.0 + 10000.0 * NBAR)
		
# split into three redshift bins #
Z_bin   = (Z[Z>zmin])[Z[Z>zmin]<zmax]
RA_bin  = (RA[Z>zmin])[Z[Z>zmin]<zmax]
DEC_bin = (DEC[Z>zmin])[Z[Z>zmin]<zmax]
W_bin   = (W[Z>zmin])[Z[Z>zmin]<zmax]

# convert degrees to radians #
RA_bin  *= np.pi / 180.0
DEC_bin *= np.pi / 180.0

# calculation of comoving distances #
CHI_bin = chi(Z_bin)
 
# convert spherical coordinates to 3D Cartesian coordinates (x,y,z) #
X_dis = CHI_bin * np.cos(DEC_bin) * np.cos(RA_bin)
Y_dis = CHI_bin * np.cos(DEC_bin) * np.sin(RA_bin)
Z_dis = CHI_bin * np.sin(DEC_bin)

SAVE = np.array([X_dis, Y_dis, Z_dis, W_bin]).T

np.savetxt("%s/%s_ZBIN%d.dat" % (output_dir, input_fname, zbin), SAVE, fmt="%.10e")

