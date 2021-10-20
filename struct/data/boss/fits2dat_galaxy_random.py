########################
# Purpose of this code #
########################
 
# The purpose of this code is to calculate the comoving distance to galaxy random particles from their redshifts and to rewrite the random distribution as 3D Cartesian coordinates.

# In order to calculate the comoving distances, we need to assume certain fiducial cosmological parameters. 
# The values of these parameters must be exactly the same as those used to calculate the theoretical model templates.

# The true cosmological parameters, which we do not know, must be different from the fiducial parameters, and these differences can produce the AP effect.

# As input parameters specific to this code, we need to set the "NS" parameter to select the NGC and SGC galaxy samples, 
# and the "zmin" and "zmax" parameters to determine the redshift bins. 

# The random particles have no observational weights.

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

args = sys.argv
NS = str(args[1])
zbin = int(args[2])
Weight = int(args[3])

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

#####################################################################################################################
# BOSS fiducial parameters #
# The values of the fiducial parameters should be exactly the same as those used to calculate the theoretical model.
#####################################################################################################################
params_cosmo = {
    'h': 0.676,
    'omega_b': 0.022,
    'omega_cdm': 0.11966256,
    'n_s': 0.96,
    'ln10^{10}A_s': 3.094,
    'tau_reio': 0.0826026,
    
    'N_ur': 2.0328,
    'N_ncdm': 1,
    'm_ncdm': 0.06,
    'T_ncdm': 0.71611,
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

# Loading galaxy random data #
with fits.open("random0_DR12v5_CMASSLOWZTOT_%s.fits" % (NS)) as hdul0:
	data0 = hdul0[1].data

with fits.open("random1_DR12v5_CMASSLOWZTOT_%s.fits" % (NS)) as hdul1:
	data1 = hdul1[1].data

RA = np.hstack([data0["RA"], data1["RA"]])
DEC = np.hstack([data0["DEC"], data1["DEC"]])
Z = np.hstack([data0["Z"], data1["Z"]])
W = np.hstack([data0["WEIGHT_FKP"], data1["WEIGHT_FKP"]])

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

# save data #
output_fname = "random_DR12v5_CMASSLOWZTOT_%s_ZBIN%d.dat" % (NS,zbin)
output_dir = "galaxy_DR12v5_CMASSLOWZTOT"
try:
    os.mkdir(output_dir)
except:
    print(output_dir, " already exists ")
 
np.savetxt("%s/%s" % (output_dir, output_fname), SAVE, fmt="%.10e")


