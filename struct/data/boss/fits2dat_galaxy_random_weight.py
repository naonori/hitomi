########################
# Purpose of this code #
########################
 
# The purpose of this code is to calculate the comoving distance to galaxy random particles from their redshifts and to rewrite the random distribution as 3D Cartesian coordinates.

# In order to calculate the comoving distances, we need to assume certain fiducial cosmological parameters. 
# The values of these parameters must be exactly the same as those used to calculate the theoretical model templates.

# The true cosmological parameters, which we do not know, must be different from the fiducial parameters, and these differences can produce the AP effect.

# As input parameters specific to this code, we need to set the "NS" parameter to select the NGC and SGC galaxy samples, 
# and the "zmin" and "zmax" parameters to determine the redshift bins. 

# In addition, to study the effect of observational systematics described by weight functions, we calculate the weights of the galaxy sample in four cases: 
# 0) when all weights are considered;
# 1) when only systematic weights are considered;
# 2) when there are no systematic weights and only fiber collision and redshift failure are taken into account;
# 3) when not all observational weights are considered.
# 4) when the observational weights are considered but the FKP weights are ignored.

# When the observational weight function is varied, it is necessary to reconstruct random particles. 
# In particular, for NGC, CMASSLOWZ, LOWZE2, and LOWZE3 have different survey regions, so it is necessary to reconstruct random particles by considering the observational weights in each region.

# The parameter that controls these observational weights is "Weight".

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

#NS = "North"
#zbin = 1
#Weight = 0

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

# Function for loading galaxy random data #
def read_samples(fname, zmin, zmax):

    with fits.open(fname) as hdul0:
    	RA  = hdul0[1].data["RA"]
    	DEC = hdul0[1].data["DEC"]
    	Z   = hdul0[1].data["Z"]
    
    Z_bin   = (Z[Z>zmin])[Z[Z>zmin]<zmax]
    RA_bin  = (RA[Z>zmin])[Z[Z>zmin]<zmax]
    DEC_bin = (DEC[Z>zmin])[Z[Z>zmin]<zmax]
    del RA, DEC, Z
    X = np.array([RA_bin, DEC_bin, Z_bin]).T
    del RA_bin, DEC_bin, Z_bin
    X = list(map(tuple, X))
    X_set = set(X)
    del X 
    return X_set

def calc_set(TOT_set, CL_set, CLE2_set, CLE3_set):

#    CL_CLE2_set = CL_set&CLE2_set
#    CL_CLE3_set = CL_set&CLE3_set
#    CLE2_CLE3_set = CLE2_set&CLE3_set
#    CLTOT_set = CL_set&CLE2_set&CLE3_set
    
#    M23_set = CLE2_CLE3_set - CLTOT_set
#    LE2_set = CLE2_set - CL_CLE2_set - M23_set
#    LE3_set = CLE3_set - CL_CLE3_set - M23_set

    CL_set = CL_set&TOT_set
    CLE2_set = CLE2_set&TOT_set
    CLE3_set = CLE3_set&TOT_set

    LE2_set = CLE2_set - CL_set&CLE2_set
    LE3_set = CLE3_set - CL_set&CLE3_set
   
    s = set()
    if CL_set == s:
        print("CL_set = NULL")
    if LE2_set == s:
        print("LE2_set = NULL")
    if LE3_set == s:
        print("LE3_set = NULL")
#    if M23_set == s:
#        print("M23_set = NULL")
    
    TOT = np.array(list(map(list, TOT_set))).T
    CL = np.array(list(map(list, CL_set))).T
    LE2 = np.array(list(map(list, LE2_set))).T
    LE3 = np.array(list(map(list, LE3_set))).T

    del TOT_set, CL_set, CLE2_set, CLE3_set

    return TOT[0], TOT[1], TOT[2],\
            CL[0],  CL[1],  CL[2],\
           LE2[0], LE2[1], LE2[2],\
           LE3[0], LE3[1], LE3[2]

###########################

# loading weight functions from galaxy samples
def read_weights(fname, zmin, zmax):

    with fits.open(fname) as hdul:
    	data = hdul[1].data
    
    RA = data["RA"]
    DEC = data["DEC"]
    Z = data["Z"]
    weight_tot = data["WEIGHT_SYSTOT"] * ( data["WEIGHT_CP"] + data["WEIGHT_NOZ"] - 1.0 )
    weight_fkp = data["WEIGHT_FKP"]

    if Weight == 0:
        W = weight_tot * weight_fkp
    elif Weight == 1:
        W = data["WEIGHT_SYSTOT"] * weight_fkp
    elif Weight == 2:
        W = ( data["WEIGHT_CP"] + data["WEIGHT_NOZ"] - 1.0 ) * weight_fkp
    elif Weight == 3:
        W = weight_fkp
    elif Weight == 4:
        W = weight_tot
    elif Weight == 5:
        W = data["WEIGHT_SYSTOT"]
    else:
        print("ERROR")
        sys.exit()

    RA_bin = (RA[Z>zmin])[Z[Z>zmin]<zmax]
    DEC_bin = (DEC[Z>zmin])[Z[Z>zmin]<zmax]
    Z_bin = (Z[Z>zmin])[Z[Z>zmin]<zmax]
    W_bin = (W[Z>zmin])[Z[Z>zmin]<zmax]
    return RA_bin, DEC_bin, Z_bin, W_bin

def assign_weights(W_G, RA_G, RA):

    W = np.zeros(len(RA))
    for i in range(len(RA)):
        W[i] = W_G[RA_G == RA[i]]
    
    return W

def assign_random_ZW(Z_G, W_G, Ngal, Nran):

    intZ_CL = np.random.choice(Ngal, Nran)
    Z_R = Z_G[intZ_CL]
    W_R = W_G[intZ_CL]
    
    return Z_R, W_R


def calc_random(RANDOM, NS, zmin, zmax):

    TOT_ran_set  = read_samples("random%d_DR12v5_CMASSLOWZTOT_%s.fits" % (RANDOM, NS), zmin, zmax)
    CL_ran_set   = read_samples("random%d_DR12v5_CMASSLOWZ_%s.fits" % (RANDOM, NS), zmin, zmax)
    CLE2_ran_set = read_samples("random%d_DR12v5_CMASSLOWZE2_%s.fits" % (RANDOM, NS), zmin, zmax)
    CLE3_ran_set = read_samples("random%d_DR12v5_CMASSLOWZE3_%s.fits" % (RANDOM, NS), zmin, zmax)
    
    RA_TOT_ran, DEC_TOT_ran, Z_TOT_ran_temp,\
     RA_CL_ran,  DEC_CL_ran,  Z_CL_ran_temp,\
    RA_LE2_ran, DEC_LE2_ran, Z_LE2_ran_temp,\
    RA_LE3_ran, DEC_LE3_ran, Z_LE3_ran_temp = calc_set(TOT_ran_set, CL_ran_set, CLE2_ran_set, CLE3_ran_set)
    
    print("random")
    print(len(RA_CL_ran), len(RA_LE2_ran), len(RA_LE3_ran))
    print(len(RA_TOT_ran), len(RA_CL_ran)+len(RA_LE2_ran)+len(RA_LE3_ran))
    
    TOT_set  = read_samples("galaxy_DR12v5_CMASSLOWZTOT_%s.fits" % NS, zmin, zmax)
    CL_set   = read_samples("galaxy_DR12v5_CMASSLOWZ_%s.fits" % NS, zmin, zmax)
    CLE2_set = read_samples("galaxy_DR12v5_CMASSLOWZE2_%s.fits" % NS, zmin, zmax)
    CLE3_set = read_samples("galaxy_DR12v5_CMASSLOWZE3_%s.fits" % NS, zmin, zmax)
    
    RA_TOT, DEC_TOT, Z_TOT,\
     RA_CL,  DEC_CL,  Z_CL,\
    RA_LE2, DEC_LE2, Z_LE2,\
    RA_LE3, DEC_LE3, Z_LE3 = calc_set(TOT_set, CL_set, CLE2_set, CLE3_set)
    
    print("galaxy")
    print(len(RA_CL), len(RA_LE2), len(RA_LE3))
    print(len(RA_TOT), len(RA_CL)+len(RA_LE2)+len(RA_LE3))
    
    RA_G, DEC_G, Z_G, W_G = read_weights("galaxy_DR12v5_CMASSLOWZTOT_%s.fits" % NS, zmin, zmax)
    
    W_CL  = assign_weights(W_G, RA_G, RA_CL)
    W_LE2 = assign_weights(W_G, RA_G, RA_LE2)
    W_LE3 = assign_weights(W_G, RA_G, RA_LE3)
    
    Z_CL_ran, W_CL_ran = assign_random_ZW(Z_CL, W_CL, len(RA_CL), len(RA_CL_ran))
    Z_LE2_ran, W_LE2_ran = assign_random_ZW(Z_LE2, W_LE2, len(RA_LE2), len(RA_LE2_ran))
    Z_LE3_ran, W_LE3_ran = assign_random_ZW(Z_LE3, W_LE3, len(RA_LE3), len(RA_LE3_ran))
    
    RA_ran = np.hstack([RA_CL_ran, RA_LE2_ran, RA_LE3_ran])
    DEC_ran = np.hstack([DEC_CL_ran, DEC_LE2_ran, DEC_LE3_ran])
    Z_ran = np.hstack([Z_CL_ran, Z_LE2_ran, Z_LE3_ran])
    W_ran = np.hstack([W_CL_ran, W_LE2_ran, W_LE3_ran])

    return RA_ran, DEC_ran, Z_ran, W_ran

def calc_random_South(RANDOM, NS, zmin, zmax):
    
    fname = "random%d_DR12v5_CMASSLOWZTOT_%s.fits" % (RANDOM, NS)
    with fits.open(fname) as hdul0:
    	RA  = hdul0[1].data["RA"]
    	DEC = hdul0[1].data["DEC"]
    	Z   = hdul0[1].data["Z"]
    
    RA_ran  = (RA[Z>zmin])[Z[Z>zmin]<zmax]
    DEC_ran = (DEC[Z>zmin])[Z[Z>zmin]<zmax]

    RA_G, DEC_G, Z_G, W_G = read_weights("galaxy_DR12v5_CMASSLOWZTOT_%s.fits" % NS, zmin, zmax)
    Z_ran, W_ran = assign_random_ZW(Z_G, W_G, len(RA_G), len(RA_ran))

    return RA_ran, DEC_ran, Z_ran, W_ran

if NS == "North":
    RA_ran0, DEC_ran0, Z_ran0, W_ran0 = calc_random(0, NS, zmin, zmax)
    RA_ran1, DEC_ran1, Z_ran1, W_ran1 = calc_random(1, NS, zmin, zmax)
elif NS == "South":
    RA_ran0, DEC_ran0, Z_ran0, W_ran0 = calc_random_South(0, NS, zmin, zmax)
    RA_ran1, DEC_ran1, Z_ran1, W_ran1 = calc_random_South(1, NS, zmin, zmax)

RA_bin = np.hstack([RA_ran0, RA_ran1])
DEC_bin = np.hstack([DEC_ran0, DEC_ran1])
Z_bin = np.hstack([Z_ran0, Z_ran1])
W_bin = np.hstack([W_ran0, W_ran1])

print(len(RA_ran0), len(RA_ran1), len(RA_bin))

#RA_ran = np.hstack([RA_CL_ran, RA_LE2_ran, RA_LE3_ran])
#DEC_ran = np.hstack([DEC_CL_ran, DEC_LE2_ran, DEC_LE3_ran])
#Z_ran = np.hstack([Z_CL_ran, Z_LE2_ran, Z_LE3_ran])
#W_ran = np.hstack([W_CL_ran, W_LE2_ran, W_LE3_ran])

# split into three redshift bins #
#Z_bin   = (Z[Z>zmin])[Z[Z>zmin]<zmax]
#RA_bin  = (RA[Z>zmin])[Z[Z>zmin]<zmax]
#DEC_bin = (DEC[Z>zmin])[Z[Z>zmin]<zmax]
#W_bin   = (W[Z>zmin])[Z[Z>zmin]<zmax]

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
if Weight == 0:
    output_fname = "random_DR12v5_CMASSLOWZTOT_%s_ZBIN%d.dat" % (NS,zbin)
elif Weight == 1:
    output_fname = "random_DR12v5_CMASSLOWZTOT_%s_ZBIN%d_Weight1_OnlySys.dat" % (NS,zbin)
elif Weight == 2:
    output_fname = "random_DR12v5_CMASSLOWZTOT_%s_ZBIN%d_Weight2_NoSys.dat" % (NS,zbin)
elif Weight == 3:
    output_fname = "random_DR12v5_CMASSLOWZTOT_%s_ZBIN%d_Weight3_NoWeight.dat" % (NS,zbin)
elif Weight == 4:
    output_fname = "random_DR12v5_CMASSLOWZTOT_%s_ZBIN%d_Weight4_NoFKP.dat" % (NS,zbin)
elif Weight == 5:
    output_fname = "random_DR12v5_CMASSLOWZTOT_%s_ZBIN%d_Weight5_OnlySys_NoFKP.dat" % (NS,zbin)

output_dir = "galaxy_DR12v5_CMASSLOWZTOT_Weight"
try:
    os.mkdir(output_dir)
except:
    print(output_dir, " already exists ")
 
np.savetxt("%s/%s" % (output_dir, output_fname), SAVE, fmt="%.10e")

