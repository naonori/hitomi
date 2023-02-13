########################
# Purpose of this code #
########################
# This code computes the parameter-decomposed 2PCFs at z=0.38 and z=0.61. In addition, it also computes the reconstructed 2PCFs.
######################################################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np

from classy import Class

# $WORK/cosmo/hitomi_theory/ #
import fiducial
import initial
import powerspec

data = "mock" # "galaxy" or "mock"

############
# Redshift #
############
# Calculates the theoretical models for the first and third BOSS galaxy samples when divided into three redshift ranges:
# 0.2 < z < 0.5, 0.4 < z < 0.6, and 0.5 < z < 0.75.
# (see $WORK/data/boss/fits2dat_galaxy.py)
##############################################
zbin = 3 # 1,3
if zbin == 1:
    redshift = 0.38 
elif zbin == 3:
    redshift = 0.61

#########
# r-bin #
#########
# Calculate the 2PCF for the same r-bins as in the case of the 2PCF measurement.
# (see /mwork0/sugiymnn/WORK/measurement/boss/2pcf/PARAMS)
##########################################################
rbin = np.linspace(30.0, 150.0, 25)

############################
# make an output directory #
############################
OUTPUT = "model_%s_zbin%d" % (data, zbin)
try:
    os.mkdir(OUTPUT)
except:
    pass

if data == "galaxy":
    ################################################################################
    # BOSS fiducial parameters #
    # The values of the fiducial parameters should be exactly the same as those used to calculate the distance to the BOSS galaxies.
    # (see $WORK/data/boss/fits2dat_galaxy.py)
    ################################################################################
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
    
        'output': 'tCl, mPk',
        'P_k_max_h/Mpc': 50.0,
        'z_max_pk': 3.0,
    }

elif data == "mock":

    ################################################################################
    # Patchy mock fiducial parameters #
    # The values of the fiducial parameters should be exactly the same as those used to calculate the distance to the Patchy mock galaxies.
    # (see $WORK/data/boss/fits2dat_galaxy.py)
    ################################################################################

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
    
        'output': 'tCl, mPk',
        'P_k_max_h/Mpc': 50.0,
        'z_max_pk': 3.0,
    
    }

else:
    print("ERROR")
    exit()

##############################################################################################
# To calculate the AP parameters, we use the fiducial cosmology parameters to calculate 
# the angular diameter distance (Da_fid) and the Hubble parameter (H_fid) at a given redshift.
##############################################################################################
fiducial_cosmo = fiducial.FiducialHubbleParameterAndAngularDiameterDistance(redshift, params_cosmo)
fiducial_cosmo.calcHubbleParameterAndAngularDiameterDistance()
Da_fid = fiducial_cosmo.getFiducialAngularDiameterDistance()
H_fid = fiducial_cosmo.getFiducialHubbleParameter()

###############
# Call Class() 
###############
cosmo = Class()
cosmo.set(params_cosmo)
cosmo.compute()

###############
# Compute the required functions:
# the linear matter power spectrum, the no-wiggle linear matter power spectrum,
# the growth factor, the growth rate, sigma8 at a given redshift, 
# and the Alcock-Paczyn'ski parameters.
###############
initial_cosmo = initial.InputPowerSpectrum(redshift, cosmo)
initial_cosmo.calcMatterPowerSpectrum()
initial_cosmo.calcNoWiggleMatterPowerSpectrum()
initial_cosmo.calcPrimordialPowerSpectrum(ln10A_s10 = params_cosmo["ln10^{10}A_s"], n_s = params_cosmo["n_s"])

k_in, pk_in = initial_cosmo.getMatterPowerSpectrum()
k_in, pk_nw_in = initial_cosmo.getNoWiggleMatterPowerSpectrum()
k_in, mk_pri_in = initial_cosmo.getTransferFunctionM()
Dz = initial_cosmo.getGrowthFactor()
fz = initial_cosmo.getGrowthRate()
alpha_perp = initial_cosmo.getAlphaPerp(Da_fid)
alpha_parallel = initial_cosmo.getAlphaParallel(H_fid)
sigma8z_norm = initial_cosmo.getSigma8ForNormalization()
sigma8z = sigma8z_norm # If the value of sigma8(z=0) is already known, such as in the case of simulated data, use sigma8z = Dz * sigma8(z=0).

print()
print("alpha_perp = ", alpha_perp)
print("alpha_parallel = ", alpha_parallel)
print("Growth factor = ", Dz)
print("Growth rate = ", fz)
print("sigma8(z) = ", sigma8z)

b1 = 2.0
fsigma8 = fz * sigma8z
b1sigma8 = b1 * sigma8z

print("fsigma8 = ",  fsigma8)
print("b1sigma8 = ", b1sigma8)
print()
###############################################
# Calculate the two-point correlation function 
###############################################
# When ignoring the RSD effect, set fz = 0
# For dark matter, set b1 = 0.
params = {
    'alpha_perp': alpha_perp,
    'alpha_parallel': alpha_parallel,
    'sigma8': sigma8z,
    'fz': fsigma8/sigma8z,
    'b1': b1sigma8/sigma8z,
}

# Call "ClassPowerSpectrum()"
P = powerspec.ClassPowerSpectrum()
# Input linear power spectra
P.set_input_pk(k_in, pk_in)
P.set_input_pk_nw(k_in, pk_nw_in)
# Set "sigma8" used to nomalize the input linear power spectrum.
P.set_normalization(sigma8z_norm)
# Set input parameters 
P.set_params(params)

######
# Calculate the parameter-decomposed 2PCF multipoles before reconstruction. #
######

name = "Tree_BAO_Template"
for ELL in [0,2,4]:
    for param_name in ["b1_b1", "b1_f", "f_f"]:

        pk_dict = P.calc_P(name=name, ELL=ELL, flag_2pcf=True, flag_BAO=True, param_name = param_name)
        xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
        np.savetxt("%s/xi%s_%s_%s.dat" % (OUTPUT, ELL, name, param_name), np.array([xi_dict["rbin"], xi_dict["2pcf"]]).T, fmt="%.5f \t %.7e" )

######
# Calculate the parameter-decomposed 2PCF multipoles after reconstruction. #
######
# Calculate the reconstructed 2PCF for the same b1_fid and R as in the case of measuring the reconstructed 2PCF.
# (see /mwork0/sugiymnn/WORK/measurement/boss/2pcf/PARAMS)
######
b1_fid = 2.0
R = 15.0
name = "Tree_BAO_Template"
for ELL in [0,2,4]:
    for param_name in ["b1_b1", "b1_f", "f_f"]:
        
        pk_dict = P.calc_P(name=name, ELL=ELL, flag_2pcf=True, flag_BAO=True, flag_recon=True, one_over_b1_fid = 1.0 / b1_fid, R = R, param_name = param_name)
        xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
        np.savetxt("%s/xi%s_%s_%s_recon.dat" % (OUTPUT, ELL, name, param_name), np.array([xi_dict["rbin"], xi_dict["2pcf"]]).T, fmt="%.5f \t %.7e" )

