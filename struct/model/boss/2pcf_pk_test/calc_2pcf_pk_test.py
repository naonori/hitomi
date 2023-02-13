########################
# Purpose of this code #
########################
# This code explains how to calculate the power spectrum and 2PCF using "$WORK/cosmo/hitomi_theory". 
# The nonlinearity is taken into account only for the BAO damping effect. 
# It is also possible to calculate the BAO damping effect after reconstruction.
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

# Redshift #
zbin = 1 # 1,3
if zbin == 1:
    redshift = 0.38 
elif zbin == 3:
    redshift = 0.61

# r-bin #
rbin = np.linspace(30.0, 150.0, 25)
# k-bin #
kbin = np.linspace(0.01,0.2,20)

##################################################################################################################################
# BOSS fiducial parameters #
# The values of the fiducial parameters should be exactly the same as those used to calculate the distance to the BOSS galaxies.
# (see $WORK/data/boss/fits2dat_galaxy.py)
##################################################################################################################################
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
print("alpha_perp =", alpha_perp)
print("alpha_parallel =", alpha_parallel)
print("Growth factor =", Dz)
print("Growth rate =", fz)
print("sigma8(z) =", sigma8z)

b1 = 2.0
fsigma8 = fz * sigma8z
b1sigma8 = b1 * sigma8z

print("fsigma8 =",  fsigma8)
print("b1sigma8 =", b1sigma8)
print()

###############################################
# Calculate the power spectrum and two-point correlation function 
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
 
# make an output directory #
output_dir = "results_test"
try:
    os.mkdir(output_dir)
except:
    pass

### ##################################
### # From here, we calculate the multipole components of the power specctrum.
### #
### # First of all, the "name" parameter specifies the type of power spectra to calculate.
### # All possible "name" parameters are the following: 
### # "Tree" means the normal linear (tree-level) power spectrum;
### # "Tree_NoWiggle" means the no-wiggle linear power spectrum;
### # "Tree_BAO" means the linear power spectrum with the damping effect of the BAO signal, proposed by Eisenstein 2007.
### # "Tree_BAO_b1_b1", "Tree_BAO_b1_f", "Tree_BAO_f_f" are decompositions of "Tree_BAO" into parts that do not depend on the "fsigma8" and "b1sigma8" parameters.
### #
### # The "kbin" parameter specifies the binning of the power spectrum to be calculated.
### # For example, kbin = np.linspace(0.01, 0.2, 20).
### #
### # The "ELL" parameter specifies the multipole component to be calculated.
### # It is only allowed for the three cases: ELL = 0, 2, 4.
### #
### 
### #####
### # The simplest example is to compute the monopole component of the linear power spectrum.
### name = "Tree"
### ELL = 0
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL)
### # Then, a dictionary object is returned.
### # Print the keys contained in "pk_dict".
### print(pk_dict.keys()) #dict_keys(['kbin', 'P', 'kbin_fft', 'P_fft', 'ELL', 'flag_2pcf', 'flag_BAO', 'flag_recon', 'sigma2_perp', 'sigma2_para'])
### # To get the resulting power spectrum corresponding to the input "kbin", use "kbin" and "P".
### k_save = pk_dict["kbin"]
### pk_save = pk_dict["P"]
### # Save the results.
### X = np.array([k_save, pk_save]).T
### np.savetxt("%s/pk%s_%s.dat" % (output_dir, ELL, name), X, fmt="%.5f \t %.7e" )
### ##
### 
### ######
### # To calculate the two-point correlation function (2PCF), add "flag_2pcf=True".
### # By default, "flag_2pcf=False".
### name = "Tree"
### ELL = 0
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True)
### # In this case, apart from the input "kbin", a power spectrum in the range 0.0003<k<10 is computed to compute the Hankel transform needed to obtain the 2PCF. 
### # The result is included in the keys "kbin_fft" and "P_fft".
### k_save = pk_dict["kbin_fft"]
### pk_save = pk_dict["P_fft"]
### # Save the results.
### X = np.array([k_save, pk_save]).T
### np.savetxt("%s/pk%s_%s_fft.dat" % (output_dir, ELL, name), X, fmt="%.5f \t %.7e" )
### 
### # The "calc_P_to_2PCF" function, which calculates the 2PCF, returns a directory object.
### # Input "rbin", e.g., rbin = np.linspace(30, 150, 25)
### xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
### # Print the keys contained in "xi_dict".
### print(xi_dict.keys()) # dict_keys(['rbin', '2pcf', 'rbin_fft', '2pcf_fft', 'ELL', 'flag_2pcf', 'flag_BAO', 'flag_recon', 'N_fftlog'])
### # To get the resulting 2PCF corresponding to the input "rbin", use "rbin" and "2pcf".
### r_save = xi_dict["rbin"]
### xi_save = xi_dict["2pcf"]
### # Save the results.
### X = np.array([r_save, xi_save]).T
### np.savetxt("%s/xi%s_%s.dat" % (output_dir, ELL, name), X, fmt="%.5f \t %.7e" )
### 
### # The data needed to compute the inverse Hankel transform, which calculates the power spectrum from the 2PCF, is included in "rbin_fft" and "2pcf_fft".
### # The range of "rbin_fft" is 0.1 < r < 3300 [Mpc/h].
### r_save = xi_dict["rbin_fft"]
### xi_save = xi_dict["2pcf_fft"]
### # Save the results.
### X = np.array([r_save, xi_save]).T
### np.savetxt("%s/xi%s_%s_fft.dat" % (output_dir, ELL, name), X, fmt="%.5f \t %.7e" )
### 
### 
### ######
### # The "calc_2PCF_to_P" function computes the power spectrum from the 2PCF. 
### # This calculation is useful when performing window function corrections to the power spectrum.
### pk_dict_from_2pcf = P.calc_2PCF_to_P(xi_dict, kbin=kbin)
### print(pk_dict_from_2pcf.keys())
### # Save the results.
### k_save = pk_dict_from_2pcf["kbin"]
### pk_save = pk_dict_from_2pcf["P"]
### X = np.array([k_save, pk_save]).T
### np.savetxt("%s/pk%s_%s_from_2pcf.dat" % (output_dir, ELL, name), X, fmt="%.5f \t %.7e" )
### # Make sure that the results "pk0_Tree_from_2pcf.dat" are in good agreement with the original "pk0_Tree.dat".
### ######
### 
### #######
### ## Calculate the quadrupole and hexadecapole of both the linear power spectrum and the 2PCF and save them.
### ELL=2
### name="Tree"
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True)
### xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
### pk_dict_from_2pcf = P.calc_2PCF_to_P(xi_dict, kbin=kbin)
### np.savetxt("%s/pk%s_%s.dat" % (output_dir, ELL, name), np.array([pk_dict["kbin"], pk_dict["P"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/pk%s_%s_fft.dat" % (output_dir, ELL, name), np.array([pk_dict["kbin_fft"], pk_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s.dat" % (output_dir, ELL, name), np.array([xi_dict["rbin"], xi_dict["2pcf"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_fft.dat" % (output_dir, ELL, name), np.array([xi_dict["rbin_fft"], xi_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/pk%s_%s_from_2pcf.dat" % (output_dir, ELL, name), np.array([pk_dict_from_2pcf["kbin"], pk_dict_from_2pcf["P"]]).T, fmt="%.5f \t %.7e" )
### 
### ELL=4
### name="Tree"
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True)
### xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
### pk_dict_from_2pcf = P.calc_2PCF_to_P(xi_dict, kbin=kbin)
### np.savetxt("%s/pk%s_%s.dat" % (output_dir, ELL, name), np.array([pk_dict["kbin"], pk_dict["P"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/pk%s_%s_fft.dat" % (output_dir, ELL, name), np.array([pk_dict["kbin_fft"], pk_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s.dat" % (output_dir, ELL, name), np.array([xi_dict["rbin"], xi_dict["2pcf"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_fft.dat" % (output_dir, ELL, name), np.array([xi_dict["rbin_fft"], xi_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/pk%s_%s_from_2pcf.dat" % (output_dir, ELL, name), np.array([pk_dict_from_2pcf["kbin"], pk_dict_from_2pcf["P"]]).T, fmt="%.5f \t %.7e" )
### 
### #######
### # Calculate the no-wiggle power spectrum multipoles 
### ELL=0
### name="Tree_NoWiggle"
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True)
### xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
### np.savetxt("%s/pk%s_%s_fft.dat" % (output_dir, ELL, name), np.array([pk_dict["kbin_fft"], pk_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_fft.dat" % (output_dir, ELL, name), np.array([xi_dict["rbin_fft"], xi_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### ELL=2
### name="Tree_NoWiggle"
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True)
### xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
### np.savetxt("%s/pk%s_%s_fft.dat" % (output_dir, ELL, name), np.array([pk_dict["kbin_fft"], pk_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_fft.dat" % (output_dir, ELL, name), np.array([xi_dict["rbin_fft"], xi_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### ELL=4
### name="Tree_NoWiggle"
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True)
### xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
### np.savetxt("%s/pk%s_%s_fft.dat" % (output_dir, ELL, name), np.array([pk_dict["kbin_fft"], pk_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_fft.dat" % (output_dir, ELL, name), np.array([xi_dict["rbin_fft"], xi_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### #######
### # Calculate the model describing the damping of the BAO signal given by Eisenstein 2007.
### # This model is given by
### # P = exp(-k**2 ( (1-mu**2)*sigma2_perp + mu**2*sigma2_para) ) * ( P_lin - P_nw ) + P_nw
### # where P_lin and P_nw are the wiggle and no-wiggle linear power spectrum including the linear Kaiser effect and the linear bias effect,
### # and mu is the cosine of the angle between wavevector and the line-of-sight.
### # The BAO damping is characterized by the two parameters: "sigma2_perp" and "sigma2_para".
### #
### # By default, the "sigma2_perp" and "sigma2_para" parameters are calculated by the linear Lagrangian perturbation theory with given cosmological parameters:
### # 
### # When calculating "Tree_BAO", be sure to add "flag_BAO=True".
### 
### ELL=0
### name="Tree_BAO"
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True)
### # Print the computed values of "sigma2_perp" and "sigma2_para".
### sigma2_perp = pk_dict["sigma2_perp"]
### sigma2_para = pk_dict["sigma2_para"]
### print("sigma2_perp = ", sigma2_perp)
### print("sigma2_para = ", sigma2_para)
### # Alternatively, you can pre-specify the sigma8 and fz parameters needed to compute sigma2_perp and sigma2_para.
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True, sigma8_fid=sigma8z, fz_fid=fz)
### sigma2_perp = pk_dict["sigma2_perp"]
### sigma2_para = pk_dict["sigma2_para"]
### print("sigma2_perp = ", sigma2_perp)
### print("sigma2_para = ", sigma2_para)
### # You can also use the values of "sigma2_perp" and "sigma2_para" themselves as input parameters.
### # In this case, be sure to add "flag_damping=True".
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True, flag_damping=True, sigma2_perp=sigma2_perp, sigma2_para=sigma2_para)
### print("sigma2_perp = ", pk_dict["sigma2_perp"])
### print("sigma2_para = ", pk_dict["sigma2_para"])
### 
### ######
### # Calculate the multipole components of both the power spectrum and the 2PCF, using Eisenstein 2007's template model.
### ELL=0
### name="Tree_BAO"
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True)
### xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
### np.savetxt("%s/pk%s_%s_fft.dat" % (output_dir, ELL, name), np.array([pk_dict["kbin_fft"], pk_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_fft.dat" % (output_dir, ELL, name), np.array([xi_dict["rbin_fft"], xi_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### ELL=2
### name="Tree_BAO"
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True)
### xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
### np.savetxt("%s/pk%s_%s_fft.dat" % (output_dir, ELL, name), np.array([pk_dict["kbin_fft"], pk_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_fft.dat" % (output_dir, ELL, name), np.array([xi_dict["rbin_fft"], xi_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### ELL=4
### name="Tree_BAO"
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True)
### xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
### np.savetxt("%s/pk%s_%s_fft.dat" % (output_dir, ELL, name), np.array([pk_dict["kbin_fft"], pk_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_fft.dat" % (output_dir, ELL, name), np.array([xi_dict["rbin_fft"], xi_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### ######
### # Compute the reconstructed model.
### # Here, we assume that only the values of "sigma2_perp" and "sigma2_para" are changed by reconstruction. 
### # The reconstructed "sigma2_perp" and "sigma2_para" are computed by linear Lagrangian perturbation theory.
### # Add "flag_recon=True".
### # The input parameters are the linear bias parameter (b1_fid) and the Gaussian smoothing parameter (R). 
### # These values must be exactly the same as those used to reconstruct the data.
### #
### 
### for R in [5, 10, 15, 20]:
###     
###     ELL=0
###     name="Tree_BAO"
###     pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True, flag_recon = True, one_over_b1_fid= 1.0/2.0, R= R)
###     xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
###     np.savetxt("%s/pk%s_%s_recon_R%02d_fft.dat" % (output_dir, ELL, name, R), np.array([pk_dict["kbin_fft"], pk_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
###     np.savetxt("%s/xi%s_%s_recon_R%02d_fft.dat" % (output_dir, ELL, name, R), np.array([xi_dict["rbin_fft"], xi_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
###     print("reconstructed sigma2_perp =", pk_dict["sigma2_perp"])
###     print("reconstructed sigma2_para =", pk_dict["sigma2_para"])
###     
###     ELL=2
###     name="Tree_BAO"
###     pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True, flag_recon = True, one_over_b1_fid= 1.0/2.0, R= R)
###     xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
###     np.savetxt("%s/pk%s_%s_recon_R%02d_fft.dat" % (output_dir, ELL, name, R), np.array([pk_dict["kbin_fft"], pk_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
###     np.savetxt("%s/xi%s_%s_recon_R%02d_fft.dat" % (output_dir, ELL, name, R), np.array([xi_dict["rbin_fft"], xi_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
###     
###     ELL=4
###     name="Tree_BAO"
###     pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True, flag_recon = True, one_over_b1_fid= 1.0/2.0, R= R)
###     xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
###     np.savetxt("%s/pk%s_%s_recon_R%02d_fft.dat" % (output_dir, ELL, name, R), np.array([pk_dict["kbin_fft"], pk_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
###     np.savetxt("%s/xi%s_%s_recon_R%02d_fft.dat" % (output_dir, ELL, name, R), np.array([xi_dict["rbin_fft"], xi_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### ######
### # Decompose the power spectrum (2PCF) computed by "Tree_BAO" into terms independent of "fsigma8" and "b1sigma8", and store them: "Tree_BAO_b1_b1", "Tree_BAO_b1_f", "Tree_BAO_f_f".
### # This way, if we want to use the power spectrum (2PCF) to estimate only "fsigma8" and "b1sigma8", in other words, if we want to ignore the AP effect, we do not need to recompute the power spectrum in "montepython".
### # Here, "name" is fixed to "Tree_BAO_Template", and a new parameter "param_name" is added.
### 
### ELL=0
### name="Tree_BAO_Template"
### param_name = "b1_b1"
### pk_b1_b1_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True, param_name = param_name)
### xi_b1_b1_dict = P.calc_P_to_2PCF(pk_b1_b1_dict, rbin=rbin)
### np.savetxt("%s/pk%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([pk_b1_b1_dict["kbin_fft"], (b1sigma8)**2 * pk_b1_b1_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([xi_b1_b1_dict["rbin_fft"], (b1sigma8)**2 * xi_b1_b1_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### param_name = "b1_f"
### pk_b1_f_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True, param_name = param_name)
### xi_b1_f_dict = P.calc_P_to_2PCF(pk_b1_f_dict, rbin=rbin)
### np.savetxt("%s/pk%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([pk_b1_f_dict["kbin_fft"], (b1sigma8) * (fsigma8) * pk_b1_f_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([xi_b1_f_dict["rbin_fft"], (b1sigma8) * (fsigma8) * xi_b1_f_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### param_name = "f_f"
### pk_f_f_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True, param_name = param_name)
### xi_f_f_dict = P.calc_P_to_2PCF(pk_f_f_dict, rbin=rbin)
### np.savetxt("%s/pk%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([pk_f_f_dict["kbin_fft"], (fsigma8)**2 * pk_f_f_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([xi_f_f_dict["rbin_fft"], (fsigma8)**2 * xi_f_f_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### # The original "Tree_BAO" power spectum and 2PCF are given by
### pk = (b1sigma8)**2 * pk_b1_b1_dict["P"]\
###    + (b1sigma8) * (fsigma8) * pk_b1_f_dict["P"]\
###    + (fsigma8)**2 * pk_f_f_dict["P"]
### 
### xi = (b1sigma8)**2 * xi_b1_b1_dict["2pcf"]\
###    + (b1sigma8) * (fsigma8) * xi_b1_f_dict["2pcf"]\
###    + (fsigma8)**2 * xi_f_f_dict["2pcf"]
### 
### # Compare the above "pk" and "xi" with those computed by "Tree_BAO" to make sure they agree well.
### name="Tree_BAO"
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True)
### xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
### print("sum(diff.) =", np.sum(pk_dict["P"] - pk))
### print("sum(diff.) =", np.sum(xi_dict["2pcf"] - xi))
### 
### ### #################
### ### # Perform the same tests for the quadrupole and hexadecapole as above.
### 
### ELL = 2
### name="Tree_BAO_Template"
### param_name = "b1_b1"
### pk_b1_b1_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True, param_name = param_name)
### xi_b1_b1_dict = P.calc_P_to_2PCF(pk_b1_b1_dict, rbin=rbin)
### np.savetxt("%s/pk%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([pk_b1_b1_dict["kbin_fft"], (b1sigma8)**2 * pk_b1_b1_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([xi_b1_b1_dict["rbin_fft"], (b1sigma8)**2 * xi_b1_b1_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### param_name = "b1_f"
### pk_b1_f_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True, param_name = param_name)
### xi_b1_f_dict = P.calc_P_to_2PCF(pk_b1_f_dict, rbin=rbin)
### np.savetxt("%s/pk%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([pk_b1_f_dict["kbin_fft"], (b1sigma8) * (fsigma8) * pk_b1_f_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([xi_b1_f_dict["rbin_fft"], (b1sigma8) * (fsigma8) * xi_b1_f_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### param_name = "f_f"
### pk_f_f_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True, param_name = param_name)
### xi_f_f_dict = P.calc_P_to_2PCF(pk_f_f_dict, rbin=rbin)
### np.savetxt("%s/pk%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([pk_f_f_dict["kbin_fft"], (fsigma8)**2 * pk_f_f_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([xi_f_f_dict["rbin_fft"], (fsigma8)**2 * xi_f_f_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### # The original "Tree_BAO" power spectum and 2PCF are given by
### pk = (b1sigma8)**2 * pk_b1_b1_dict["P"]\
###    + (b1sigma8) * (fsigma8) * pk_b1_f_dict["P"]\
###    + (fsigma8)**2 * pk_f_f_dict["P"]
### 
### xi = (b1sigma8)**2 * xi_b1_b1_dict["2pcf"]\
###    + (b1sigma8) * (fsigma8) * xi_b1_f_dict["2pcf"]\
###    + (fsigma8)**2 * xi_f_f_dict["2pcf"]
### 
### # Compare the above "pk" and "xi" with those computed by "Tree_BAO" to make sure they agree well.
### name="Tree_BAO"
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True)
### xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
### print("sum(diff.) =", np.sum(pk_dict["P"] - pk))
### print("sum(diff.) =", np.sum(xi_dict["2pcf"] - xi))
### 
### #############################
### #############################
### 
### ELL = 4
### name="Tree_BAO_Template"
### param_name = "b1_b1"
### pk_b1_b1_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True, param_name = param_name)
### xi_b1_b1_dict = P.calc_P_to_2PCF(pk_b1_b1_dict, rbin=rbin)
### np.savetxt("%s/pk%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([pk_b1_b1_dict["kbin_fft"], (b1sigma8)**2 * pk_b1_b1_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([xi_b1_b1_dict["rbin_fft"], (b1sigma8)**2 * xi_b1_b1_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### param_name = "b1_f"
### pk_b1_f_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True, param_name = param_name)
### xi_b1_f_dict = P.calc_P_to_2PCF(pk_b1_f_dict, rbin=rbin)
### np.savetxt("%s/pk%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([pk_b1_f_dict["kbin_fft"], (b1sigma8) * (fsigma8) * pk_b1_f_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([xi_b1_f_dict["rbin_fft"], (b1sigma8) * (fsigma8) * xi_b1_f_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### param_name = "f_f"
### pk_f_f_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True, param_name = param_name)
### xi_f_f_dict = P.calc_P_to_2PCF(pk_f_f_dict, rbin=rbin)
### np.savetxt("%s/pk%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([pk_f_f_dict["kbin_fft"], (fsigma8)**2 * pk_f_f_dict["P_fft"]]).T, fmt="%.5f \t %.7e" )
### np.savetxt("%s/xi%s_%s_%s_fft.dat" % (output_dir, ELL, name, param_name),\
###         np.array([xi_f_f_dict["rbin_fft"], (fsigma8)**2 * xi_f_f_dict["2pcf_fft"]]).T, fmt="%.5f \t %.7e" )
### 
### # The original "Tree_BAO" power spectum and 2PCF are given by
### pk = (b1sigma8)**2 * pk_b1_b1_dict["P"]\
###    + (b1sigma8) * (fsigma8) * pk_b1_f_dict["P"]\
###    + (fsigma8)**2 * pk_f_f_dict["P"]
### 
### xi = (b1sigma8)**2 * xi_b1_b1_dict["2pcf"]\
###    + (b1sigma8) * (fsigma8) * xi_b1_f_dict["2pcf"]\
###    + (fsigma8)**2 * xi_f_f_dict["2pcf"]
### 
### # Compare the above "pk" and "xi" with those computed by "Tree_BAO" to make sure they agree well.
### name="Tree_BAO"
### pk_dict = P.calc_P(name=name, kbin = kbin, ELL=ELL, flag_2pcf=True, flag_BAO=True)
### xi_dict = P.calc_P_to_2PCF(pk_dict, rbin=rbin)
### print("sum(diff.) =", np.sum(pk_dict["P"] - pk))
### print("sum(diff.) =", np.sum(xi_dict["2pcf"] - xi))
### 
