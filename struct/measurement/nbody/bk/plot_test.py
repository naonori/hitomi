########################
# Purpose of this code #
########################
#
# This code aims to plot the measured bispectra.
#
######################################################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate

params = {
#    'font.family': 'Times New Roman',
#    'text.usetex': 'True',
    'font.size': 11.0,
}

plt.rcParams.update(params)
#fig = plt.figure(figsize=(210.0/25.4, 264.0/25.4/2.3))
#fig = plt.figure(figsize=(1.8*210.0/25.4, 1.8*210.0/25.4/3.0))
fig = plt.figure(figsize=(7.0,7.0))
ax = []

##--------------------------
total_h = 0.94
total_w = 0.99

dw = 0.07
dh = 0.00
left = 0.15
bottom = 0.11
height = (total_h - (bottom + dh))
width  = (total_w - (left + dw))

#left bottom, width, height
rect1 = [left, bottom, width, height]
ax.append(fig.add_axes(rect1))

### ################
### ### Dark Matter
### ################
### for k1_i, k1 in enumerate(np.linspace(0.02, 0.2, 10)):
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
###     
###     ell1 = 0
###     ell2 = 0
###     ELL = 0
###     zbin = 3
###     ax[0].axvline(k1, color="black", ls="--")
###         
###     k, bk = np.loadtxt("results_test_Gadget_zbin%d_RSDTrue/0001/bk%d%d%d_%02d.dat" % (zbin, ell1, ell2, ELL, k1_i), usecols=(1,2), unpack=True)
### 
###     ax[0].plot(k, k*bk, "o-", label="bk%d%d%d w/ RSDs" % (ell1, ell2, ELL))
###     
###     ax[0].set_xlim(0.01, 0.21)
###     ax[0].set_xticks(np.linspace(0.02, 0.2, 10))
###     ax[0].set_xlabel("k2 [h/Mpc]")
###     ax[0].set_ylabel(r"$k2\, $ bispectrum(k1, k2)")
###     ax[0].legend(loc=0)
###     ax[0].set_title("k1 = %0.2f [h/Mpc]" % k1)
###     plt.savefig("figure/bk%d%d%d_%02d_matter.png" % (ell1, ell2, ELL, k1_i))
###     os.system("display figure/bk%d%d%d_%02d_matter.png" % (ell1, ell2, ELL, k1_i))
###
### ################
### ### Halos
### ################
### for k1_i, k1 in enumerate(np.linspace(0.02, 0.2, 10)):
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
###     
###     ell1 = 0
###     ell2 = 0
###     ELL = 0
###     zbin = 3
###     ax[0].axvline(k1, color="black", ls="--")
###         
###     k, bk = np.loadtxt("results_test_Rockstar_Mmin13.0_Mmax13.5_zbin%d_RSDTrue/0001/bk%d%d%d_%02d.dat" % (zbin, ell1, ell2, ELL, k1_i), usecols=(1,2), unpack=True)
### 
###     ax[0].plot(k, k*bk, "o-", label="bk%d%d%d: halo w/ RSDs" % (ell1, ell2, ELL))
###     
###     ax[0].set_xlim(0.01, 0.21)
###     ax[0].set_xticks(np.linspace(0.02, 0.2, 10))
###     ax[0].set_xlabel("k2 [h/Mpc]")
###     ax[0].set_ylabel(r"$k2\, $ bispectrum(k1, k2)")
###     ax[0].legend(loc=0)
###     ax[0].set_title("k1 = %0.2f [h/Mpc]" % k1)
###     plt.savefig("figure/bk%d%d%d_%02d_halo.png" % (ell1, ell2, ELL, k1_i))
###     os.system("display figure/bk%d%d%d_%02d_halo.png" % (ell1, ell2, ELL, k1_i))
###
### ################
### ### Reconstructed halos
### ################
### for k1_i, k1 in enumerate(np.linspace(0.02, 0.2, 10)):
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
###     
###     ell1 = 0
###     ell2 = 0
###     ELL = 0
###     zbin = 3
###     ax[0].axvline(k1, color="black", ls="--")
### 
###     k, bk = np.loadtxt("results_test_Rockstar_Mmin13.0_Mmax13.5_zbin%d_RSDTrue/0001/bk%d%d%d_%02d.dat" % (zbin, ell1, ell2, ELL, k1_i), usecols=(1,2), unpack=True)
###     ax[0].plot(k, k*bk, "o-", label="bk%d%d%d: pre-recon halo w/ RSDs" % (ell1, ell2, ELL))
###     k, bk = np.loadtxt("results_test_Rockstar_Mmin13.0_Mmax13.5_zbin%d_RSDTrue_recon_R15/0001/bk%d%d%d_%02d_recon.dat" % (zbin, ell1, ell2, ELL, k1_i), usecols=(1,2), unpack=True)
###     ax[0].plot(k, k*bk, "s--", label="bk%d%d%d: post-recon halo w/ RSDs" % (ell1, ell2, ELL))
###     
###     ax[0].set_xlim(0.01, 0.21)
###     ax[0].set_xticks(np.linspace(0.02, 0.2, 10))
###     ax[0].set_xlabel("k2 [h/Mpc]")
###     ax[0].set_ylabel(r"$k2\, $ bispectrum(k1, k2)")
###     ax[0].legend(loc=0)
###     ax[0].set_title("k1 = %0.2f [h/Mpc]" % k1)
###     plt.savefig("figure/bk%d%d%d_%02d_halo_recon.png" % (ell1, ell2, ELL, k1_i))
###     os.system("display figure/bk%d%d%d_%02d_halo_recon.png" % (ell1, ell2, ELL, k1_i))
### 
