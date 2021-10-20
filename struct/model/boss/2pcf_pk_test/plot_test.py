########################
# Purpose of this code #
########################
#
# This code aims to plot the computed power spectra and 2PCFs.
#
######################################################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
import matplotlib.pyplot as plt

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
left = 0.13
bottom = 0.11
height = (total_h - (bottom + dh))
width  = (total_w - (left + dw))

#left bottom, width, height
rect1 = [left, bottom, width, height]
ax.append(fig.add_axes(rect1))

### ##########
### ## Linear power spectra
### 
### for ELL in [0,2,4]:
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
### 
###     k, pk = np.loadtxt("results_test/pk%d_Tree_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(k, k**1.5*pk, "-", color="royalblue", label="pk%d" % ELL)
###     ax[0].set_xlim(0.01, 0.2)
###     ax[0].set_xlabel("k [h/Mpc]")
###     ax[0].set_ylabel(r"$k^{1.5}$ power spectrum(k)")
###     ax[0].legend(loc=0)
###     plt.savefig("figure/pk%d_Tree.png" % ELL)
###     os.system("display figure/pk%d_Tree.png" % ELL)
### 
### ##########
### ## Linear 2PCFs
### 
### for ELL in [0,2,4]:
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
### 
###     r, xi = np.loadtxt("results_test/xi%d_Tree_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(r, r**2*xi, "-", color="royalblue", label="xi%d" % ELL)
###     ax[0].set_xlim(40, 150)
###     ax[0].set_xlabel("r [Mpc/h]")
###     ax[0].set_ylabel(r"$r^2$ 2PCF(r)")
###     ax[0].legend(loc=0)
###     plt.savefig("figure/xi%d_Tree.png" % ELL)
###     os.system("display figure/xi%d_Tree.png" % ELL)
### 
### 
### ##########
### ## Linear no-wiggle power spectra
### ###
### for ELL in [0,2,4]:
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
### 
###     k, pk = np.loadtxt("results_test/pk%d_Tree_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(k, k**1.5*pk, "-", color="royalblue", label="pk%d: wiggle" % ELL)
###     k, pk = np.loadtxt("results_test/pk%d_Tree_NoWiggle_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(k, k**1.5*pk, "-", color="forestgreen", label="pk%d: no-wiggle" % ELL)
###     ax[0].set_xlim(0.01, 0.2)
###     ax[0].set_xlabel("k [h/Mpc]")
###     ax[0].set_ylabel(r"$k^{1.5}$ power spectrum(k)")
###     ax[0].legend(loc=0)
###     plt.savefig("figure/pk%d_Tree_NoWiggle.png" % ELL)
###     os.system("display figure/pk%d_Tree_NoWiggle.png" % ELL)
### 
### ##########
### ## Linear no-wiggle 2PCFs
### 
### for ELL in [0,2,4]:
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
### 
###     r, xi = np.loadtxt("results_test/xi%d_Tree_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(r, r**2*xi, "-", color="royalblue", label="xi%d: wiggle" % ELL)
###     r, xi = np.loadtxt("results_test/xi%d_Tree_NoWiggle_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(r, r**2*xi, "-", color="forestgreen", label="xi%d: no-wiggle" % ELL)
###     ax[0].set_xlim(40, 150)
###     ax[0].set_xlabel("r [Mpc/h]")
###     ax[0].set_ylabel(r"$r^2$ 2PCF(r)")
###     ax[0].legend(loc=0)
###     plt.savefig("figure/xi%d_Tree_NoWiggle.png" % ELL)
###     os.system("display figure/xi%d_Tree_NoWiggle.png" % ELL)
### 
###
### ##########
### ## Linear power spectra with non-linear BAO
### ###
### for ELL in [0,2,4]:
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
### 
###     k, pk = np.loadtxt("results_test/pk%d_Tree_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(k, k**1.5*pk, "-", color="royalblue", label="pk%d: wiggle" % ELL)
###     k, pk = np.loadtxt("results_test/pk%d_Tree_NoWiggle_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(k, k**1.5*pk, "-", color="forestgreen", label="pk%d: no-wiggle" % ELL)
###     k, pk = np.loadtxt("results_test/pk%d_Tree_BAO_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(k, k**1.5*pk, "-", color="magenta", label="pk%d: NL BAO" % ELL)
###     ax[0].set_xlim(0.01, 0.2)
###     ax[0].set_xlabel("k [h/Mpc]")
###     ax[0].set_ylabel(r"$k^{1.5}$ power spectrum(k)")
###     ax[0].legend(loc=0)
###     plt.savefig("figure/pk%d_Tree_BAO.png" % ELL)
###     os.system("display figure/pk%d_Tree_BAO.png" % ELL)
### 
### ##########
### ## Linear 2PCFs with non-linear BAO
### 
### for ELL in [0,2,4]:
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
### 
###     r, xi = np.loadtxt("results_test/xi%d_Tree_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(r, r**2*xi, "-", color="royalblue", label="xi%d: wiggle" % ELL)
###     r, xi = np.loadtxt("results_test/xi%d_Tree_NoWiggle_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(r, r**2*xi, "-", color="forestgreen", label="xi%d: no-wiggle" % ELL)
###     r, xi = np.loadtxt("results_test/xi%d_Tree_BAO_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(r, r**2*xi, "-", color="magenta", label="xi%d: NL BAO" % ELL)
###     ax[0].set_xlim(40, 150)
###     ax[0].set_xlabel("r [Mpc/h]")
###     ax[0].set_ylabel(r"$r^2$ 2PCF(r)")
###     ax[0].legend(loc=0)
###     plt.savefig("figure/xi%d_Tree_BAO.png" % ELL)
###     os.system("display figure/xi%d_Tree_BAO.png" % ELL)
### 
### 
### 
###
### ##########
### ## Linear power spectra with reconstructed non-linear BAO
### ###
### for ELL in [0,2,4]:
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
### 
###     k, pk = np.loadtxt("results_test/pk%d_Tree_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(k, k**1.5*pk, "-", color="royalblue", label="pk%d: wiggle" % ELL)
###     k, pk = np.loadtxt("results_test/pk%d_Tree_NoWiggle_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(k, k**1.5*pk, "-", color="forestgreen", label="pk%d: no-wiggle" % ELL)
###     k, pk = np.loadtxt("results_test/pk%d_Tree_BAO_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(k, k**1.5*pk, "-", color="magenta", label="pk%d: NL BAO" % ELL)
### 
###     for R in [5, 10, 15, 20]:
###         k, pk = np.loadtxt("results_test/pk%d_Tree_BAO_recon_R%02d_fft.dat" % (ELL, R), usecols=(0,1), unpack=True)
###         ax[0].plot(k, k**1.5*pk, "-", label="pk%d: R = %02d" % (ELL, R))
### 
###     ax[0].set_xlim(0.01, 0.2)
###     ax[0].set_xlabel("k [h/Mpc]")
###     ax[0].set_ylabel(r"$k^{1.5}$ power spectrum(k)")
###     ax[0].legend(loc=0)
###     plt.savefig("figure/pk%d_Tree_BAO_recon.png" % ELL)
###     os.system("display figure/pk%d_Tree_BAO_recon.png" % ELL)
### 
### ##########
### ## Linear 2PCFs with reconstructed non-linear BAO
### 
### for ELL in [0,2,4]:
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
### 
###     r, xi = np.loadtxt("results_test/xi%d_Tree_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(r, r**2*xi, "-", color="royalblue", label="xi%d: wiggle" % ELL)
###     r, xi = np.loadtxt("results_test/xi%d_Tree_NoWiggle_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(r, r**2*xi, "-", color="forestgreen", label="xi%d: no-wiggle" % ELL)
###     r, xi = np.loadtxt("results_test/xi%d_Tree_BAO_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(r, r**2*xi, "-", color="magenta", label="xi%d: NL BAO" % ELL)
### 
###     for R in [5, 10, 15, 20]:
###         r, xi = np.loadtxt("results_test/xi%d_Tree_BAO_recon_R%02d_fft.dat" % (ELL, R), usecols=(0,1), unpack=True)
###         ax[0].plot(r, r**2*xi, "-", label="xi%d: R = %02d" % (ELL, R))
### 
###     ax[0].set_xlim(40, 150)
###     ax[0].set_xlabel("r [Mpc/h]")
###     ax[0].set_ylabel(r"$r^2$ 2PCF(r)")
###     ax[0].legend(loc=0)
###     plt.savefig("figure/xi%d_Tree_BAO_recon.png" % ELL)
###     os.system("display figure/xi%d_Tree_BAO_recon.png" % ELL)
### 
### 
### ##########
### ## Linear power spectra with non-linear BAO, decomposed by parameters.
### ###
### for ELL in [0,2,4]:
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
### 
###     k, pk = np.loadtxt("results_test/pk%d_Tree_BAO_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(k, k**1.5*pk, "-", color="magenta", label="pk%d: total" % ELL)
###     k, pk = np.loadtxt("results_test/pk%d_Tree_BAO_Template_b1_b1_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(k, k**1.5*pk, "-", label="pk%d: b1_b1" % ELL)
###     k, pk = np.loadtxt("results_test/pk%d_Tree_BAO_Template_b1_f_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(k, k**1.5*pk, "-", label="pk%d: b1_f" % ELL)
###     k, pk = np.loadtxt("results_test/pk%d_Tree_BAO_Template_f_f_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(k, k**1.5*pk, "-", label="pk%d: f_f" % ELL)
###     ax[0].set_xlim(0.01, 0.2)
###     ax[0].set_xlabel("k [h/Mpc]")
###     ax[0].set_ylabel(r"$k^{1.5}$ power spectrum(k)")
###     ax[0].legend(loc=0)
###     plt.savefig("figure/pk%d_Tree_BAO_decon.png" % ELL)
###     os.system("display figure/pk%d_Tree_BAO_decon.png" % ELL)
### 
### ##########
### ## Linear 2PCFs with non-linear BAO, decomposed by parameters.
### 
### for ELL in [0,2,4]:
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
### 
###     r, xi = np.loadtxt("results_test/xi%d_Tree_BAO_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(r, r**2*xi, "-", color="magenta", label="xi%d: total" % ELL)
###     r, xi = np.loadtxt("results_test/xi%d_Tree_BAO_Template_b1_b1_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(r, r**2*xi, "-", label="xi%d: b1_b1" % ELL)
###     r, xi = np.loadtxt("results_test/xi%d_Tree_BAO_Template_b1_f_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(r, r**2*xi, "-", label="xi%d: b1_f" % ELL)
###     r, xi = np.loadtxt("results_test/xi%d_Tree_BAO_Template_f_f_fft.dat" % ELL, usecols=(0,1), unpack=True)
###     ax[0].plot(r, r**2*xi, "-", label="xi%d: f_f" % ELL)
###     ax[0].set_xlim(40, 150)
###     ax[0].set_xlabel("r [Mpc/h]")
###     ax[0].set_ylabel(r"$r^2$ 2PCF(r)")
###     ax[0].legend(loc=0)
###     plt.savefig("figure/xi%d_Tree_BAO_decon.png" % ELL)
###     os.system("display figure/xi%d_Tree_BAO_decon.png" % ELL)
### # 
### 
