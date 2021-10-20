########################
# Purpose of this code #
########################
#
# This code aims to plot the measured 2PCFs.
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
left = 0.15
bottom = 0.11
height = (total_h - (bottom + dh))
width  = (total_w - (left + dw))

#left bottom, width, height
rect1 = [left, bottom, width, height]
ax.append(fig.add_axes(rect1))

## multipoles
#r, xi0 = np.loadtxt("results_test/xi0.dat", usecols=(0,1), unpack=True)
#ax[0].plot(r, xi0*r**2, "o", color="royalblue", label="xi0")
#r, xi2 = np.loadtxt("results_test/xi2.dat", usecols=(0,1), unpack=True)
#ax[0].plot(r, xi2*r**2, "o", color="magenta", label="xi2")
#r, xi4 = np.loadtxt("results_test/xi4.dat", usecols=(0,1), unpack=True)
#ax[0].plot(r, xi4*r**2, "o", color="forestgreen", label="xi4")
#
#ax[0].legend(loc=0)
#ax[0].set_xlabel("r [Mpc/h]")
#ax[0].set_ylabel("2PCFs")
#plt.savefig("figure/xi_test.png")
#os.system("display figure/xi_test.png")

## wider rbins for xi0
#r, xi0 = np.loadtxt("results_test/xi0.dat", usecols=(0,1), unpack=True)
#ax[0].plot(r, xi0*r**2, "o", color="royalblue", label="dr = 5 [Mpc/h]")
#r, xi0 = np.loadtxt("results_test_wider_rbin/xi0.dat", usecols=(0,1), unpack=True)
#ax[0].plot(r, xi0*r**2, "o", color="magenta", label="dr = 10 [Mpc/h]")
#ax[0].legend(loc=0)
#ax[0].set_xlabel("r [Mpc/h]")
#ax[0].set_ylabel("2PCFs")
#plt.savefig("figure/xi_test_wider_rbin.png")
#os.system("display figure/xi_test_wider_rbin.png")

## window 2PCFs
#r, xi0 = np.loadtxt("results_test_window/xi0_window.dat", usecols=(0,1), unpack=True)
#r, xi2 = np.loadtxt("results_test_window/xi2_window.dat", usecols=(0,1), unpack=True)
#ax[0].plot(r, xi2/xi0, "-o", color="magenta", label="window xi2/xi0")
#r, xi4 = np.loadtxt("results_test_window/xi4_window.dat", usecols=(0,1), unpack=True)
#ax[0].plot(r, xi4/xi0, "-o", color="forestgreen", label="window xi4/xi0")
#
#ax[0].legend(loc=0)
#ax[0].set_xlabel("r [Mpc/h]")
#ax[0].set_ylabel("2PCFs")
#plt.savefig("figure/xi_window.png")
#os.system("display figure/xi_window.png")
#
#
## reconstruction xi0 or xi2: R = 0, 5, 10, 15 20
#
#ELL = 0
#
#for R in [0,5,10,15,20]:
#    r, xi = np.loadtxt("results_test_recon_R%02d/xi%d_recon.dat" % (R,ELL), usecols=(0,1), unpack=True)
#    ax[0].plot(r, xi*r**2, "o-", label="xi%d: R = %02d" % (ELL, R), markerfacecolor="white", ms=3.5)
#r, xi = np.loadtxt("results_test/xi%d.dat" % ELL, usecols=(0,1), unpack=True)
#ax[0].plot(r, xi*r**2, "o-", label="xi%d: pre-recon" % ELL, markerfacecolor="white", ms=3.5)
#
#ax[0].legend(loc=0, ncol=2)
#ax[0].set_xlabel("r [Mpc/h]")
#ax[0].set_ylabel("2PCF")
#plt.savefig("figure/xi%d_recon.png" % ELL)
#os.system("display figure/xi%d_recon.png" % ELL)
#
#

## recon. window 2PCF: R = 0, 5, 10, 15 20
#
#ELL = 2
#
#for R in [0,5,10,15,20]:
#    r, xi = np.loadtxt("results_test_recon_R%02d_window/xi%d_recon_window.dat" % (R,ELL), usecols=(0,1), unpack=True)
#    ax[0].plot(r, xi, "o-", label="xi%d: R = %02d" % (ELL, R), markerfacecolor="white", ms=3.5)
#r, xi = np.loadtxt("results_test_window/xi%d_window.dat" % ELL, usecols=(0,1), unpack=True)
#ax[0].plot(r, xi, "o-", label="xi%d: pre-recon" % ELL, markerfacecolor="white", ms=3.5)
#
#ax[0].legend(loc=0, ncol=2)
#ax[0].set_xlabel("r [Mpc/h]")
#ax[0].set_ylabel("window 2PCFs")
#plt.savefig("figure/xi%d_recon_window.png" % ELL)
#os.system("display figure/xi%d_recon_window.png" % ELL)
#
#
