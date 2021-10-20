########################
# Purpose of this code #
########################
#
# This code aims to plot the measured power spectra.
#
######################################################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
import matplotlib.pyplot as plt

NS = "North" # North or South

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
left = 0.10
bottom = 0.11
height = (total_h - (bottom + dh))
width  = (total_w - (left + dw))

#left bottom, width, height
rect1 = [left, bottom, width, height]
ax.append(fig.add_axes(rect1))

## multipoles
#k, pk0 = np.loadtxt("results_test/pk0.dat", usecols=(0,1), unpack=True)
#ax[0].plot(k, pk0*k**1.5, "o", color="royalblue", label="pk0")
#k, pk2 = np.loadtxt("results_test/pk2.dat", usecols=(0,1), unpack=True)
#ax[0].plot(k, pk2*k**1.5, "o", color="magenta", label="pk2")
#k, pk4 = np.loadtxt("results_test/pk4.dat", usecols=(0,1), unpack=True)
#ax[0].plot(k, pk4*k**1.5, "o", color="forestgreen", label="pk4")
#
#ax[0].legend(loc=0)
#ax[0].set_xlabel("k [h/Mpc]")
#ax[0].set_ylabel("power spectra")
#plt.savefig("figure/pk_test.png")
#os.system("display figure/pk_test.png")

## wider kbins for pk0
#k, pk0 = np.loadtxt("results_test/pk0.dat", usecols=(0,1), unpack=True)
#ax[0].plot(k, pk0*k**1.5, "o", color="royalblue", label="dk = 0.01 [h/Mpc]")
#k, pk2 = np.loadtxt("results_test_wider_kbin/pk0.dat", usecols=(0,1), unpack=True)
#ax[0].plot(k, pk2*k**1.5, "o", color="magenta", label="dk = 0.02 [h/Mpc]")
#ax[0].legend(loc=0)
#ax[0].set_xlabel("k [h/Mpc]")
#ax[0].set_ylabel("power spectra")
#plt.savefig("figure/pk_test_wider_kbin.png")
#os.system("display figure/pk_test_wider_kbin.png")


## reconstruction pk0 or pk2: R = 0, 5, 10, 15 20
#
#ELL = 2
#
#for R in [0,5,10,15,20]:
#    k, pk = np.loadtxt("results_test_recon_R%02d/pk%d_recon.dat" % (R,ELL), usecols=(0,1), unpack=True)
#    ax[0].plot(k, pk*k**1.5, "o-", label="pk%d: R = %02d" % (ELL, R), markerfacecolor="white", ms=3.5)
#k, pk = np.loadtxt("results_test/pk%d.dat" % ELL, usecols=(0,1), unpack=True)
#ax[0].plot(k, pk*k**1.5, "o-", label="pk%d: pre-recon" % ELL, markerfacecolor="white", ms=3.5)
#
#ax[0].legend(loc=0, ncol=2)
#ax[0].set_xlabel("k [h/Mpc]")
#ax[0].set_ylabel("power spectra")
#plt.savefig("figure/pk%d_recon.png" % ELL)
#os.system("display figure/pk%d_recon.png" % ELL)
#

