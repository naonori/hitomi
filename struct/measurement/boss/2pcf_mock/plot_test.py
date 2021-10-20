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

#r, xi = np.loadtxt("results_test/0001/xi0.dat", usecols=(0,1), unpack=True)
#ax[0].plot(r, xi*r**2, "-o", label="xi0: # 0001", markerfacecolor="white")
#r, xi = np.loadtxt("results_test/0100/xi0.dat", usecols=(0,1), unpack=True)
#ax[0].plot(r, xi*r**2, "-o", label="xi0: # 0100", markerfacecolor="white")
#r, xi = np.loadtxt("results_test/1000/xi0.dat", usecols=(0,1), unpack=True)
#ax[0].plot(r, xi*r**2, "-o", label="xi0: # 1000", markerfacecolor="white")
#r, xi = np.loadtxt("results_test/2000/xi0.dat", usecols=(0,1), unpack=True)
#ax[0].plot(r, xi*r**2, "-o", label="xi0: # 2000", markerfacecolor="white")
#
#ax[0].legend(loc=0)
#ax[0].set_xlabel("r [Mpc/h]")
#ax[0].set_ylabel("2PCFs")
#plt.savefig("figure/xi_test_mock.png")
#os.system("display figure/xi_test_mock.png")
#
#
