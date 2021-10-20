########################
# Purpose of this code #
########################
#
# This code aims to plot the measured 3PCFs.
#
######################################################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
import matplotlib.pyplot as plt

for i in range(13):

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
    
    r, zeta = np.loadtxt("results_test/0001/zeta000_%02d.dat" % i, usecols=(1,2), unpack=True)
    ax[0].plot(r, zeta*r**2, "-o", label="zeta000: # 0001", markerfacecolor="white")
    r, zeta = np.loadtxt("results_test/0101/zeta000_%02d.dat" % i, usecols=(1,2), unpack=True)
    ax[0].plot(r, zeta*r**2, "-o", label="zeta000: # 0101", markerfacecolor="white")
    r, zeta = np.loadtxt("results_test/1001/zeta000_%02d.dat" % i, usecols=(1,2), unpack=True)
    ax[0].plot(r, zeta*r**2, "-o", label="zeta000: # 1001", markerfacecolor="white")
    r, zeta = np.loadtxt("results_test/2001/zeta000_%02d.dat" % i, usecols=(1,2), unpack=True)
    ax[0].plot(r, zeta*r**2, "-o", label="zeta000: # 2001", markerfacecolor="white")
    
    ax[0].legend(loc=0)
    ax[0].set_xlabel("r2 [Mpc/h]")
    ax[0].set_ylabel("3PCFs")
    ax[0].set_title("r1 = %2.0f [Mpc/h]" % r[i])
    plt.savefig("figure/zeta_test_mock_%02d.png" % i)
    os.system("display figure/zeta_test_mock_%02d.png" % i)


