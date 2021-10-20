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

ell1 = 0
ell2 = 0
ELL  = 0

params = {
#    'font.family': 'Times New Roman',
#    'text.usetex': 'True',
    'font.size': 11.0,
}

plt.rcParams.update(params)

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

##############
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
colors = [(0, 0, 1), (0, 1, 0), (1, 0, 0)]  # R -> G -> B
cmap_name = 'my_list'
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=13)
##############

## 3PCF ##

#r = np.loadtxt("results_test/zeta%d%d%d_00.dat" % (ell1, ell2, ELL), unpack=True, usecols=(1,))
#Nr = len(r)
#zeta = np.zeros((Nr, Nr))
#for i in range(Nr):
#    zeta = np.loadtxt("results_test/zeta%d%d%d_%02d.dat" % (ell1, ell2, ELL, i), unpack=True, usecols=(2,))
#    if i == 0:
#        ax[0].plot(r, zeta*r**2, "o-", label="zeta%d%d%d" % (ell1,ell2,ELL), ms=4.5, color=cm(float(Nr-i-1)/Nr), markerfacecolor="white")
#    else:
#        ax[0].plot(r, zeta*r**2, "o-", ms=4.5, color=cm(float(Nr-i-1)/Nr), markerfacecolor="white")
#
#X, Y = np.meshgrid(np.linspace(0.0,1.0e-10,Nr),np.linspace(0.0,1.0e-10,Nr))
#Z = np.zeros((Nr,Nr))
#SC = ax[0].pcolor(X, Y, Z, cmap=cm, vmin=25.0, vmax=155.0, shading='auto')
#cbar = plt.colorbar(SC, ax=ax[0], ticks=np.linspace(30,150,13), pad=0.02, location="right")
#cbar.set_label(r"$r_1\ \left[{\rm Mpc}/h \right]$",labelpad= 30.0,rotation=270)
#cbar.set_ticklabels(["150", "140", "130", "120", "110", "100", "90", "80", "70", "60", "50", "40", "30"])
#
#plt.legend(loc=0)
#plt.ylabel(r"$r_2^{2}$ 3PCF(r1,r2)")
#plt.xlabel(r"$r_2 [{\rm Mpc}/h]$")
#plt.savefig("figure/zeta%d%d%d.png" % (ell1, ell2, ELL))
#os.system("display figure/zeta%d%d%d.png" % (ell1, ell2, ELL))
#

## window 3PCF ##

#r = np.loadtxt("results_test_window/zeta%d%d%d_00_window.dat" % (ell1, ell2, ELL), unpack=True, usecols=(1,))
#Nr = len(r)
#zeta = np.zeros((Nr, Nr))
#for i in range(Nr):
#    w000 = np.loadtxt("results_test_window/zeta000_%02d_window.dat" % (i), unpack=True, usecols=(2,))
#    zeta = np.loadtxt("results_test_window/zeta%d%d%d_%02d_window.dat" % (ell1, ell2, ELL, i), unpack=True, usecols=(2,))
#    zeta = zeta / w000
#    if i == 0:
#        ax[0].plot(r, zeta, "o-", label="window zeta%d%d%d/zeta000" % (ell1,ell2,ELL), ms=4.5, color=cm(float(i)/Nr), markerfacecolor="white")
#    else:
#        ax[0].plot(r, zeta, "o-", ms=4.5, color=cm(float(i)/Nr), markerfacecolor="white")
#
#X, Y = np.meshgrid(np.linspace(0.0,1.0e-10,Nr),np.linspace(0.0,1.0e-10,Nr))
#Z = np.zeros((Nr,Nr))
#SC = ax[0].pcolor(X, Y, Z, cmap=cm, vmin=25.0, vmax=155.0, shading='auto')
#cbar = plt.colorbar(SC, ax=ax[0], ticks=np.linspace(30,150,13), pad=0.02, location="right")
#cbar.set_label(r"$r_1\ \left[{\rm Mpc}/h \right]$",labelpad= 30.0,rotation=270)
##cbar.set_ticklabels(["150", "140", "130", "120", "110", "100", "90", "80", "70", "60", "50", "40", "30"])
#cbar.set_ticklabels(["30","40","50","60","70","80","90","100","110","120","130","140","150"])
#
#plt.legend(loc=0)
#plt.ylabel(r"window 3PCF(r1,r2)")
#plt.xlabel(r"$r_2 [{\rm Mpc}/h]$")
#plt.savefig("figure/zeta%d%d%d_window.png" % (ell1, ell2, ELL))
#os.system("display figure/zeta%d%d%d_window.png" % (ell1, ell2, ELL))
#

## reconstruction: R = 10, 15
#
#for i in range(13):
#
#    fig = plt.figure(figsize=(7.0,7.0))
#    ax = []
#    
#    ##--------------------------
#    total_h = 0.94
#    total_w = 0.99
#    
#    dw = 0.07
#    dh = 0.00
#    left = 0.15
#    bottom = 0.11
#    height = (total_h - (bottom + dh))
#    width  = (total_w - (left + dw))
#    
#    #left bottom, width, height
#    rect1 = [left, bottom, width, height]
#    ax.append(fig.add_axes(rect1))
#
#    r = np.loadtxt("results_test/zeta%d%d%d_00.dat" % (ell1, ell2, ELL), unpack=True, usecols=(1,))
#    zeta = np.loadtxt("results_test/zeta%d%d%d_%02d.dat" % (ell1, ell2, ELL, i), unpack=True, usecols=(2,))
#    ax[0].plot(r, zeta*r**2, "o-", label="zeta%d%d%d: pre-recon" % (ell1,ell2,ELL), ms=4.5, markerfacecolor="white")
#    zeta = np.loadtxt("results_test_recon_R15/zeta%d%d%d_%02d_recon.dat" % (ell1, ell2, ELL, i), unpack=True, usecols=(2,))
#    ax[0].plot(r, zeta*r**2, "o-", label="zeta%d%d%d: R = 15 Mpc/h" % (ell1,ell2,ELL), ms=4.5, markerfacecolor="white")
#    zeta = np.loadtxt("results_test_recon_R10/zeta%d%d%d_%02d_recon.dat" % (ell1, ell2, ELL, i), unpack=True, usecols=(2,))
#    ax[0].plot(r, zeta*r**2, "o-", label="zeta%d%d%d: R = 10 Mpc/h" % (ell1,ell2,ELL), ms=4.5, markerfacecolor="white")
#
#    ax[0].set_title("r1 = %2.0f [Mpc/h]" % r[i])
#    plt.legend(loc=0)
#    plt.ylabel(r"$r_2^{2}$ 3PCF(r1,r2)")
#    plt.xlabel(r"$r_2 [{\rm Mpc}/h]$")
#    plt.savefig("figure/zeta%d%d%d_recon_%02d.png" % (ell1, ell2, ELL, i))
#    os.system("display figure/zeta%d%d%d_recon_%02d.png" % (ell1, ell2, ELL, i))
#    
#    
#
## reconstructed window 3PCF: R = 10, 15
#for i in range(13):
#
#    fig = plt.figure(figsize=(7.0,7.0))
#    ax = []
#    
#    ##--------------------------
#    total_h = 0.94
#    total_w = 0.99
#    
#    dw = 0.07
#    dh = 0.00
#    left = 0.15
#    bottom = 0.11
#    height = (total_h - (bottom + dh))
#    width  = (total_w - (left + dw))
#    
#    #left bottom, width, height
#    rect1 = [left, bottom, width, height]
#    ax.append(fig.add_axes(rect1))
#
#    r = np.loadtxt("results_test_window/zeta%d%d%d_00_window.dat" % (ell1, ell2, ELL), unpack=True, usecols=(1,))
#    zeta = np.loadtxt("results_test_window/zeta%d%d%d_%02d_window.dat" % (ell1, ell2, ELL, i), unpack=True, usecols=(2,))
#    ax[0].plot(r, zeta, "o-", label="zeta%d%d%d: pre-recon" % (ell1,ell2,ELL), ms=4.5, markerfacecolor="white")
#    zeta = np.loadtxt("results_test_recon_R15_window/zeta%d%d%d_%02d_recon_window.dat" % (ell1, ell2, ELL, i), unpack=True, usecols=(2,))
#    ax[0].plot(r, zeta, "o-", label="zeta%d%d%d: R = 15 Mpc/h" % (ell1,ell2,ELL), ms=4.5, markerfacecolor="white")
#    zeta = np.loadtxt("results_test_recon_R10_window/zeta%d%d%d_%02d_recon_window.dat" % (ell1, ell2, ELL, i), unpack=True, usecols=(2,))
#    ax[0].plot(r, zeta, "o-", label="zeta%d%d%d: R = 10 Mpc/h" % (ell1,ell2,ELL), ms=4.5, markerfacecolor="white")
#
#    ax[0].set_title("r1 = %2.0f [Mpc/h]" % r[i])
#    plt.legend(loc=0)
#    plt.ylabel(r"window 3PCF(r1,r2)")
#    plt.xlabel(r"$r_2 [{\rm Mpc}/h]$")
#    plt.savefig("figure/zeta%d%d%d_%02d_window.png" % (ell1, ell2, ELL, i))
#    os.system("display figure/zeta%d%d%d_%02d_window.png" % (ell1, ell2, ELL, i))
#    
#    
#
