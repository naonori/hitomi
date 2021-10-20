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
left = 0.10
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
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=10)
##############

#k = np.loadtxt("results_test/bk%d%d%d_00.dat" % (ell1, ell2, ELL), unpack=True, usecols=(1,))
#Nk = len(k)
#bk = np.zeros((Nk, Nk))
#for i in range(Nk):
#    bk = np.loadtxt("results_test/bk%d%d%d_%02d.dat" % (ell1, ell2, ELL, i), unpack=True, usecols=(2,))
#    if i == 0:
#        ax[0].plot(k, bk*k**(1.5), "o-", label="bk%d%d%d" % (ell1,ell2,ELL), ms=4.5, color=cm(float(9-i)/Nk), markerfacecolor="white")
#    else:
#        ax[0].plot(k, bk*k**(1.5), "o-", ms=4.5, color=cm(float(9-i)/Nk), markerfacecolor="white")
#
#X, Y = np.meshgrid(np.linspace(0.0,0.001,Nk),np.linspace(0.0,0.001,Nk))
#Z = np.zeros((Nk,Nk))
##colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]  # R -> G -> B
##cmap_name = 'my_list'
##cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=Nk)
##
#SC = ax[0].pcolor(X, Y, Z, cmap=cm, vmin=0.01, vmax=0.21, shading='auto')
#cbar = plt.colorbar(SC, ax=ax[0], ticks=np.linspace(0.20,0.02,10), pad=0.02, location="right")
#cbar.set_label(r"$k_1\ \left[h\, {\rm Mpc}^{-1} \right]$",labelpad= 30.0,rotation=270)
#cbar.set_ticklabels(["0.02","0.04","0.06","0.08","0.10","0.12","0.14","0.16","0.18","0.20"])
##cbar.set_ticklabels(["0.20","0.18","0.16","0.14","0.12","0.10","0.08","0.06","0.04","0.02"])
#
#plt.legend(loc=0)
#plt.ylabel(r"$k_2^{1.5}$ Bispectrum(k1,k2)")
#plt.xlabel(r"$k_2 [h/{\rm Mpc}]$")
#plt.savefig("figure/bk%d%d%d.png" % (ell1, ell2, ELL))
#os.system("display figure/bk%d%d%d.png" % (ell1, ell2, ELL))
#
#
## reconstruction: R = 15, 10
#
#for i in range(10):
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
#    left = 0.10
#    bottom = 0.11
#    height = (total_h - (bottom + dh))
#    width  = (total_w - (left + dw))
#    
#    #left bottom, width, height
#    rect1 = [left, bottom, width, height]
#    ax.append(fig.add_axes(rect1))
#    
#    k = np.loadtxt("results_test/bk%d%d%d_00.dat" % (ell1, ell2, ELL), unpack=True, usecols=(1,))
#    bk = np.loadtxt("results_test/bk%d%d%d_%02d.dat" % (ell1, ell2, ELL, i), unpack=True, usecols=(2,))
#    ax[0].plot(k, bk*k**(1.5), "o-", label="bk%d%d%d: pre-recon" % (ell1,ell2,ELL), ms=4.5, markerfacecolor="white")
#    bk = np.loadtxt("results_test_recon_R15/bk%d%d%d_%02d_recon.dat" % (ell1, ell2, ELL, i), unpack=True, usecols=(2,))
#    ax[0].plot(k, bk*k**(1.5), "o-", label="bk%d%d%d: R = 15 Mpc/h" % (ell1,ell2,ELL), ms=4.5, markerfacecolor="white")
#    bk = np.loadtxt("results_test_recon_R10/bk%d%d%d_%02d_recon.dat" % (ell1, ell2, ELL, i), unpack=True, usecols=(2,))
#    ax[0].plot(k, bk*k**(1.5), "o-", label="bk%d%d%d: R = 10 Mpc/h" % (ell1,ell2,ELL), ms=4.5, markerfacecolor="white")
#    
#    plt.legend(loc=0)
#    ax[0].set_title("k_1 = %1.2f [h/Mpc]" % k[i])
#    plt.ylabel(r"$k_2^{1.5}$ Bispectrum(k1,k2)")
#    plt.xlabel(r"$k_2 [h/{\rm Mpc}]$")
#    plt.savefig("figure/bk%d%d%d_recon_%02d.png" % (ell1, ell2, ELL, i))
#    os.system("display figure/bk%d%d%d_recon_%02d.png" % (ell1, ell2, ELL, i))
#    
#
