########################
# Purpose of this code #
########################
#
# This code aims to plot the location of the particles and intuitively know how the galaxies are distributed. The horizontal axis is the comoving distance to the galaxy, and the vertical axis is the x-, y-, or z-axis in Cartesian coordinates.
#
######################################################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import random

NS = "South" # North or South

params = {
#    'font.family': 'Times New Roman',
#    'text.usetex': 'True',
    'font.size': 11.0,
}

plt.rcParams.update(params)
#fig = plt.figure(figsize=(210.0/25.4, 264.0/25.4/2.3))
fig = plt.figure(figsize=(1.8*210.0/25.4, 1.8*210.0/25.4/3.0))
ax = []

##--------------------------
total_h = 0.94
total_w = 0.99

dw = 0.07
dh = 0.00
left = 0.07
bottom = 0.11
height = (total_h - (bottom + dh))
width  = (total_w - (left + 2*dw)) / 3.0

#left bottom, width, height
rect1 = [left + 0*(width + dw), bottom, width, height]
rect2 = [left + 1*(width + dw), bottom, width, height]
rect3 = [left + 2*(width + dw), bottom, width, height]

ax.append(fig.add_axes(rect1))
ax.append(fig.add_axes(rect2))
ax.append(fig.add_axes(rect3))

#for zbin in [1,3]:
#    
#    x, y, z = np.loadtxt("galaxy_DR12v5_CMASSLOWZTOT/galaxy_DR12v5_CMASSLOWZTOT_%s_ZBIN%d.dat" % (NS, zbin), unpack=True, usecols=(0,1,2))
#    
#    N = 1000
#    xR = np.zeros(N)
#    yR = np.zeros(N)
#    zR = np.zeros(N)
#    cR = np.zeros(N)
#    for i in range(N):
#        ran = random.randint(0, len(x))
#        xR[i] = x[ran]
#        yR[i] = y[ran]
#        zR[i] = z[ran]
#        cR[i] = np.sqrt(xR[i]**2 + yR[i]**2 + zR[i]**2)
# 
#    ax[0].plot(cR, xR, ".", label="zbin = %d" % zbin)
#    ax[1].plot(cR, yR, ".")
#    ax[2].plot(cR, zR, ".")
#
#ax[0].set_ylabel("x [Mpc/h]")
#ax[1].set_ylabel("y [Mpc/h]")
#ax[2].set_ylabel("z [Mpc/h]")
#ax[0].set_xlabel("comoving distance [Mpc/h]")
#ax[1].set_xlabel("comoving distance [Mpc/h]")
#ax[2].set_xlabel("comoving distance [Mpc/h]")
#ax[0].legend(loc=0)
#ax[1].set_title("galaxy %s" % NS)
#
#plt.savefig("figure/galaxy_%s.png" % NS)
#os.system("display figure/galaxy_%s.png" % NS)

##########
## random particles
##########
## Check that the random particles are distributed in the same region as the galaxy sample.
#for zbin in [1,3]:
#    
#    x, y, z = np.loadtxt("galaxy_DR12v5_CMASSLOWZTOT/random_DR12v5_CMASSLOWZTOT_%s_ZBIN%d.dat" % (NS, zbin), unpack=True, usecols=(0,1,2))
#    
#    N = 3000
#    xR = np.zeros(N)
#    yR = np.zeros(N)
#    zR = np.zeros(N)
#    cR = np.zeros(N)
#    for i in range(N):
#        ran = random.randint(0, len(x))
#        xR[i] = x[ran]
#        yR[i] = y[ran]
#        zR[i] = z[ran]
#        cR[i] = np.sqrt(xR[i]**2 + yR[i]**2 + zR[i]**2)
# 
#    ax[0].plot(cR, xR, ".", color="darkgray")
#    ax[1].plot(cR, yR, ".", color="darkgray")
#    ax[2].plot(cR, zR, ".", color="darkgray")
#
#plt.savefig("figure/random_%s.png" % NS)
#os.system("display figure/random_%s.png" % NS)

##########
## mock catalogues
#########
## See the difference between two mock particle distributions with different realizations.
#for zbin in [1,3]:
#    for realization in [1, 2048]:
#
#        if NS == "North":
#            x, y, z = np.loadtxt("Patchy-Mocks-DR12NGC-COMPSAM_V6C_ZBIN%d/Patchy-Mocks-DR12NGC-COMPSAM_V6C_ZBIN%d_%04d.dat" % (zbin, zbin, realization), unpack=True, usecols=(0,1,2))
#        elif NS == "South":                                              
#            x, y, z = np.loadtxt("Patchy-Mocks-DR12SGC-COMPSAM_V6C_ZBIN%d/Patchy-Mocks-DR12SGC-COMPSAM_V6C_ZBIN%d_%04d.dat" % (zbin, zbin, realization), unpack=True, usecols=(0,1,2))
#        
#        N = 1000
#        xR = np.zeros(N)
#        yR = np.zeros(N)
#        zR = np.zeros(N)
#        cR = np.zeros(N)
#        for i in range(N):
#            ran = random.randint(0, len(x))
#            xR[i] = x[ran]
#            yR[i] = y[ran]
#            zR[i] = z[ran]
#            cR[i] = np.sqrt(xR[i]**2 + yR[i]**2 + zR[i]**2)
#     
#        if realization == 1:
#            color = "royalblue"
#        else:
#            color = "magenta"
#
#        if zbin == 1:
#            ax[0].plot(cR, xR, ".", color=color, label="# %04d" % realization)
#        else:
#            ax[0].plot(cR, xR, ".", color=color)
#        ax[1].plot(cR, yR, ".", color=color)
#        ax[2].plot(cR, zR, ".", color=color)
#
#ax[0].set_ylabel("x [Mpc/h]")
#ax[1].set_ylabel("y [Mpc/h]")
#ax[2].set_ylabel("z [Mpc/h]")
#ax[0].set_xlabel("comoving distance [Mpc/h]")
#ax[1].set_xlabel("comoving distance [Mpc/h]")
#ax[2].set_xlabel("comoving distance [Mpc/h]")
#ax[0].legend(loc=0)
#ax[1].set_title("mock %s" % NS)

#plt.savefig("figure/mock_%s.png" % NS)
#os.system("display figure/mock_%s.png" % NS)

##########
## mock random catalogues
#########
#for zbin in [1,3]:
#
#    if NS == "North":
#        x, y, z = np.loadtxt("Patchy-Mocks-DR12NGC-COMPSAM_V6C_ZBIN%d/Patchy-Mocks-Randoms-DR12NGC-COMPSAM_V6C_x100_ZBIN%d.dat" % (zbin, zbin), unpack=True, usecols=(0,1,2))
#    elif NS == "South":                                              
#        x, y, z = np.loadtxt("Patchy-Mocks-DR12SGC-COMPSAM_V6C_ZBIN%d/Patchy-Mocks-Randoms-DR12SGC-COMPSAM_V6C_x100_ZBIN%d.dat" % (zbin, zbin), unpack=True, usecols=(0,1,2))
# 
#    N = 3000
#    xR = np.zeros(N)
#    yR = np.zeros(N)
#    zR = np.zeros(N)
#    cR = np.zeros(N)
#    for i in range(N):
#        ran = random.randint(0, len(x))
#        xR[i] = x[ran]
#        yR[i] = y[ran]
#        zR[i] = z[ran]
#        cR[i] = np.sqrt(xR[i]**2 + yR[i]**2 + zR[i]**2)
#    
#    color = "darkgray"
#    if zbin == 1:
#        ax[0].plot(cR, xR, ".", color=color, label="random")
#    else:
#        ax[0].plot(cR, xR, ".", color=color)
#    ax[1].plot(cR, yR, ".", color=color)
#    ax[2].plot(cR, zR, ".", color=color)
#    
#    ax[0].set_ylabel("x [Mpc/h]")
#    ax[1].set_ylabel("y [Mpc/h]")
#    ax[2].set_ylabel("z [Mpc/h]")
#    ax[0].set_xlabel("comoving distance [Mpc/h]")
#    ax[1].set_xlabel("comoving distance [Mpc/h]")
#    ax[2].set_xlabel("comoving distance [Mpc/h]")
#    ax[0].legend(loc=0)
#    ax[1].set_title("mock %s" % NS)
#
#plt.savefig("figure/mock_random_%s.png" % NS)
#os.system("display figure/mock_random_%s.png" % NS)
#
#
