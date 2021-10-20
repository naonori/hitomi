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
fig = plt.figure(figsize=(12.0,7.0))
ax = []

##--------------------------
total_h = 0.94
total_w = 0.99

dw = 0.35
dh = 0.00
left = 0.13
bottom = 0.11
height = (total_h - (bottom + dh))
width  = (total_w - (left + dw))

#left bottom, width, height
rect1 = [left, bottom, width, height]
ax.append(fig.add_axes(rect1))

### ##############
### ## Dark matter
### #############
### ELL = 0
### RSD = "False"
### if RSD == "False":
###     label_RSD = "w/o RSD"
### elif RSD == "True":
###     label_RSD = "w/  RSD"
### 
### for zbin, z in enumerate([2.0, 1.0, 0.5, 0.35, 0.0]):
###     r, xi = np.loadtxt("results_test_Gadget_zbin%d_RSD%s/0001/xi%d.dat" % (zbin, RSD, ELL), usecols=(0,1), unpack=True)
###     ax[0].plot(r, xi*r**2, "o", label="xi%d: z = %2.1f %s" % (ELL, z, label_RSD))
### 
### ax[0].legend(loc=0)
### ax[0].set_xlabel("r [h/Mpc]")
### ax[0].set_ylabel("$r^{2}$ 2PCF(r)")
### plt.savefig("figure/xi%d_RSD%s_Gadget.png" % (ELL, RSD))
### os.system("display figure/xi%d_RSD%s_Gadget.png" % (ELL, RSD))
### 
### ##############
### ## Halos
### #############
### ELL = 0
### RSD = "False"
### if RSD == "False":
###     label_RSD = "w/o RSD"
### elif RSD == "True":
###     label_RSD = "w/  RSD"
### 
### color = ["royalblue", "orange", "forestgreen", "mediumorchid", "orangered"]
### ms = ["o", "^", "s"]
### 
### for zbin, z in enumerate([2.0, 1.0, 0.5, 0.35, 0.0]):
### 
###     for i_M, (Mmin, Mmax) in enumerate([(12.5, 13.0), (13.0, 13.5), (13.5, 14.0)]):
###         if zbin == 0 and (Mmin, Mmax) == (13.5, 14.0):
###             continue
### 
###         r, xi = np.loadtxt("results_test_Rockstar_Mmin%s_Mmax%s_zbin%d_RSD%s/0001/xi%d.dat" % (Mmin, Mmax, zbin, RSD, ELL), usecols=(0,1), unpack=True)
###         ax[0].plot(r, xi*r**2, "-%s" % ms[i_M], color=color[zbin], label="xi%d: z = %2.1f %s" % (ELL, z, label_RSD))
### 
### ax[0].legend(loc="upper left", bbox_to_anchor=(1,1))
### ax[0].set_xlabel("r [h/Mpc]")
### ax[0].set_ylabel("$r^{2}$ 2PCF(r)")
### plt.savefig("figure/xi%d_RSD%s_Rockstar.png" % (ELL, RSD))
### os.system("display figure/xi%d_RSD%s_Rockstar.png" % (ELL, RSD))
### 
### ##############
### ## Reconstructed Halos
### #############
### ELL = 0
### RSD = "False"
### if RSD == "False":
###     label_RSD = "w/o RSD"
### elif RSD == "True":
###     label_RSD = "w/  RSD"
### 
### color = ["royalblue", "orange", "forestgreen", "mediumorchid", "orangered"]
### ms = ["o", "^", "s"]
### 
### for zbin, z in enumerate([2.0, 1.0, 0.5, 0.35, 0.0]):
### 
###     for i_M, (Mmin, Mmax) in enumerate([(12.5, 13.0), (13.0, 13.5), (13.5, 14.0)]):
###         if zbin == 0 and (Mmin, Mmax) == (13.5, 14.0):
###             continue
### 
###         r, xi = np.loadtxt("results_test_Rockstar_Mmin%s_Mmax%s_zbin%d_RSD%s_recon_R15/0001/xi%d_recon.dat" % (Mmin, Mmax, zbin, RSD, ELL), usecols=(0,1), unpack=True)
###         ax[0].plot(r, xi*r**2, "-%s" % ms[i_M], color=color[zbin], label="xi%d: z = %2.1f %s" % (ELL, z, label_RSD))
### 
### ax[0].legend(loc="upper left", bbox_to_anchor=(1,1))
### ax[0].set_xlabel("r [h/Mpc]")
### ax[0].set_ylabel("$r^{2}$ 2PCF(r)")
### plt.savefig("figure/xi%d_RSD%s_Rockstar_recon_R15.png" % (ELL, RSD))
### os.system("display figure/xi%d_RSD%s_Rockstar_recon_R15.png" % (ELL, RSD))
### 
### 
