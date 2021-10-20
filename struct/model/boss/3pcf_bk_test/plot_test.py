########################
# Purpose of this code #
########################
#
# This code aims to plot the computed bispectra and 3PCFs.
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

### ##########
### ## Tree-level bispectra
### for k1_i, k1 in enumerate(np.linspace(0.02, 0.2, 10)):
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
###     
###     ell1 = 0
###     ell2 = 0
###     ELL = 0
###     ax[0].axvline(k1, color="black", ls="--")
###     color = ["royalblue", "forestgreen", "magenta"]
###     label = ["Tree", "Tree_NoWiggle", "Tree_BAO"]
### #    label = ["Tree"]
### #    label = ["Tree", "Tree_NoWiggle"]
###     for i_name, name in enumerate(label):
###         
###         k, bk = np.loadtxt("results_test/bk%d%d%d_%s_fft.dat" % (ell1, ell2, ELL, name), usecols=(1,2), unpack=True)
###         k_len = int(np.sqrt(len(k)))
###         k = k.reshape(k_len, k_len)
###         bk = bk.reshape(k_len, k_len)
###         kbin = k[0,:]
###         f_bk = interpolate.interp2d(kbin, kbin, bk, kind="cubic")
###         kbin = np.linspace(0.01, 0.3, 200)
###         bk = f_bk(k1, kbin)[:,0]
###         ax[0].plot(kbin, kbin*bk, "-", color=color[i_name], label="bk%d%d%d: %s" % (ell1, ell2, ELL, label[i_name]))
###     
###     ax[0].set_xlim(0.01, 0.21)
###     ax[0].set_xticks(np.linspace(0.02, 0.2, 10))
###     ax[0].set_xlabel("k2 [h/Mpc]")
###     ax[0].set_ylabel(r"$k2\, $ bispectrum(k1, k2)")
###     ax[0].legend(loc=0)
###     ax[0].set_title("k1 = %0.2f [h/Mpc]" % k1)
### #    plt.savefig("figure/bk%d%d%d_Tree_%02d.png" % (ell1, ell2, ELL, k1_i))
### #    os.system("display figure/bk%d%d%d_Tree_%02d.png" % (ell1, ell2, ELL, k1_i))
### #    plt.savefig("figure/bk%d%d%d_Tree_NoWiggle_%02d.png" % (ell1, ell2, ELL, k1_i))
### #    os.system("display figure/bk%d%d%d_Tree_NoWiggle_%02d.png" % (ell1, ell2, ELL, k1_i))
###     plt.savefig("figure/bk%d%d%d_Tree_BAO_%02d.png" % (ell1, ell2, ELL, k1_i))
###     os.system("display figure/bk%d%d%d_Tree_BAO_%02d.png" % (ell1, ell2, ELL, k1_i))

### ##########
### ## Compare the original bispectrum with the one computed from 3PCF.
### for k1_i, k1 in enumerate(np.linspace(0.02, 0.2, 10)):
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
###     
###     ell1 = 0
###     ell2 = 0
###     ELL = 0
###     ax[0].axvline(k1, color="black", ls="--")
###     color = ["royalblue"]
###     label = ["Tree"]
###     for i_name, name in enumerate(label):
###         
###         k, bk1 = np.loadtxt("results_test/bk%d%d%d_%s.dat" % (ell1, ell2, ELL, name), usecols=(1,2), unpack=True)
###         k, bk2 = np.loadtxt("results_test/bk%d%d%d_%s_from_3pcf.dat" % (ell1, ell2, ELL, name), usecols=(1,2), unpack=True)
###         k_len = int(np.sqrt(len(k)))
###         k = k.reshape(k_len, k_len)
###         bk1 = bk1.reshape(k_len, k_len)
###         bk2 = bk2.reshape(k_len, k_len)
###         k = k[k1_i, :]
###         bk1 = bk1[k1_i, :]
###         bk2 = bk2[k1_i, :]
###         ax[0].plot(k, abs(100.0*(bk1 - bk2)/bk1), "-", color=color[i_name], label="bk%d%d%d: %s" % (ell1, ell2, ELL, label[i_name]))
###     
###     ax[0].set_xlim(0.01, 0.21)
###     ax[0].set_xticks(np.linspace(0.02, 0.2, 10))
###     ax[0].set_xlabel("k [h/Mpc]")
###     ax[0].set_ylabel(r"diff. [$\%$]")
###     ax[0].legend(loc=0)
###     ax[0].set_title("k1 = %0.2f [h/Mpc]" % k1)
###     plt.savefig("figure/bk%d%d%d_Tree_%02d_diff.png" % (ell1, ell2, ELL, k1_i))
###     os.system("display figure/bk%d%d%d_Tree_%02d_diff.png" % (ell1, ell2, ELL, k1_i))
### 
### ##########
### ## Tree-level 3PCF
### ell1 = 0
### ell2 = 0
### ELL = 0
### f_zeta = {}
### for i_name, name in enumerate(["Tree", "Tree_NoWiggle", "Tree_BAO"]):
### #for i_name, name in enumerate(["Tree"]):
### #for i_name, name in enumerate(["Tree", "Tree_NoWiggle"]):
###     
###     r, zeta = np.loadtxt("results_test/zeta%d%d%d_%s_fft.dat" % (ell1, ell2, ELL, name), usecols=(1,2), unpack=True)
###     r_len = int(np.sqrt(len(r)))
###     r = r.reshape(r_len, r_len)
###     zeta = zeta.reshape(r_len, r_len)
###     rbin = r[0,:]
###     f_zeta_temp = interpolate.interp2d(rbin, rbin, zeta, kind="cubic")
###     rbin = np.linspace(30, 150, 13)
###     zeta = f_zeta_temp(rbin, rbin)
###     f_zeta.update({"%s" % name: interpolate.interp2d(rbin, rbin, zeta, kind="cubic")}) 
### 
### for r1_i, r1 in enumerate(np.linspace(30, 150, 13)):
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
###     
###     ax[0].axvline(r1, color="black", ls="--")
### 
###     for i_name, name in enumerate(["Tree", "Tree_NoWiggle", "Tree_BAO"]):
### #    for i_name, name in enumerate(["Tree"]):
### #    for i_name, name in enumerate(["Tree", "Tree_NoWiggle"]):
###         
###         color = ["royalblue", "forestgreen", "magenta"]
###         label = ["tree", "no-wiggle", "NL BAO"]
###         rbin = np.linspace(30, 150, 200)
###         zeta = f_zeta["%s" % name](r1, rbin)[:,0]
###         ax[0].plot(rbin, rbin**2*zeta, "-", color=color[i_name], label="zeta%d%d%d: %s" % (ell1, ell2, ELL, label[i_name]))
###     
###     ax[0].set_xlim(0.25, 155)
###     ax[0].set_xticks(np.linspace(30,150,13))
###     ax[0].set_xlabel("r2 [Mpc/h]")
###     ax[0].set_ylabel(r"$r2^2\, $ 3PCF(r1, r2)")
###     ax[0].legend(loc=0)
###     ax[0].set_title("r1 = %0.2f [h/Mpc]" % r1)
### #    plt.savefig("figure/zeta%d%d%d_Tree_%02d.png" % (ell1, ell2, ELL, r1_i))
### #    os.system("display figure/zeta%d%d%d_Tree_%02d.png" % (ell1, ell2, ELL, r1_i))
### #    plt.savefig("figure/zeta%d%d%d_Tree_NoWiggle_%02d.png" % (ell1, ell2, ELL, r1_i))
### #    os.system("display figure/zeta%d%d%d_Tree_NoWiggle_%02d.png" % (ell1, ell2, ELL, r1_i))
###     plt.savefig("figure/zeta%d%d%d_Tree_BAO_%02d.png" % (ell1, ell2, ELL, r1_i))
###     os.system("display figure/zeta%d%d%d_Tree_BAO_%02d.png" % (ell1, ell2, ELL, r1_i))
### 
###
### ##########
### ## Reconstructed tree-level bispectra with NL BAO
### for k1_i, k1 in enumerate(np.linspace(0.02, 0.2, 10)):
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
###     
###     ell1 = 0
###     ell2 = 0
###     ELL = 0
###     ax[0].axvline(k1, color="black", ls="--")
###     color = ["royalblue", "orange", "forestgreen", "tomato", "black"]
###     label = ["R = 5", "R = 10", "R = 15", "R = 20", "pre-recon"]
###     for i_R, R in enumerate([5, 10, 15, 20, "pre-recon"]):
###         
###         if R == "pre-recon":
###             k, bk = np.loadtxt("results_test/bk%d%d%d_Tree_BAO_fft.dat" % (ell1, ell2, ELL), usecols=(1,2), unpack=True)
###         else:
###             k, bk = np.loadtxt("results_test/bk%d%d%d_Tree_BAO_Reconstructed_R%02d_fft.dat" % (ell1, ell2, ELL, R), usecols=(1,2), unpack=True)
### 
###         k_len = int(np.sqrt(len(k)))
###         k = k.reshape(k_len, k_len)
###         bk = bk.reshape(k_len, k_len)
###         kbin = k[0,:]
###         f_bk = interpolate.interp2d(kbin, kbin, bk, kind="cubic")
###         kbin = np.linspace(0.01, 0.3, 200)
###         bk = f_bk(k1, kbin)[:,0]
###         ax[0].plot(kbin, kbin*bk, ls="-", color=color[i_R], label="bk%d%d%d: %s" % (ell1, ell2, ELL, label[i_R]))
###     
###     ax[0].set_xlim(0.01, 0.21)
###     ax[0].set_xticks(np.linspace(0.02, 0.2, 10))
###     ax[0].set_xlabel("k2 [h/Mpc]")
###     ax[0].set_ylabel(r"$k2\, $ bispectrum(k1, k2)")
###     ax[0].legend(loc=0)
###     ax[0].set_title("k1 = %0.2f [h/Mpc]" % k1)
###     plt.savefig("figure/bk%d%d%d_Tree_BAO_Reconstructed_%02d.png" % (ell1, ell2, ELL, k1_i))
###     os.system("display figure/bk%d%d%d_Tree_BAO_Reconstructed_%02d.png" % (ell1, ell2, ELL, k1_i))
### 
### 
### 
### ##########
### ## Reconstructed tree-level 3PCF with NL BAO
### ell1 = 0
### ell2 = 0
### ELL = 0
### f_zeta = {}
### label = ["R = 5", "R = 10", "R = 15", "R = 20", "pre-recon"]
### color = ["royalblue", "orange", "forestgreen", "tomato", "black"]
### for i_R, R in enumerate([5, 10, 15, 20, "pre-recon"]):
###     
###     if R == "pre-recon":
###         r, zeta = np.loadtxt("results_test/zeta%d%d%d_Tree_BAO_fft.dat" % (ell1, ell2, ELL), usecols=(1,2), unpack=True)
###     else:
###         r, zeta = np.loadtxt("results_test/zeta%d%d%d_Tree_BAO_Reconstructed_R%02d_fft.dat" % (ell1, ell2, ELL, R), usecols=(1,2), unpack=True)
###     
###     r_len = int(np.sqrt(len(r)))
###     r = r.reshape(r_len, r_len)
###     zeta = zeta.reshape(r_len, r_len)
###     rbin = r[0,:]
###     f_zeta_temp = interpolate.interp2d(rbin, rbin, zeta, kind="cubic")
###     rbin = np.linspace(30, 150, 13)
###     zeta = f_zeta_temp(rbin, rbin)
###     f_zeta.update({"%s" % label[i_R]: interpolate.interp2d(rbin, rbin, zeta, kind="cubic")}) 
### 
### for r1_i, r1 in enumerate(np.linspace(30, 150, 13)):
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
###     
###     ax[0].axvline(r1, color="black", ls="--")
### 
###     for i_R, R in enumerate([5, 10, 15, 20, "pre-recon"]):
###         
###         rbin = np.linspace(30, 150, 200)
###         zeta = f_zeta["%s" % label[i_R]](r1, rbin)[:,0]
###         ax[0].plot(rbin, rbin**2*zeta, "-", color=color[i_R], label="zeta%d%d%d: %s" % (ell1, ell2, ELL, label[i_R]))
###     
###     ax[0].set_xlim(0.25, 155)
###     ax[0].set_xticks(np.linspace(30,150,13))
###     ax[0].set_xlabel("r2 [Mpc/h]")
###     ax[0].set_ylabel(r"$r2^2\, $ 3PCF(r1, r2)")
###     ax[0].legend(loc=0)
###     ax[0].set_title("r1 = %0.2f [h/Mpc]" % r1)
### 
###     plt.savefig("figure/zeta%d%d%d_Tree_BAO_Reconstructed_%02d.png" % (ell1, ell2, ELL, r1_i))
###     os.system("display figure/zeta%d%d%d_Tree_BAO_Reconstructed_%02d.png" % (ell1, ell2, ELL, r1_i))
###
### ##########
### ## Decomposed tree-level bispectra with NL BAO
### 
### ell1 = 0
### ell2 = 0
### ELL = 0
### 
### label = []
### label.append("total")
### 
### label.append("FG_b3_f0")
### label.append("FS_b3_f0")
### label.append("FT_b3_f0")
### 
### label.append("FG_b2_f1")
### label.append("FS_b2_f1")
### label.append("FT_b2_f1")
### 
### label.append("FG_b1_f2")
### label.append("FS_b1_f2")
### label.append("FT_b1_f2")
### 
### label.append("GG_b2_f1")
### label.append("GS_b2_f1")
### label.append("GT_b2_f1")
### 
### label.append("GG_b1_f2")
### label.append("GS_b1_f2")
### label.append("GT_b1_f2")
### 
### label.append("GG_b0_f3")
### label.append("GS_b0_f3")
### label.append("GT_b0_f3")
### 
### label.append("b3_f1")
### label.append("b2_f2")
### label.append("b1_f3")
### label.append("b0_f4")
### 
### f_bk = {}
### for name in label:
###     
###     if name == "total":
###         k, bk = np.loadtxt("results_test/bk%d%d%d_Tree_BAO_fft.dat" % (ell1, ell2, ELL), usecols=(1,2), unpack=True)
###     else:
###         k, bk = np.loadtxt("results_test/bk%d%d%d_Tree_BAO_Template_%s_fft.dat" % (ell1, ell2, ELL, name), usecols=(1,2), unpack=True)
### 
###     k_len = int(np.sqrt(len(k)))
###     k = k.reshape(k_len, k_len)
###     bk = bk.reshape(k_len, k_len)
###     kbin = k[0,:]
###     f_bk.update({"%s" % name: interpolate.interp2d(kbin, kbin, bk, kind="cubic")})
### 
### kbin = np.linspace(0.01, 0.3, 200)
### 
### for k1_i, k1 in enumerate(np.linspace(0.02, 0.2, 10)):
### 
###     fig = plt.figure(figsize=(12.0,7.0))
### 
###     ##--------------------------
###     total_h = 0.94
###     total_w = 0.99
###     
###     dw = 0.4
###     dh = 0.00
###     left = 0.1
###     bottom = 0.11
###     height = (total_h - (bottom + dh))
###     width  = (total_w - (left + dw))
###     
###     #left bottom, width, height
###     rect1 = [left, bottom, width, height]
### 
###     ax = []
###     ax.append(fig.add_axes(rect1))
###     
###     ax[0].axvline(k1, color="black", ls="--")
### 
###     for i_name, name in enumerate(label[0:1]):
###         
###         bk = f_bk["%s" % name](k1, kbin)[:,0]
###         ax[0].plot(kbin, kbin*bk, "-", lw = 2.5, color="magenta", label="bk%d%d%d: %s" % (ell1, ell2, ELL, label[i_name]))
### 
###     for i_name, name in enumerate(label[1:10]):
###         
###         bk = f_bk["%s" % name](k1, kbin)[:,0]
###         ax[0].plot(kbin, kbin*bk, "-", label="bk%d%d%d: %s" % (ell1, ell2, ELL, label[i_name+1]))
### 
###     for i_name, name in enumerate(label[10:19]):
###         
###         bk = f_bk["%s" % name](k1, kbin)[:,0]
###         ax[0].plot(kbin, kbin*bk, "--", label="bk%d%d%d: %s" % (ell1, ell2, ELL, label[i_name+10]))
###  
###     for i_name, name in enumerate(label[19:23]):
###         
###         bk = f_bk["%s" % name](k1, kbin)[:,0]
###         ax[0].plot(kbin, kbin*bk, ":", label="bk%d%d%d: %s" % (ell1, ell2, ELL, label[i_name+19]))
###  
###     ax[0].set_xlim(0.01, 0.21)
###     ax[0].set_xticks(np.linspace(0.02, 0.2, 10))
###     ax[0].set_xlabel("k2 [h/Mpc]")
###     ax[0].set_ylabel(r"$k2\, $ bispectrum(k1, k2)")
###     ax[0].legend(loc="upper left", bbox_to_anchor=(1, 1), ncol=2)
###     ax[0].set_title("k1 = %0.2f [h/Mpc]" % k1)
###     plt.savefig("figure/bk%d%d%d_Tree_BAO_Decomposed_%02d.png" % (ell1, ell2, ELL, k1_i))
###     os.system("display figure/bk%d%d%d_Tree_BAO_Decomposed_%02d.png" % (ell1, ell2, ELL, k1_i))
### 
### ##########
### ## Decomposed tree-level 3PCF with NL BAO
### 
### ell1 = 0
### ell2 = 0
### ELL = 0
### 
### label = []
### label.append("total")
### 
### label.append("FG_b3_f0")
### label.append("FS_b3_f0")
### label.append("FT_b3_f0")
### 
### label.append("FG_b2_f1")
### label.append("FS_b2_f1")
### label.append("FT_b2_f1")
### 
### label.append("FG_b1_f2")
### label.append("FS_b1_f2")
### label.append("FT_b1_f2")
### 
### label.append("GG_b2_f1")
### label.append("GS_b2_f1")
### label.append("GT_b2_f1")
### 
### label.append("GG_b1_f2")
### label.append("GS_b1_f2")
### label.append("GT_b1_f2")
### 
### label.append("GG_b0_f3")
### label.append("GS_b0_f3")
### label.append("GT_b0_f3")
### 
### label.append("b3_f1")
### label.append("b2_f2")
### label.append("b1_f3")
### label.append("b0_f4")
### 
### f_zeta = {}
### for name in label:
###     
###     if name == "total":
###         r, zeta = np.loadtxt("results_test/zeta%d%d%d_Tree_BAO_fft.dat" % (ell1, ell2, ELL), usecols=(1,2), unpack=True)
###     else:
###         r, zeta = np.loadtxt("results_test/zeta%d%d%d_Tree_BAO_Template_%s_fft.dat" % (ell1, ell2, ELL, name), usecols=(1,2), unpack=True)
### 
###     r_len = int(np.sqrt(len(r)))
###     r = r.reshape(r_len, r_len)
###     zeta = zeta.reshape(r_len, r_len)
###     rbin = r[0,:]
###     f_zeta_temp = interpolate.interp2d(rbin, rbin, zeta, kind="cubic")
###     rbin = np.linspace(20, 160, 15)
###     zeta = f_zeta_temp(rbin, rbin)
###     f_zeta.update({"%s" % name: interpolate.interp2d(rbin, rbin, zeta, kind="cubic")}) 
### 
### rbin = np.linspace(20, 160, 200)
### 
### for r1_i, r1 in enumerate(np.linspace(30, 150, 13)):
### 
###     fig = plt.figure(figsize=(12.0,7.0))
### 
###     ##--------------------------
###     total_h = 0.94
###     total_w = 0.99
###     
###     dw = 0.4
###     dh = 0.00
###     left = 0.1
###     bottom = 0.11
###     height = (total_h - (bottom + dh))
###     width  = (total_w - (left + dw))
###     
###     #left bottom, width, height
###     rect1 = [left, bottom, width, height]
### 
###     ax = []
###     ax.append(fig.add_axes(rect1))
###     
###     ax[0].axvline(r1, color="black", ls="--")
### 
###     for i_name, name in enumerate(label[0:1]):
###         
###         zeta = f_zeta["%s" % name](r1, rbin)[:,0]
###         ax[0].plot(rbin, rbin**2*zeta, "-", lw=2.5, color="magenta", label="zeta%d%d%d: %s" % (ell1, ell2, ELL, label[i_name]))
### 
###     for i_name, name in enumerate(label[1:10]):
###         
###         zeta = f_zeta["%s" % name](r1, rbin)[:,0]
###         ax[0].plot(rbin, rbin**2*zeta, "-", label="zeta%d%d%d: %s" % (ell1, ell2, ELL, label[i_name+1]))
### 
###     for i_name, name in enumerate(label[10:19]):
###         
###         zeta = f_zeta["%s" % name](r1, rbin)[:,0]
###         ax[0].plot(rbin, rbin**2*zeta, "--", label="zeta%d%d%d: %s" % (ell1, ell2, ELL, label[i_name+10]))
### 
###     for i_name, name in enumerate(label[19:23]):
###         
###         zeta = f_zeta["%s" % name](r1, rbin)[:,0]
###         ax[0].plot(rbin, rbin**2*zeta, ":", label="zeta%d%d%d: %s" % (ell1, ell2, ELL, label[i_name+19]))
### 
###     ax[0].set_xlim(20, 160)
###     ax[0].set_xticks(np.linspace(20, 160, 15))
###     ax[0].set_xlabel("r2 [h/Mpc]")
###     ax[0].set_ylabel(r"$r2^2\, $ 3PCF(r1, r2)")
###     ax[0].legend(loc="upper left", bbox_to_anchor=(1, 1), ncol=2)
###     ax[0].set_title("r1 = %2.1f [Mpc/h]" % r1)
###     plt.savefig("figure/zeta%d%d%d_Tree_BAO_Decomposed_%02d.png" % (ell1, ell2, ELL, r1_i))
###     os.system("display figure/zeta%d%d%d_Tree_BAO_Decomposed_%02d.png" % (ell1, ell2, ELL, r1_i))
### 
### ##########
### ## Tree-level bispectra generated by primordial non-Gaussianities
### for k1_i, k1 in enumerate(np.linspace(0.02, 0.2, 10)):
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
###     
###     ell1 = 0
###     ell2 = 0
###     ELL = 0
###     ax[0].axvline(k1, color="black", ls="--")
###     color = ["royalblue", "forestgreen", "tomato"]
###     label = ["Tree_NonGaussian_Local", "Tree_NonGaussian_Equilateral", "Tree_NonGaussian_Orthogonal"]
###     for i_name, name in enumerate(label):
###         
###         k, bk = np.loadtxt("results_test/bk%d%d%d_%s_fft.dat" % (ell1, ell2, ELL, name), usecols=(1,2), unpack=True)
### 
###         k_len = int(np.sqrt(len(k)))
###         k = k.reshape(k_len, k_len)
###         bk = bk.reshape(k_len, k_len)
###         kbin = k[0,:]
###         f_bk = interpolate.interp2d(kbin, kbin, bk, kind="cubic")
###         kbin = np.linspace(0.01, 0.3, 200)
###         bk = f_bk(k1, kbin)[:,0]
###         ax[0].plot(kbin, kbin*bk, "-", color=color[i_name], label="bk%d%d%d: %s" % (ell1, ell2, ELL, label[i_name]))
###     
###     ax[0].set_xlim(0.01, 0.21)
###     ax[0].set_xticks(np.linspace(0.02, 0.2, 10))
###     ax[0].set_xlabel("k2 [h/Mpc]")
###     ax[0].set_ylabel(r"$k2\, $ bispectrum(k1, k2)")
###     ax[0].legend(loc=0)
###     ax[0].set_title("k1 = %0.2f [h/Mpc]" % k1)
###     plt.savefig("figure/bk%d%d%d_Tree_PNG_%02d.png" % (ell1, ell2, ELL, k1_i))
###     os.system("display figure/bk%d%d%d_Tree_PNG_%02d.png" % (ell1, ell2, ELL, k1_i))
### 
### ##########
### ## Tree-level 3PCF generated by PNGs
### ell1 = 0
### ell2 = 0
### ELL = 0
### f_zeta = {}
### for i_name, name in enumerate(["Tree_NonGaussian_Local", "Tree_NonGaussian_Equilateral", "Tree_NonGaussian_Orthogonal"]):
###     
###     r, zeta = np.loadtxt("results_test/zeta%d%d%d_%s_fft.dat" % (ell1, ell2, ELL, name), usecols=(1,2), unpack=True)
###     r_len = int(np.sqrt(len(r)))
###     r = r.reshape(r_len, r_len)
###     zeta = zeta.reshape(r_len, r_len)
###     rbin = r[0,:]
###     f_zeta_temp = interpolate.interp2d(rbin, rbin, zeta, kind="cubic")
###     rbin = np.linspace(30, 150, 13)
###     zeta = f_zeta_temp(rbin, rbin)
###     f_zeta.update({"%s" % name: interpolate.interp2d(rbin, rbin, zeta, kind="cubic")}) 
### 
### for r1_i, r1 in enumerate(np.linspace(30, 150, 13)):
### 
###     fig = plt.figure(figsize=(7.0,7.0))
###     ax = []
###     ax.append(fig.add_axes(rect1))
###     
###     ax[0].axvline(r1, color="black", ls="--")
### 
###     for i_name, name in enumerate(["Tree_NonGaussian_Local", "Tree_NonGaussian_Equilateral", "Tree_NonGaussian_Orthogonal"]):
###         color = ["royalblue", "forestgreen", "tomato"]
###         label = ["Tree_NonGaussian_Local", "Tree_NonGaussian_Equilateral", "Tree_NonGaussian_Orthogonal"]
###         rbin = np.linspace(30, 150, 200)
###         zeta = f_zeta["%s" % name](r1, rbin)[:,0]
###         ax[0].plot(rbin, rbin**2*zeta, "-", color=color[i_name], label="zeta%d%d%d: %s" % (ell1, ell2, ELL, label[i_name]))
###         
###     ax[0].set_xlim(0.25, 155)
###     ax[0].set_xticks(np.linspace(30,150,13))
###     ax[0].set_xlabel("r2 [Mpc/h]")
###     ax[0].set_ylabel(r"$r2^2\, $ 3PCF(r1, r2)")
###     ax[0].legend(loc=0)
###     ax[0].set_title("r1 = %0.2f [h/Mpc]" % r1)
###     plt.savefig("figure/zeta%d%d%d_Tree_PNG_%02d.png" % (ell1, ell2, ELL, r1_i))
###     os.system("display figure/zeta%d%d%d_Tree_PNG_%02d.png" % (ell1, ell2, ELL, r1_i))
###   


