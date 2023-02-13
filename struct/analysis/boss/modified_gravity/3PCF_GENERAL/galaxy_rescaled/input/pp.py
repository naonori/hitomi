#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


for (NS, zbin) in [("North", 1), ("South", 1), ("North", 3), ("South", 3)]:

    print(NS,  zbin)

    dir_name = "/mwork2/sugiymnn/WORK/analysis/boss/modified_gravity/fisher_analysis/results_fisher/"
    fname_fid = "fiducial_values_%s_zbin%d.dat" % (NS, zbin)
    fname_fisher = "fisher_GENERAL_%s_zbin%d_case7.dat" % (NS, zbin)
    
    param_fid = np.loadtxt(dir_name+fname_fid)
    param_fisher = np.loadtxt(dir_name+fname_fisher)[5,:]
    
    print("[%.2f, %.2f, %.2f, %.2f]" % (param_fid[0], param_fid[0]-5.0 * param_fisher[1], param_fid[0]+5.0 * param_fisher[1], param_fisher[1]))
    print("[%.2f, %.2f, %.2f, %.2f]" % (param_fid[1], param_fid[1]-5.0 * param_fisher[2], param_fid[1]+5.0 * param_fisher[2], param_fisher[2]))
    print("[%.2f, %.2f, %.2f, %.2f]" % (param_fid[2], param_fid[2]-5.0 * param_fisher[3], param_fid[2]+5.0 * param_fisher[3], param_fisher[3]))
    print("[%.2f, %.2f, %.2f, %.2f]" % (param_fid[3], param_fid[3]-5.0 * param_fisher[4], param_fid[3]+5.0 * param_fisher[4], param_fisher[4]))
    print("[%.2f, %.2f, %.2f, %.2f]" % (param_fid[4], param_fid[4]-5.0 * param_fisher[5], param_fid[4]+5.0 * param_fisher[5], param_fisher[5]))
    
    print("[%.2f, %.2f, %.2f, %.2f]" % ((13.0/21.0)*param_fid[3], (13.0/21.0)*param_fid[3] - 5.0 * param_fisher[6],\
                                        (13.0/21.0)*param_fid[3] + 5.0 * param_fisher[6], param_fisher[6]))
    
    print("[%.2f, %.2f, %.2f, %.2f]" % (1.0, 1.0-5.0 * param_fisher[7], 1.0+5.0 * param_fisher[7], param_fisher[7]))
    
    print("[%.2f, %.2f, %.2f, %.2f]" % (param_fid[6], param_fid[6]-5.0 * param_fisher[8], param_fid[6]+5.0 * param_fisher[8], param_fisher[8]))



