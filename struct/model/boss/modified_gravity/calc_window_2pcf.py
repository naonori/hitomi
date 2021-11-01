########################
# Purpose of this code #
########################
# If we ignore the AP effect, we can precompute the window function correction to the 2PCF.
# The window function correction is applied to each of the terms in the decomposed 2PCF.
######################################################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
import scipy as sp
from sympy.physics.wigner import wigner_9j, wigner_3j
import pickle

for data in ["galaxy", "mock"]:
    for NS in ["North", "South"]:
        for zbin in [1,3]:
            for recon in [False, True]:

                # H factor 
                H = lambda x, y, z: float(wigner_3j(x,y,z,0,0,0))
                
                # make an output directory #
                OUTPUT = "model_%s_%s_zbin%d" % (data, NS, zbin)
                try:
                    os.mkdir(OUTPUT)
                except:
                    pass
                
                # reading window 2PCFs measured from data
                WORK = "/mwork0/sugiymnn/WORK"
                
                if data == "galaxy":
                    MEASUREMENT_2PCF = WORK + "/measurement/boss/2pcf"
                elif data == "mock":
                    MEASUREMENT_2PCF = WORK + "/measurement/boss/2pcf_mock"
                
                if recon == False:
                    DATA = "%s_%s_zbin%d_Window" % (data, NS, zbin)
                elif recon == True:
                    DATA = "%s_%s_zbin%d_recon_R15_Window" % (data, NS, zbin)
                
                w_2pcf = {}
                for ELL in [0,2,4]:
                    if recon == False:
                        fname = "%s/%s/xi%d_window.dat" % (MEASUREMENT_2PCF, DATA, ELL)
                    elif recon == True:
                        fname = "%s/%s/xi%d_recon_window.dat" % (MEASUREMENT_2PCF, DATA, ELL)
                
                    w_2pcf.update({"xi%d" % ELL: np.loadtxt(fname, usecols=(1,), unpack=True)})
                
                w_2pcf["xi2"] = w_2pcf["xi2"] / w_2pcf["xi0"]
                w_2pcf["xi4"] = w_2pcf["xi4"] / w_2pcf["xi0"]
                w_2pcf["xi0"] = w_2pcf["xi0"] / w_2pcf["xi0"]
                
                # reading parameter-decomposed 2PCFs calculated using perturbation theory
                PARAM_NAME_2pcf = []
                PARAM_NAME_2pcf.append("b1_b1")
                PARAM_NAME_2pcf.append("b1_f")
                PARAM_NAME_2pcf.append("f_f")
                
                xi = {}
                for param_name in PARAM_NAME_2pcf:
                    for ELL in [0,2,4]:
                        if recon == False:
                            fname = "model_%s_zbin%d/xi%d_Tree_BAO_Template_%s.dat" % (data, zbin, ELL, param_name)
                        elif recon == True:
                            fname = "model_%s_zbin%d/xi%d_Tree_BAO_Template_%s_recon.dat" % (data, zbin, ELL, param_name)
                
                        xi.update({"xi%d_%s" % (ELL, param_name) : np.loadtxt(fname, usecols=(1,), unpack=True)})


                # window function corrections 
                rbin = np.linspace(30.0, 150.0, 25)
                r_len = 25
                for param_name in PARAM_NAME_2pcf:
                    for ELL in [0,2,4]:
                        xi_w = np.zeros(r_len)
                        for ELL1 in [0,2,4]:
                            for ELL2 in [0,2,4]:
                                factor = (2.0 * float(ELL) + 1.0) * H(ELL1, ELL2, ELL)**2
                                xi_w = xi_w + factor * w_2pcf["xi%d" % ELL1] * xi["xi%d_%s" % (ELL2, param_name)]
                        
                
                        if recon == False:
                            fname = "%s/xi%d_Tree_BAO_Template_%s.dat" % (OUTPUT, ELL, param_name)
                        elif recon == True:
                            fname = "%s/xi%d_Tree_BAO_Template_%s_recon.dat" % (OUTPUT, ELL, param_name)
                
                        np.savetxt(fname, np.array([rbin, xi_w]).T, fmt="%.5f \t %.7e" )
                
                
                
                
                
