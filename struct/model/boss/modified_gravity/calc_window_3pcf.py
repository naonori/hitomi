########################
# Purpose of this code #
########################
# If we ignore the AP effect, we can precompute the window function correction to the 3PCF.
# The window function correction is applied to each of the terms in the decomposed 3PCF.
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
            print(data, NS, zbin)

            rbin = np.linspace(30.0, 150.0, 13)
            r_len = 13
            
            # H factor 
            H = lambda x, y, z: float(wigner_3j(x,y,z,0,0,0))
            
            # make an output directory #
            OUTPUT = "model_%s_%s_zbin%d" % (data, NS, zbin)
            try:
                os.mkdir(OUTPUT)
            except:
                pass
            
            # reading window 3PCFs measured from data
            WORK = "/mwork0/sugiymnn/WORK"
            
            if data == "galaxy":
                MEASUREMENT_3PCF = WORK + "/measurement/boss/3pcf"
            elif data == "mock":
                MEASUREMENT_3PCF = WORK + "/measurement/boss/3pcf_mock"
            
            DATA = "%s_%s_zbin%d_Window" % (data, NS, zbin)
            
            MULTIPOLE = []
            MONO = [(0,0,0),(1,1,0),(2,2,0),(3,3,0),(4,4,0)]
            QUAD = [(2,0,2),(1,1,2),(0,2,2),(3,1,2),(2,2,2),(1,3,2),(4,2,2),(3,3,2),(2,4,2)]
            MULTIPOLE.extend(MONO)
            MULTIPOLE.extend(QUAD)
            
            w_3pcf = {}
            for (ell1, ell2, ELL) in MULTIPOLE:
                w_3pcf_temp = np.zeros(0)
                rbin1 = np.zeros(0)
                rbin2 = np.zeros(0)
                for r1_i in range(r_len):
            
                    fname000 = "%s/%s/zeta000_%02d_window.dat" % (MEASUREMENT_3PCF, DATA, r1_i)
                    w000 = np.loadtxt(fname000, usecols=(2,), unpack=True)
                    
                    fname = "%s/%s/zeta%d%d%d_%02d_window.dat" % (MEASUREMENT_3PCF, DATA, ell1, ell2, ELL, r1_i)
                    r1, r2, w = np.loadtxt(fname, usecols=(0,1,2), unpack=True)
            
                    w = w / w000
                
                    rbin1 = np.hstack((rbin1, r1))
                    rbin2 = np.hstack((rbin2, r2))
                    w_3pcf_temp = np.hstack((w_3pcf_temp, w))
            
                w_3pcf.update({"zeta%d%d%d" % (ell1, ell2, ELL) : w_3pcf_temp})
            
            # reading parameter-decomposed 3PCFs calculated using perturbation theory
            PARAM_NAME_3pcf = []
            PARAM_NAME_3pcf.append("FG_b3_f0")
            PARAM_NAME_3pcf.append("FS_b3_f0")
            PARAM_NAME_3pcf.append("FT_b3_f0")
            PARAM_NAME_3pcf.append("FG_b2_f1")
            PARAM_NAME_3pcf.append("FS_b2_f1")
            PARAM_NAME_3pcf.append("FT_b2_f1")
            PARAM_NAME_3pcf.append("FG_b1_f2")
            PARAM_NAME_3pcf.append("FS_b1_f2")
            PARAM_NAME_3pcf.append("FT_b1_f2")
            
            PARAM_NAME_3pcf.append("GG_b2_f1")
            PARAM_NAME_3pcf.append("GS_b2_f1")
            PARAM_NAME_3pcf.append("GT_b2_f1")
            PARAM_NAME_3pcf.append("GG_b1_f2")
            PARAM_NAME_3pcf.append("GS_b1_f2")
            PARAM_NAME_3pcf.append("GT_b1_f2")
            PARAM_NAME_3pcf.append("GG_b0_f3")
            PARAM_NAME_3pcf.append("GS_b0_f3")
            PARAM_NAME_3pcf.append("GT_b0_f3")
            
            PARAM_NAME_3pcf.append("b3_f1")
            PARAM_NAME_3pcf.append("b2_f2")
            PARAM_NAME_3pcf.append("b1_f3")
            PARAM_NAME_3pcf.append("b0_f4")
            
            zeta = {}
            for param_name in PARAM_NAME_3pcf:
                for (ell1, ell2, ELL) in MULTIPOLE:
                    fname = "model_%s_zbin%d/zeta%d%d%d_Tree_BAO_Template_%s.dat" % (data, zbin, ell1, ell2, ELL, param_name)
                    zeta.update({"zeta%d%d%d_%s" % (ell1, ell2, ELL, param_name) : np.loadtxt(fname, usecols=(2,), unpack=True)})
            
            
            
            # window function corrections 
            for param_name in PARAM_NAME_3pcf:
                for (ell1, ell2, ELL) in [(0,0,0), (1,1,0), (2,0,2), (1,1,2)]:
                    zeta_w = np.zeros(r_len**2)
                    for (ell1_d, ell2_d, ELL_d) in MULTIPOLE:
                        for (ell1_dd, ell2_dd, ELL_dd) in MULTIPOLE:
            
                            factor_W9 = wigner_9j(ell1_dd, ell2_dd, ELL_dd, ell1_d, ell2_d, ELL_d, ell1, ell2, ELL)
                            factor_H = H(ell1, ell2, ELL) * H(ell1, ell1_d, ell1_dd) * H(ell2, ell2_d, ell2_dd) * H(ELL, ELL_d, ELL_dd)\
                                     / ( H(ell1_d, ell2_d, ELL_d) * H(ell1_dd, ell2_dd, ELL_dd) )
                            factor_N = (2 * ell1 + 1) * (2 * ell2 + 1) * (2 * ELL + 1)
                            factor = float(factor_W9) * float(factor_H) * float(factor_N)
            
                            zeta_w = zeta_w + factor * w_3pcf["zeta%d%d%d" % (ell1_dd, ell2_dd, ELL_dd)] * zeta["zeta%d%d%d_%s" % (ell1_d, ell2_d, ELL_d, param_name)]
                        
                    fname = "%s/zeta%d%d%d_Tree_BAO_Template_%s.dat" % (OUTPUT, ell1, ell2, ELL, param_name)
                    np.savetxt(fname, np.array([rbin1, rbin2, zeta_w]).T, fmt="%.5f \t %.5f \t %.7e")
            
