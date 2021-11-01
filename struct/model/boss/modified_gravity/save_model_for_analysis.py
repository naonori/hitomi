########################
# Purpose of this code #
########################
# The code re-stores the theoretically computed 2PCF and 3PCF results in a form that facilitates parameter estimation in "Montepython". 
# (see also $WORK/covariance/boss/modified_gravity/save_data_and_cov_for_analysis.py)
# 
# In this study of the modified gravity theory, we consider seven cases of data combination:
# (1) xi0+xi2;
# (2) xi0+xi2+zeta000+zeta110;
# (3) xi0+xi2+zeta202+zeta112;
# (4) xi0+xi2+zeta000+zeta202;
# (5) xi0+xi2+zeta000+zeta202+zeta110;
# (6) xi0+xi2+zeta000+zeta202+zeta112;
# (7) xi0+xi2+zeta000+zeta202+zeta110+zeta112.
#
# In addition, we test the following four cases:
#
# TEST (1) xi0+xi2+zeta000;
# TEST (2) xi0+xi2+zeta202;
# TEST (3) xi0+xi2+zeta110;
# TEST (4) xi0+xi2+zeta112;
#
# We only focus on the large scale of 80<=r<=150 [Mpc/h] for parameter estimation in "montepython".
#
######################################################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import numpy as np
from distutils.util import strtobool
import pickle

NS = "North"
zbin = 1
recon = False

for NS in ["North", "South"]:
    for zbin in [1,3]:
        for recon in [False, True]:

            # Make an output directory #
            OUTPUT = "model_for_analysis"
            try:
                os.mkdir(OUTPUT)
            except:
                pass
            
            MODEL_case1 = {}
            MODEL_case2 = {}
            MODEL_case3 = {}
            MODEL_case4 = {}
            MODEL_case5 = {}
            MODEL_case6 = {}
            MODEL_case7 = {}
            
            MODEL_test1 = {}
            MODEL_test2 = {}
            MODEL_test3 = {}
            MODEL_test4 = {}
            
            rbin_2pcf = np.linspace(30,150,25)
            rbin_3pcf = np.linspace(30,150,13)
            rbin_min_2pcf = 10
            rbin_min_3pcf = 5
            rbin_max_2pcf = 25
            rbin_max_3pcf = 13
            
            for data in ["galaxy", "mock"]:
            
                # INPUT directory #
                INPUT = "model_%s_%s_zbin%d" % (data, NS, zbin)
                
                PARAM_NAME_2pcf = []
                PARAM_NAME_2pcf.append("b1_b1")
                PARAM_NAME_2pcf.append("b1_f")
                PARAM_NAME_2pcf.append("f_f")
                
                for param_name in PARAM_NAME_2pcf:
                 
                    if recon == False:
                        xi0 = np.loadtxt("%s/xi0_Tree_BAO_Template_%s.dat" % (INPUT, param_name), usecols=(1,), unpack=True)
                        xi2 = np.loadtxt("%s/xi2_Tree_BAO_Template_%s.dat" % (INPUT, param_name), usecols=(1,), unpack=True)
                    elif recon == True:
                        xi0 = np.loadtxt("%s/xi0_Tree_BAO_Template_%s_recon.dat" % (INPUT, param_name), usecols=(1,), unpack=True)
                        xi2 = np.loadtxt("%s/xi2_Tree_BAO_Template_%s_recon.dat" % (INPUT, param_name), usecols=(1,), unpack=True)
                
                    xi0 = xi0[rbin_min_2pcf:rbin_max_2pcf]
                    xi2 = xi2[rbin_min_2pcf:rbin_max_2pcf]
                
                    xi_M = np.hstack((xi0, xi2))
                
                    MODEL_case1.update({"xi0xi2_%s_%s" % (param_name, data): xi_M})
                    MODEL_case2.update({"xi0xi2_%s_%s" % (param_name, data): xi_M})
                    MODEL_case3.update({"xi0xi2_%s_%s" % (param_name, data): xi_M})
                    MODEL_case4.update({"xi0xi2_%s_%s" % (param_name, data): xi_M})
                    MODEL_case5.update({"xi0xi2_%s_%s" % (param_name, data): xi_M})
                    MODEL_case6.update({"xi0xi2_%s_%s" % (param_name, data): xi_M})
                    MODEL_case7.update({"xi0xi2_%s_%s" % (param_name, data): xi_M})
                                                                   
                    MODEL_test1.update({"xi0xi2_%s_%s" % (param_name, data): xi_M})
                    MODEL_test2.update({"xi0xi2_%s_%s" % (param_name, data): xi_M})
                    MODEL_test3.update({"xi0xi2_%s_%s" % (param_name, data): xi_M})
                    MODEL_test4.update({"xi0xi2_%s_%s" % (param_name, data): xi_M})
                
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
                
                for param_name in PARAM_NAME_3pcf:
                    
                    zeta000 = np.loadtxt("%s/zeta000_Tree_BAO_Template_%s.dat" % (INPUT, param_name), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                    zeta110 = np.loadtxt("%s/zeta110_Tree_BAO_Template_%s.dat" % (INPUT, param_name), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                    zeta202 = np.loadtxt("%s/zeta202_Tree_BAO_Template_%s.dat" % (INPUT, param_name), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                    zeta112 = np.loadtxt("%s/zeta112_Tree_BAO_Template_%s.dat" % (INPUT, param_name), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                
                    # for zeta000, zeta110 and zeta112, r1>=r2. 
                    zeta000 = np.array([zeta000[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                    zeta202 = np.array([zeta202[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf)])
                    zeta110 = np.array([zeta110[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                    zeta112 = np.array([zeta112[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                
                    zeta_M = np.hstack((zeta000, zeta110))
                    MODEL_case2.update({"zeta000zeta110_%s_%s" % (param_name, data): zeta_M})
                
                    zeta_M = np.hstack((zeta202, zeta112))
                    MODEL_case3.update({"zeta202zeta112_%s_%s" % (param_name, data): zeta_M})
                
                    zeta_M = np.hstack((zeta000, zeta202))
                    MODEL_case4.update({"zeta000zeta202_%s_%s" % (param_name, data): zeta_M})
                
                    zeta_M = np.hstack((zeta000, zeta202, zeta110))
                    MODEL_case5.update({"zeta000zeta202zeta110_%s_%s" % (param_name, data): zeta_M})
                
                    zeta_M = np.hstack((zeta000, zeta202, zeta112))
                    MODEL_case6.update({"zeta000zeta202zeta112_%s_%s" % (param_name, data): zeta_M})
                
                    zeta_M = np.hstack((zeta000, zeta202, zeta110, zeta112))
                    MODEL_case7.update({"zeta000zeta202zeta110zeta112_%s_%s" % (param_name, data): zeta_M})
                
                    ######
                    MODEL_test1.update({"zeta000_%s_%s" % (param_name, data): zeta000})
                    MODEL_test2.update({"zeta202_%s_%s" % (param_name, data): zeta202})
                    MODEL_test3.update({"zeta110_%s_%s" % (param_name, data): zeta110})
                    MODEL_test4.update({"zeta112_%s_%s" % (param_name, data): zeta112})
                
                
            # save data
            if recon == False:
            
                with open("%s/MODEL_%s_zbin%d_case1.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_case1, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_case2.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_case2, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_case3.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_case3, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_case4.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_case4, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_case5.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_case5, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_case6.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_case6, pp)               
                
                with open("%s/MODEL_%s_zbin%d_case7.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_case7, pp)               
                                                               
                with open("%s/MODEL_%s_zbin%d_test1.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_test1, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_test2.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_test2, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_test3.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_test3, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_test4.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_test4, pp)
                 
            elif recon == True:
            
                with open("%s/MODEL_%s_zbin%d_recon_case1.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_case1, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_recon_case2.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_case2, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_recon_case3.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_case3, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_recon_case4.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_case4, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_recon_case5.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_case5, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_recon_case6.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_case6, pp)               
                
                with open("%s/MODEL_%s_zbin%d_recon_case7.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_case7, pp)               
                                                               
                with open("%s/MODEL_%s_zbin%d_recon_test1.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_test1, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_recon_test2.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_test2, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_recon_test3.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_test3, pp)               
                                                              
                with open("%s/MODEL_%s_zbin%d_recon_test4.p" % (OUTPUT, NS, zbin), "wb") as pp:
                    pickle.dump(MODEL_test4, pp)
                 
            
            
