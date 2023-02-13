########################
# Purpose of this code #
########################
# This code saves the data for six minimum scales of rmin = 30, 40, 50, 60, 70, and 80 Mpc/h, while keeping the same settings as "save_data_and_cov_for_analysis.py". This data can be used in the future when a 3PCF model is created that is applicable down to small scales. It can also be used to simulate the behavior of small scales using Fisher analysis.
######################################################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import numpy as np

for rbin_min_3pcf in [0,1,2,3,4,5]:
    for NS in ["North", "South"]:
        for zbin in [1,3]:
            for recon in [False, True]:
    
                print(NS, "zbin = %d" % zbin, "recon = %s" % recon)
    
                DATA_case1 = {}
                DATA_case2 = {}
                DATA_case3 = {}
                DATA_case4 = {}
                DATA_case5 = {}
                DATA_case6 = {}
                DATA_case7 = {}
                DATA_case8 = {}
                
                DATA_test1 = {}
                DATA_test2 = {}
                DATA_test3 = {}
                DATA_test4 = {}
                
                # INPUT directory #
                if recon == False:
                    INPUT = "data_and_cov_%s_zbin%d" % (NS, zbin)
                elif recon == True:
                    INPUT = "data_and_cov_%s_zbin%d_recon_R15" % (NS, zbin)
                
                # r-bin for the 2pcf and 3pcf. #
                # Use the range 80<= r <= 150 #
                rbin_2pcf = np.linspace(30,150,25)
                rbin_3pcf = np.linspace(30,150,13)
                rbin_min_2pcf = 2*rbin_min_3pcf
                rbin_max_2pcf = 25
                rbin_max_3pcf = 13
                
                # Covariance matrixes of the 2PCF #
                cov_xi_xi = {}
                for ell in [0,2,4]:
                    for ell_d in [0,2,4]:
                        cov_temp = np.loadtxt("%s/cov_xi%d_xi%d.dat" % (INPUT, ell, ell_d))
                        cov_temp = cov_temp[rbin_min_2pcf:rbin_max_2pcf, rbin_min_2pcf:rbin_max_2pcf]
                        cov_xi_xi.update({"cov_xi%d_xi%d" % (ell, ell_d): cov_temp})
                
                # Cross-covariance matrixes between the 2PCF and 3PCF #
                cov_xi_zeta = {}
                for ell in [0,2,4]:
                    for (ell1_d, ell2_d, ELL_d) in [(0,0,0),(1,1,0),(2,0,2),(1,1,2)]:
                        cov_temp = np.loadtxt("%s/cov_xi%d_zeta%d%d%d.dat" % (INPUT, ell, ell1_d, ell2_d, ELL_d))
                        cov_temp = cov_temp[rbin_min_2pcf:rbin_max_2pcf, :]
                        cov_temp = cov_temp.reshape(rbin_max_2pcf - rbin_min_2pcf, rbin_max_3pcf, rbin_max_3pcf)
                        # IMPORTANT:  r1 >= r2 for zeta000, zeta110,and  zeta112. #
                        if ell1_d == ell2_d:
                            cov_temp = np.array([cov_temp[:, i, j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j ]).T
                        else:
                            cov_temp = np.array([cov_temp[:, i, j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf)]).T
                        
                        cov_xi_zeta.update({"cov_xi%d_zeta%d%d%d" % (ell, ell1_d, ell2_d, ELL_d): cov_temp})
                
                
                # Cross-covariance matrixes between the 2PCF and 3PCF #
                cov_zeta_xi = {}
                for (ell1, ell2, ELL) in [(0,0,0),(1,1,0),(2,0,2),(1,1,2)]:
                    for ell_d in [0,2,4]:
                        cov_temp = np.loadtxt("%s/cov_zeta%d%d%d_xi%d.dat" % (INPUT, ell1, ell2, ELL, ell_d))
                        cov_temp = cov_temp[:, rbin_min_2pcf:rbin_max_2pcf]
                        cov_temp = cov_temp.reshape(rbin_max_3pcf, rbin_max_3pcf, rbin_max_2pcf - rbin_min_2pcf)
                        # IMPORTANT:  r1 >= r2 for zeta000, zeta110,and  zeta112. #
                        if ell1 == ell2:
                            cov_temp = np.array([cov_temp[i, j, :] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j ])
                        else:
                            cov_temp = np.array([cov_temp[i, j, :] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf)])
                        
                        cov_zeta_xi.update({"cov_zeta%d%d%d_xi%d" % (ell1, ell2, ELL, ell_d): cov_temp})
                
                
                # Covariance matrixes of the 3PCF #
                cov_zeta_zeta = {}
                for (ell1, ell2, ELL) in [(0,0,0),(1,1,0),(2,0,2),(1,1,2)]:
                    for (ell1_d, ell2_d, ELL_d) in [(0,0,0),(1,1,0),(2,0,2),(1,1,2)]:
                        cov_temp = np.loadtxt("%s/cov_zeta%d%d%d_zeta%d%d%d.dat" % (INPUT, ell1, ell2, ELL, ell1_d, ell2_d, ELL_d))
                        cov_temp = cov_temp.reshape(rbin_max_3pcf, rbin_max_3pcf, rbin_max_3pcf, rbin_max_3pcf)
                        # IMPORTANT:  r1 >= r2 for zeta000, zeta110,and  zeta112. #
                        if ell1 == ell2:
                            cov_temp = np.array([cov_temp[i, j, :, :] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j ])
                        else:
                            cov_temp = np.array([cov_temp[i, j, :, :] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf)])
                
                        if ell1_d == ell2_d:
                            cov_temp = np.array([cov_temp[:, i, j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j ]).T
                        else:
                            cov_temp = np.array([cov_temp[:, i, j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf)]).T
                
                        cov_zeta_zeta.update({"cov_zeta%d%d%d_zeta%d%d%d" % (ell1, ell2, ELL, ell1_d, ell2_d, ELL_d): cov_temp})
                
                
                # Covariance matrix for case 1 #
                C0 = np.hstack((cov_xi_xi["cov_xi0_xi0"], cov_xi_xi["cov_xi0_xi2"]))
                C1 = np.hstack((cov_xi_xi["cov_xi2_xi0"], cov_xi_xi["cov_xi2_xi2"]))
                C = np.vstack((C0,C1))
                DATA_case1.update({"xi0xi2_cov": C})
                
                # Covariance matrix for case 2 #
                C0 = np.hstack((cov_xi_xi["cov_xi0_xi0"], cov_xi_xi["cov_xi0_xi2"],\
                                cov_xi_zeta["cov_xi0_zeta000"], cov_xi_zeta["cov_xi0_zeta110"]))
                C1 = np.hstack((cov_xi_xi["cov_xi2_xi0"], cov_xi_xi["cov_xi2_xi2"],\
                                cov_xi_zeta["cov_xi2_zeta000"], cov_xi_zeta["cov_xi2_zeta110"]))
                C2 = np.hstack((cov_zeta_xi["cov_zeta000_xi0"], cov_zeta_xi["cov_zeta000_xi2"],\
                                cov_zeta_zeta["cov_zeta000_zeta000"], cov_zeta_zeta["cov_zeta000_zeta110"]))
                C3 = np.hstack((cov_zeta_xi["cov_zeta110_xi0"], cov_zeta_xi["cov_zeta110_xi2"],\
                                cov_zeta_zeta["cov_zeta110_zeta000"], cov_zeta_zeta["cov_zeta110_zeta110"]))
                C = np.vstack((C0,C1,C2,C3))
                DATA_case2.update({"xi0xi2zeta000zeta110_cov": C})
                
                # Covariance matrix for case 3 #
                C0 = np.hstack((cov_xi_xi["cov_xi0_xi0"], cov_xi_xi["cov_xi0_xi2"],\
                                cov_xi_zeta["cov_xi0_zeta202"], cov_xi_zeta["cov_xi0_zeta112"]))
                C1 = np.hstack((cov_xi_xi["cov_xi2_xi0"], cov_xi_xi["cov_xi2_xi2"],\
                                cov_xi_zeta["cov_xi2_zeta202"], cov_xi_zeta["cov_xi2_zeta112"]))
                C2 = np.hstack((cov_zeta_xi["cov_zeta202_xi0"], cov_zeta_xi["cov_zeta202_xi2"],\
                                cov_zeta_zeta["cov_zeta202_zeta202"], cov_zeta_zeta["cov_zeta202_zeta112"]))
                C3 = np.hstack((cov_zeta_xi["cov_zeta112_xi0"], cov_zeta_xi["cov_zeta112_xi2"],\
                                cov_zeta_zeta["cov_zeta112_zeta202"], cov_zeta_zeta["cov_zeta112_zeta112"]))
                C = np.vstack((C0,C1,C2,C3))
                DATA_case3.update({"xi0xi2zeta202zeta112_cov": C})
                
                # Covariance matrix for case 4 #
                C0 = np.hstack((cov_xi_xi["cov_xi0_xi0"], cov_xi_xi["cov_xi0_xi2"],\
                                cov_xi_zeta["cov_xi0_zeta000"], cov_xi_zeta["cov_xi0_zeta202"]))
                C1 = np.hstack((cov_xi_xi["cov_xi2_xi0"], cov_xi_xi["cov_xi2_xi2"],\
                                cov_xi_zeta["cov_xi2_zeta000"], cov_xi_zeta["cov_xi2_zeta202"]))
                C2 = np.hstack((cov_zeta_xi["cov_zeta000_xi0"], cov_zeta_xi["cov_zeta000_xi2"],\
                                cov_zeta_zeta["cov_zeta000_zeta000"], cov_zeta_zeta["cov_zeta000_zeta202"]))
                C3 = np.hstack((cov_zeta_xi["cov_zeta202_xi0"], cov_zeta_xi["cov_zeta202_xi2"],\
                                cov_zeta_zeta["cov_zeta202_zeta000"], cov_zeta_zeta["cov_zeta202_zeta202"]))
                C = np.vstack((C0,C1,C2,C3))
                DATA_case4.update({"xi0xi2zeta000zeta202_cov": C})
                
                # Covariance matrix for case 5 #
                C0 = np.hstack((cov_xi_xi["cov_xi0_xi0"], cov_xi_xi["cov_xi0_xi2"],\
                                cov_xi_zeta["cov_xi0_zeta000"], cov_xi_zeta["cov_xi0_zeta202"], cov_xi_zeta["cov_xi0_zeta110"]))
                C1 = np.hstack((cov_xi_xi["cov_xi2_xi0"], cov_xi_xi["cov_xi2_xi2"],\
                                cov_xi_zeta["cov_xi2_zeta000"], cov_xi_zeta["cov_xi2_zeta202"], cov_xi_zeta["cov_xi2_zeta110"]))
                C2 = np.hstack((cov_zeta_xi["cov_zeta000_xi0"], cov_zeta_xi["cov_zeta000_xi2"],\
                                cov_zeta_zeta["cov_zeta000_zeta000"], cov_zeta_zeta["cov_zeta000_zeta202"], cov_zeta_zeta["cov_zeta000_zeta110"]))
                C3 = np.hstack((cov_zeta_xi["cov_zeta202_xi0"], cov_zeta_xi["cov_zeta202_xi2"],\
                                cov_zeta_zeta["cov_zeta202_zeta000"], cov_zeta_zeta["cov_zeta202_zeta202"], cov_zeta_zeta["cov_zeta202_zeta110"]))
                C4 = np.hstack((cov_zeta_xi["cov_zeta110_xi0"], cov_zeta_xi["cov_zeta110_xi2"],\
                                cov_zeta_zeta["cov_zeta110_zeta000"], cov_zeta_zeta["cov_zeta110_zeta202"], cov_zeta_zeta["cov_zeta110_zeta110"]))
                C = np.vstack((C0,C1,C2,C3,C4))
                DATA_case5.update({"xi0xi2zeta000zeta202zeta110_cov": C})
                
                # Covariance matrix for case 6 #
                C0 = np.hstack((cov_xi_xi["cov_xi0_xi0"], cov_xi_xi["cov_xi0_xi2"],\
                                cov_xi_zeta["cov_xi0_zeta000"], cov_xi_zeta["cov_xi0_zeta202"], cov_xi_zeta["cov_xi0_zeta112"]))
                C1 = np.hstack((cov_xi_xi["cov_xi2_xi0"], cov_xi_xi["cov_xi2_xi2"],\
                                cov_xi_zeta["cov_xi2_zeta000"], cov_xi_zeta["cov_xi2_zeta202"], cov_xi_zeta["cov_xi2_zeta112"]))
                C2 = np.hstack((cov_zeta_xi["cov_zeta000_xi0"], cov_zeta_xi["cov_zeta000_xi2"],\
                                cov_zeta_zeta["cov_zeta000_zeta000"], cov_zeta_zeta["cov_zeta000_zeta202"], cov_zeta_zeta["cov_zeta000_zeta112"]))
                C3 = np.hstack((cov_zeta_xi["cov_zeta202_xi0"], cov_zeta_xi["cov_zeta202_xi2"],\
                                cov_zeta_zeta["cov_zeta202_zeta000"], cov_zeta_zeta["cov_zeta202_zeta202"], cov_zeta_zeta["cov_zeta202_zeta112"]))
                C4 = np.hstack((cov_zeta_xi["cov_zeta112_xi0"], cov_zeta_xi["cov_zeta112_xi2"],\
                                cov_zeta_zeta["cov_zeta112_zeta000"], cov_zeta_zeta["cov_zeta112_zeta202"], cov_zeta_zeta["cov_zeta112_zeta112"]))
                C = np.vstack((C0,C1,C2,C3,C4))
                DATA_case6.update({"xi0xi2zeta000zeta202zeta112_cov": C})
                
                
                # Covariance matrix for case 7 #
                C0 = np.hstack((cov_xi_xi["cov_xi0_xi0"], cov_xi_xi["cov_xi0_xi2"],\
                                cov_xi_zeta["cov_xi0_zeta000"], cov_xi_zeta["cov_xi0_zeta202"],\
                                cov_xi_zeta["cov_xi0_zeta110"], cov_xi_zeta["cov_xi0_zeta112"]))
                C1 = np.hstack((cov_xi_xi["cov_xi2_xi0"], cov_xi_xi["cov_xi2_xi2"],\
                                cov_xi_zeta["cov_xi2_zeta000"], cov_xi_zeta["cov_xi2_zeta202"],\
                                cov_xi_zeta["cov_xi2_zeta110"], cov_xi_zeta["cov_xi2_zeta112"]))
                C2 = np.hstack((cov_zeta_xi["cov_zeta000_xi0"], cov_zeta_xi["cov_zeta000_xi2"],\
                                cov_zeta_zeta["cov_zeta000_zeta000"], cov_zeta_zeta["cov_zeta000_zeta202"],\
                                cov_zeta_zeta["cov_zeta000_zeta110"], cov_zeta_zeta["cov_zeta000_zeta112"]))
                C3 = np.hstack((cov_zeta_xi["cov_zeta202_xi0"], cov_zeta_xi["cov_zeta202_xi2"],\
                                cov_zeta_zeta["cov_zeta202_zeta000"], cov_zeta_zeta["cov_zeta202_zeta202"],\
                                cov_zeta_zeta["cov_zeta202_zeta110"], cov_zeta_zeta["cov_zeta202_zeta112"]))
                C4 = np.hstack((cov_zeta_xi["cov_zeta110_xi0"], cov_zeta_xi["cov_zeta110_xi2"],\
                                cov_zeta_zeta["cov_zeta110_zeta000"], cov_zeta_zeta["cov_zeta110_zeta202"],\
                                cov_zeta_zeta["cov_zeta110_zeta110"], cov_zeta_zeta["cov_zeta110_zeta112"]))
                C5 = np.hstack((cov_zeta_xi["cov_zeta112_xi0"], cov_zeta_xi["cov_zeta112_xi2"],\
                                cov_zeta_zeta["cov_zeta112_zeta000"], cov_zeta_zeta["cov_zeta112_zeta202"],\
                                cov_zeta_zeta["cov_zeta112_zeta110"], cov_zeta_zeta["cov_zeta112_zeta112"]))
                
                C = np.vstack((C0,C1,C2,C3,C4,C5))
                DATA_case7.update({"xi0xi2zeta000zeta202zeta110zeta112_cov": C})
    
                # Covariance matrix for case 8 #
                C0 = np.hstack((cov_xi_xi["cov_xi0_xi0"], cov_xi_xi["cov_xi0_xi2"],\
                                cov_xi_zeta["cov_xi0_zeta110"], cov_xi_zeta["cov_xi0_zeta112"]))
                C1 = np.hstack((cov_xi_xi["cov_xi2_xi0"], cov_xi_xi["cov_xi2_xi2"],\
                                cov_xi_zeta["cov_xi2_zeta110"], cov_xi_zeta["cov_xi2_zeta112"]))
                C2 = np.hstack((cov_zeta_xi["cov_zeta110_xi0"], cov_zeta_xi["cov_zeta110_xi2"],\
                                cov_zeta_zeta["cov_zeta110_zeta110"], cov_zeta_zeta["cov_zeta110_zeta112"]))
                C3 = np.hstack((cov_zeta_xi["cov_zeta112_xi0"], cov_zeta_xi["cov_zeta112_xi2"],\
                                cov_zeta_zeta["cov_zeta112_zeta110"], cov_zeta_zeta["cov_zeta112_zeta112"]))
                C = np.vstack((C0,C1,C2,C3))
                DATA_case8.update({"xi0xi2zeta110zeta112_cov": C})
                
                # Covariance matrix for test 1 #
                C0 = np.hstack((cov_xi_xi["cov_xi0_xi0"], cov_xi_xi["cov_xi0_xi2"], cov_xi_zeta["cov_xi0_zeta000"]))
                C1 = np.hstack((cov_xi_xi["cov_xi2_xi0"], cov_xi_xi["cov_xi2_xi2"], cov_xi_zeta["cov_xi2_zeta000"]))
                C2 = np.hstack((cov_zeta_xi["cov_zeta000_xi0"], cov_zeta_xi["cov_zeta000_xi2"], cov_zeta_zeta["cov_zeta000_zeta000"]))
                C = np.vstack((C0,C1,C2))
                DATA_test1.update({"xi0xi2zeta000_cov": C})
                
                # Covariance matrix for test 2 #
                C0 = np.hstack((cov_xi_xi["cov_xi0_xi0"], cov_xi_xi["cov_xi0_xi2"], cov_xi_zeta["cov_xi0_zeta202"]))
                C1 = np.hstack((cov_xi_xi["cov_xi2_xi0"], cov_xi_xi["cov_xi2_xi2"], cov_xi_zeta["cov_xi2_zeta202"]))
                C2 = np.hstack((cov_zeta_xi["cov_zeta202_xi0"], cov_zeta_xi["cov_zeta202_xi2"], cov_zeta_zeta["cov_zeta202_zeta202"]))
                C = np.vstack((C0,C1,C2))
                DATA_test2.update({"xi0xi2zeta202_cov": C})
                
                # Covariance matrix for test 3 #
                C0 = np.hstack((cov_xi_xi["cov_xi0_xi0"], cov_xi_xi["cov_xi0_xi2"], cov_xi_zeta["cov_xi0_zeta110"]))
                C1 = np.hstack((cov_xi_xi["cov_xi2_xi0"], cov_xi_xi["cov_xi2_xi2"], cov_xi_zeta["cov_xi2_zeta110"]))
                C2 = np.hstack((cov_zeta_xi["cov_zeta110_xi0"], cov_zeta_xi["cov_zeta110_xi2"], cov_zeta_zeta["cov_zeta110_zeta110"]))
                C = np.vstack((C0,C1,C2))
                DATA_test3.update({"xi0xi2zeta110_cov": C})
                
                # Covariance matrix for test 4 #
                C0 = np.hstack((cov_xi_xi["cov_xi0_xi0"], cov_xi_xi["cov_xi0_xi2"], cov_xi_zeta["cov_xi0_zeta112"]))
                C1 = np.hstack((cov_xi_xi["cov_xi2_xi0"], cov_xi_xi["cov_xi2_xi2"], cov_xi_zeta["cov_xi2_zeta112"]))
                C2 = np.hstack((cov_zeta_xi["cov_zeta112_xi0"], cov_zeta_xi["cov_zeta112_xi2"], cov_zeta_zeta["cov_zeta112_zeta112"]))
                C = np.vstack((C0,C1,C2))
                DATA_test4.update({"xi0xi2zeta112_cov": C})
                
                ###########################################################
                # Measurements of the 2PCF and 3PCF #
                xi0_galaxy = np.loadtxt("%s/xi0_galaxy.dat" % (INPUT), usecols=(1,), unpack=True)
                xi2_galaxy = np.loadtxt("%s/xi2_galaxy.dat" % (INPUT), usecols=(1,), unpack=True)
                zeta000_galaxy = np.loadtxt("%s/zeta000_galaxy.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                zeta202_galaxy = np.loadtxt("%s/zeta202_galaxy.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                zeta110_galaxy = np.loadtxt("%s/zeta110_galaxy.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                zeta112_galaxy = np.loadtxt("%s/zeta112_galaxy.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                
                xi0_galaxy_W1 = np.loadtxt("%s/xi0_galaxy_W1.dat" % (INPUT), usecols=(1,), unpack=True)
                xi2_galaxy_W1 = np.loadtxt("%s/xi2_galaxy_W1.dat" % (INPUT), usecols=(1,), unpack=True)
                zeta000_galaxy_W1 = np.loadtxt("%s/zeta000_galaxy_W1.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                zeta202_galaxy_W1 = np.loadtxt("%s/zeta202_galaxy_W1.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                zeta110_galaxy_W1 = np.loadtxt("%s/zeta110_galaxy_W1.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                zeta112_galaxy_W1 = np.loadtxt("%s/zeta112_galaxy_W1.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                
                xi0_galaxy_W2 = np.loadtxt("%s/xi0_galaxy_W2.dat" % (INPUT), usecols=(1,), unpack=True)
                xi2_galaxy_W2 = np.loadtxt("%s/xi2_galaxy_W2.dat" % (INPUT), usecols=(1,), unpack=True)
                zeta000_galaxy_W2 = np.loadtxt("%s/zeta000_galaxy_W2.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                zeta202_galaxy_W2 = np.loadtxt("%s/zeta202_galaxy_W2.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                zeta110_galaxy_W2 = np.loadtxt("%s/zeta110_galaxy_W2.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                zeta112_galaxy_W2 = np.loadtxt("%s/zeta112_galaxy_W2.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                
                xi0_galaxy_W3 = np.loadtxt("%s/xi0_galaxy_W3.dat" % (INPUT), usecols=(1,), unpack=True)
                xi2_galaxy_W3 = np.loadtxt("%s/xi2_galaxy_W3.dat" % (INPUT), usecols=(1,), unpack=True)
                zeta000_galaxy_W3 = np.loadtxt("%s/zeta000_galaxy_W3.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                zeta202_galaxy_W3 = np.loadtxt("%s/zeta202_galaxy_W3.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                zeta110_galaxy_W3 = np.loadtxt("%s/zeta110_galaxy_W3.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                zeta112_galaxy_W3 = np.loadtxt("%s/zeta112_galaxy_W3.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                
                xi0_mean = np.loadtxt("%s/xi0_mean.dat" % (INPUT), usecols=(1,), unpack=True)
                xi2_mean = np.loadtxt("%s/xi2_mean.dat" % (INPUT), usecols=(1,), unpack=True)
                zeta000_mean = np.loadtxt("%s/zeta000_mean.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                zeta202_mean = np.loadtxt("%s/zeta202_mean.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                zeta110_mean = np.loadtxt("%s/zeta110_mean.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                zeta112_mean = np.loadtxt("%s/zeta112_mean.dat" % (INPUT), usecols=(2,), unpack=True).reshape(rbin_max_3pcf, rbin_max_3pcf)
                
                ############
                ############
                ############
                
                xi0_galaxy = xi0_galaxy[rbin_min_2pcf:rbin_max_2pcf]
                xi2_galaxy = xi2_galaxy[rbin_min_2pcf:rbin_max_2pcf]
                # for zeta000, zeta110 and zeta112, r1>=r2. 
                zeta000_galaxy = np.array([zeta000_galaxy[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                zeta202_galaxy = np.array([zeta202_galaxy[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf)])
                zeta110_galaxy = np.array([zeta110_galaxy[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                zeta112_galaxy = np.array([zeta112_galaxy[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                
                xi0_galaxy_W1 = xi0_galaxy_W1[rbin_min_2pcf:rbin_max_2pcf]
                xi2_galaxy_W1 = xi2_galaxy_W1[rbin_min_2pcf:rbin_max_2pcf]
                # for zeta000, zeta110 and zeta112, r1>=r2. 
                zeta000_galaxy_W1 = np.array([zeta000_galaxy_W1[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                zeta202_galaxy_W1 = np.array([zeta202_galaxy_W1[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf)])
                zeta110_galaxy_W1 = np.array([zeta110_galaxy_W1[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                zeta112_galaxy_W1 = np.array([zeta112_galaxy_W1[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                
                xi0_galaxy_W2 = xi0_galaxy_W2[rbin_min_2pcf:rbin_max_2pcf]
                xi2_galaxy_W2 = xi2_galaxy_W2[rbin_min_2pcf:rbin_max_2pcf]
                # for zeta000, zeta110 and zeta112, r1>=r2. 
                zeta000_galaxy_W2 = np.array([zeta000_galaxy_W2[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                zeta202_galaxy_W2 = np.array([zeta202_galaxy_W2[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf)])
                zeta110_galaxy_W2 = np.array([zeta110_galaxy_W2[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                zeta112_galaxy_W2 = np.array([zeta112_galaxy_W2[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                
                xi0_galaxy_W3 = xi0_galaxy_W3[rbin_min_2pcf:rbin_max_2pcf]
                xi2_galaxy_W3 = xi2_galaxy_W3[rbin_min_2pcf:rbin_max_2pcf]
                # for zeta000, zeta110 and zeta112, r1>=r2. 
                zeta000_galaxy_W3 = np.array([zeta000_galaxy_W3[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                zeta202_galaxy_W3 = np.array([zeta202_galaxy_W3[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf)])
                zeta110_galaxy_W3 = np.array([zeta110_galaxy_W3[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                zeta112_galaxy_W3 = np.array([zeta112_galaxy_W3[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                
                xi0_mean = xi0_mean[rbin_min_2pcf:rbin_max_2pcf]
                xi2_mean = xi2_mean[rbin_min_2pcf:rbin_max_2pcf]
                # for zeta000, zeta110 and zeta112, r1>=r2. 
                zeta000_mean = np.array([zeta000_mean[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                zeta202_mean = np.array([zeta202_mean[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf)])
                zeta110_mean = np.array([zeta110_mean[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                zeta112_mean = np.array([zeta112_mean[i,j] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j])
                
                # Measurements for case 1
                D_galaxy = np.hstack((xi0_galaxy, xi2_galaxy))
                DATA_case1.update({"xi0xi2_galaxy": D_galaxy})
                D_galaxy_W1 = np.hstack((xi0_galaxy_W1, xi2_galaxy_W1))
                DATA_case1.update({"xi0xi2_galaxy_W1": D_galaxy_W1})
                D_galaxy_W2 = np.hstack((xi0_galaxy_W2, xi2_galaxy_W2))
                DATA_case1.update({"xi0xi2_galaxy_W2": D_galaxy_W2})
                D_galaxy_W3 = np.hstack((xi0_galaxy_W3, xi2_galaxy_W3))
                DATA_case1.update({"xi0xi2_galaxy_W3": D_galaxy_W3})
                D_mean = np.hstack((xi0_mean, xi2_mean))
                DATA_case1.update({"xi0xi2_mean": D_mean})
                
                # Measurements for case 2
                D_galaxy = np.hstack((xi0_galaxy, xi2_galaxy, zeta000_galaxy, zeta110_galaxy))
                DATA_case2.update({"xi0xi2zeta000zeta110_galaxy": D_galaxy})
                D_galaxy_W1 = np.hstack((xi0_galaxy_W1, xi2_galaxy_W1, zeta000_galaxy_W1, zeta110_galaxy_W1))
                DATA_case2.update({"xi0xi2zeta000zeta110_galaxy_W1": D_galaxy_W1})
                D_galaxy_W2 = np.hstack((xi0_galaxy_W2, xi2_galaxy_W2, zeta000_galaxy_W2, zeta110_galaxy_W2))
                DATA_case2.update({"xi0xi2zeta000zeta110_galaxy_W2": D_galaxy_W2})
                D_galaxy_W3 = np.hstack((xi0_galaxy_W3, xi2_galaxy_W3, zeta000_galaxy_W3, zeta110_galaxy_W3))
                DATA_case2.update({"xi0xi2zeta000zeta110_galaxy_W3": D_galaxy_W3})
                D_mean = np.hstack((xi0_mean, xi2_mean, zeta000_mean, zeta110_mean))
                DATA_case2.update({"xi0xi2zeta000zeta110_mean": D_mean})
                
                # Measurements for case 3
                D_galaxy = np.hstack((xi0_galaxy, xi2_galaxy, zeta202_galaxy, zeta112_galaxy))
                DATA_case3.update({"xi0xi2zeta202zeta112_galaxy": D_galaxy})
                D_galaxy_W1 = np.hstack((xi0_galaxy_W1, xi2_galaxy_W1, zeta202_galaxy_W1, zeta112_galaxy_W1))
                DATA_case3.update({"xi0xi2zeta202zeta112_galaxy_W1": D_galaxy_W1})
                D_galaxy_W2 = np.hstack((xi0_galaxy_W2, xi2_galaxy_W2, zeta202_galaxy_W2, zeta112_galaxy_W2))
                DATA_case3.update({"xi0xi2zeta202zeta112_galaxy_W2": D_galaxy_W2})
                D_galaxy_W3 = np.hstack((xi0_galaxy_W3, xi2_galaxy_W3, zeta202_galaxy_W3, zeta112_galaxy_W3))
                DATA_case3.update({"xi0xi2zeta202zeta112_galaxy_W3": D_galaxy_W3})
                D_mean = np.hstack((xi0_mean, xi2_mean, zeta202_mean, zeta112_mean))
                DATA_case3.update({"xi0xi2zeta202zeta112_mean": D_mean})
                
                # Measurements for case 4
                D_galaxy = np.hstack((xi0_galaxy, xi2_galaxy, zeta000_galaxy, zeta202_galaxy))
                DATA_case4.update({"xi0xi2zeta000zeta202_galaxy": D_galaxy})
                D_galaxy_W1 = np.hstack((xi0_galaxy_W1, xi2_galaxy_W1, zeta000_galaxy_W1, zeta202_galaxy_W1))
                DATA_case4.update({"xi0xi2zeta000zeta202_galaxy_W1": D_galaxy_W1})
                D_galaxy_W2 = np.hstack((xi0_galaxy_W2, xi2_galaxy_W2, zeta000_galaxy_W2, zeta202_galaxy_W2))
                DATA_case4.update({"xi0xi2zeta000zeta202_galaxy_W2": D_galaxy_W2})
                D_galaxy_W3 = np.hstack((xi0_galaxy_W3, xi2_galaxy_W3, zeta000_galaxy_W3, zeta202_galaxy_W3))
                DATA_case4.update({"xi0xi2zeta000zeta202_galaxy_W3": D_galaxy_W3})
                D_mean = np.hstack((xi0_mean, xi2_mean, zeta000_mean, zeta202_mean))
                DATA_case4.update({"xi0xi2zeta000zeta202_mean": D_mean})
                
                # Measurements for case 5
                D_galaxy = np.hstack((xi0_galaxy, xi2_galaxy, zeta000_galaxy, zeta202_galaxy, zeta110_galaxy))
                DATA_case5.update({"xi0xi2zeta000zeta202zeta110_galaxy": D_galaxy})
                D_galaxy_W1 = np.hstack((xi0_galaxy_W1, xi2_galaxy_W1, zeta000_galaxy_W1, zeta202_galaxy_W1, zeta110_galaxy_W1))
                DATA_case5.update({"xi0xi2zeta000zeta202zeta110_galaxy_W1": D_galaxy_W1})
                D_galaxy_W2 = np.hstack((xi0_galaxy_W2, xi2_galaxy_W2, zeta000_galaxy_W2, zeta202_galaxy_W2, zeta110_galaxy_W2))
                DATA_case5.update({"xi0xi2zeta000zeta202zeta110_galaxy_W2": D_galaxy_W2})
                D_galaxy_W3 = np.hstack((xi0_galaxy_W3, xi2_galaxy_W3, zeta000_galaxy_W3, zeta202_galaxy_W3, zeta110_galaxy_W3))
                DATA_case5.update({"xi0xi2zeta000zeta202zeta110_galaxy_W3": D_galaxy_W3})
                D_mean = np.hstack((xi0_mean, xi2_mean, zeta000_mean, zeta202_mean, zeta110_mean))
                DATA_case5.update({"xi0xi2zeta000zeta202zeta110_mean": D_mean})
                
                # Measurements for case 6
                D_galaxy = np.hstack((xi0_galaxy, xi2_galaxy, zeta000_galaxy, zeta202_galaxy, zeta112_galaxy))
                DATA_case6.update({"xi0xi2zeta000zeta202zeta112_galaxy": D_galaxy})
                D_galaxy_W1 = np.hstack((xi0_galaxy_W1, xi2_galaxy_W1, zeta000_galaxy_W1, zeta202_galaxy_W1, zeta112_galaxy_W1))
                DATA_case6.update({"xi0xi2zeta000zeta202zeta112_galaxy_W1": D_galaxy_W1})
                D_galaxy_W2 = np.hstack((xi0_galaxy_W2, xi2_galaxy_W2, zeta000_galaxy_W2, zeta202_galaxy_W2, zeta112_galaxy_W2))
                DATA_case6.update({"xi0xi2zeta000zeta202zeta112_galaxy_W2": D_galaxy_W2})
                D_galaxy_W3 = np.hstack((xi0_galaxy_W3, xi2_galaxy_W3, zeta000_galaxy_W3, zeta202_galaxy_W3, zeta112_galaxy_W3))
                DATA_case6.update({"xi0xi2zeta000zeta202zeta112_galaxy_W3": D_galaxy_W3})
                D_mean = np.hstack((xi0_mean, xi2_mean, zeta000_mean, zeta202_mean, zeta112_mean))
                DATA_case6.update({"xi0xi2zeta000zeta202zeta112_mean": D_mean})
                
                # Measurements for case 7
                D_galaxy = np.hstack((xi0_galaxy, xi2_galaxy, zeta000_galaxy, zeta202_galaxy, zeta110_galaxy, zeta112_galaxy))
                DATA_case7.update({"xi0xi2zeta000zeta202zeta110zeta112_galaxy": D_galaxy})
                D_galaxy_W1 = np.hstack((xi0_galaxy_W1, xi2_galaxy_W1, zeta000_galaxy_W1, zeta202_galaxy_W1, zeta110_galaxy_W1, zeta112_galaxy_W1))
                DATA_case7.update({"xi0xi2zeta000zeta202zeta110zeta112_galaxy_W1": D_galaxy_W1})
                D_galaxy_W2 = np.hstack((xi0_galaxy_W2, xi2_galaxy_W2, zeta000_galaxy_W2, zeta202_galaxy_W2, zeta110_galaxy_W2, zeta112_galaxy_W2))
                DATA_case7.update({"xi0xi2zeta000zeta202zeta110zeta112_galaxy_W2": D_galaxy_W2})
                D_galaxy_W3 = np.hstack((xi0_galaxy_W3, xi2_galaxy_W3, zeta000_galaxy_W3, zeta202_galaxy_W3, zeta110_galaxy_W3, zeta112_galaxy_W3))
                DATA_case7.update({"xi0xi2zeta000zeta202zeta110zeta112_galaxy_W3": D_galaxy_W3})
                D_mean = np.hstack((xi0_mean, xi2_mean, zeta000_mean, zeta202_mean, zeta110_mean, zeta112_mean))
                DATA_case7.update({"xi0xi2zeta000zeta202zeta110zeta112_mean": D_mean})
    
                # Measurements for case 8
                D_galaxy = np.hstack((xi0_galaxy, xi2_galaxy, zeta110_galaxy, zeta112_galaxy))
                DATA_case8.update({"xi0xi2zeta110zeta112_galaxy": D_galaxy})
                D_galaxy_W1 = np.hstack((xi0_galaxy_W1, xi2_galaxy_W1, zeta110_galaxy_W1, zeta112_galaxy_W1))
                DATA_case8.update({"xi0xi2zeta110zeta112_galaxy_W1": D_galaxy_W1})
                D_galaxy_W2 = np.hstack((xi0_galaxy_W2, xi2_galaxy_W2, zeta110_galaxy_W2, zeta112_galaxy_W2))
                DATA_case8.update({"xi0xi2zeta110zeta112_galaxy_W2": D_galaxy_W2})
                D_galaxy_W3 = np.hstack((xi0_galaxy_W3, xi2_galaxy_W3, zeta110_galaxy_W3, zeta112_galaxy_W3))
                DATA_case8.update({"xi0xi2zeta110zeta112_galaxy_W3": D_galaxy_W3})
                D_mean = np.hstack((xi0_mean, xi2_mean, zeta110_mean, zeta112_mean))
                DATA_case8.update({"xi0xi2zeta110zeta112_mean": D_mean})
                
                # Measurements for test 1
                D_galaxy = np.hstack((xi0_galaxy, xi2_galaxy, zeta000_galaxy))
                DATA_test1.update({"xi0xi2zeta000_galaxy": D_galaxy})
                D_galaxy_W1 = np.hstack((xi0_galaxy_W1, xi2_galaxy_W1, zeta000_galaxy_W1))
                DATA_test1.update({"xi0xi2zeta000_galaxy_W1": D_galaxy_W1})
                D_galaxy_W2 = np.hstack((xi0_galaxy_W2, xi2_galaxy_W2, zeta000_galaxy_W2))
                DATA_test1.update({"xi0xi2zeta000_galaxy_W2": D_galaxy_W2})
                D_galaxy_W3 = np.hstack((xi0_galaxy_W3, xi2_galaxy_W3, zeta000_galaxy_W3))
                DATA_test1.update({"xi0xi2zeta000_galaxy_W3": D_galaxy_W3})
                D_mean = np.hstack((xi0_mean, xi2_mean, zeta000_mean))
                DATA_test1.update({"xi0xi2zeta000_mean": D_mean})
                
                # Measurements for test 2
                D_galaxy = np.hstack((xi0_galaxy, xi2_galaxy, zeta202_galaxy))
                DATA_test2.update({"xi0xi2zeta202_galaxy": D_galaxy})
                D_galaxy_W1 = np.hstack((xi0_galaxy_W1, xi2_galaxy_W1, zeta202_galaxy_W1))
                DATA_test2.update({"xi0xi2zeta202_galaxy_W1": D_galaxy_W1})
                D_galaxy_W2 = np.hstack((xi0_galaxy_W2, xi2_galaxy_W2, zeta202_galaxy_W2))
                DATA_test2.update({"xi0xi2zeta202_galaxy_W2": D_galaxy_W2})
                D_galaxy_W3 = np.hstack((xi0_galaxy_W3, xi2_galaxy_W3, zeta202_galaxy_W3))
                DATA_test2.update({"xi0xi2zeta202_galaxy_W3": D_galaxy_W3})
                D_mean = np.hstack((xi0_mean, xi2_mean, zeta202_mean))
                DATA_test2.update({"xi0xi2zeta202_mean": D_mean})
                
                # Measurements for test 3
                D_galaxy = np.hstack((xi0_galaxy, xi2_galaxy, zeta110_galaxy))
                DATA_test3.update({"xi0xi2zeta110_galaxy": D_galaxy})
                D_galaxy_W1 = np.hstack((xi0_galaxy_W1, xi2_galaxy_W1, zeta110_galaxy_W1))
                DATA_test3.update({"xi0xi2zeta110_galaxy_W1": D_galaxy_W1})
                D_galaxy_W2 = np.hstack((xi0_galaxy_W2, xi2_galaxy_W2, zeta110_galaxy_W2))
                DATA_test3.update({"xi0xi2zeta110_galaxy_W2": D_galaxy_W2})
                D_galaxy_W3 = np.hstack((xi0_galaxy_W3, xi2_galaxy_W3, zeta110_galaxy_W3))
                DATA_test3.update({"xi0xi2zeta110_galaxy_W3": D_galaxy_W3})
                D_mean = np.hstack((xi0_mean, xi2_mean, zeta110_mean))
                DATA_test3.update({"xi0xi2zeta110_mean": D_mean})
                
                # Measurements for test 4
                D_galaxy = np.hstack((xi0_galaxy, xi2_galaxy, zeta112_galaxy))
                DATA_test4.update({"xi0xi2zeta112_galaxy": D_galaxy})
                D_galaxy_W1 = np.hstack((xi0_galaxy_W1, xi2_galaxy_W1, zeta112_galaxy_W1))
                DATA_test4.update({"xi0xi2zeta112_galaxy_W1": D_galaxy_W1})
                D_galaxy_W2 = np.hstack((xi0_galaxy_W2, xi2_galaxy_W2, zeta112_galaxy_W2))
                DATA_test4.update({"xi0xi2zeta112_galaxy_W2": D_galaxy_W2})
                D_galaxy_W3 = np.hstack((xi0_galaxy_W3, xi2_galaxy_W3, zeta112_galaxy_W3))
                DATA_test4.update({"xi0xi2zeta112_galaxy_W3": D_galaxy_W3})
                D_mean = np.hstack((xi0_mean, xi2_mean, zeta112_mean))
                DATA_test4.update({"xi0xi2zeta112_mean": D_mean})
                
                ##############################################
                # Measurements for each mock simulation data #
                ##############################################
                
                import pickle
                with open("%s/xi_R.p" % (INPUT), "rb") as pp:
                    xi_R = pickle.load(pp)
                
                with open("%s/zeta_R.p" % (INPUT), "rb") as pp:
                    zeta_R = pickle.load(pp)
                
                Rmax = xi_R["xi0"].shape[1]
                if Rmax != 2048:
                    print("Rmax = ", Rmax)
                    print("ERROR")
                    exit()
                
                xi0_RR = xi_R["xi0"]
                xi2_RR = xi_R["xi2"]
                zeta000_RR = zeta_R["zeta000"].reshape(rbin_max_3pcf, rbin_max_3pcf, Rmax)
                zeta202_RR = zeta_R["zeta202"].reshape(rbin_max_3pcf, rbin_max_3pcf, Rmax)
                zeta110_RR = zeta_R["zeta110"].reshape(rbin_max_3pcf, rbin_max_3pcf, Rmax)
                zeta112_RR = zeta_R["zeta112"].reshape(rbin_max_3pcf, rbin_max_3pcf, Rmax)
                
                xi0_RR = xi0_RR[rbin_min_2pcf:rbin_max_2pcf,:]
                xi2_RR = xi2_RR[rbin_min_2pcf:rbin_max_2pcf,:]
                # for zeta000, zeta110 and zeta112, r1>=r2. 
                zeta000_RR = np.array([zeta000_RR[i,j,:] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j ])
                zeta202_RR = np.array([zeta202_RR[i,j,:] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf)])
                zeta110_RR = np.array([zeta110_RR[i,j,:] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j ])
                zeta112_RR = np.array([zeta112_RR[i,j,:] for i in range(rbin_min_3pcf, rbin_max_3pcf) for j in range(rbin_min_3pcf, rbin_max_3pcf) if i>=j ])
                
                for RR in range(Rmax):
                
                    # Measurements for case 1
                    D_RR = np.hstack((xi0_RR[:,RR], xi2_RR[:,RR]))
                    DATA_case1.update({"xi0xi2_mock_%04d" % (RR+1): D_RR})
                
                    # Measurements for case 2
                    D_RR = np.hstack((xi0_RR[:,RR], xi2_RR[:,RR], zeta000_RR[:,RR], zeta110_RR[:,RR]))
                    DATA_case2.update({"xi0xi2zeta000zeta110_mock_%04d" % (RR+1): D_RR})
                
                    # Measurements for case 3
                    D_RR = np.hstack((xi0_RR[:,RR], xi2_RR[:,RR], zeta202_RR[:,RR], zeta112_RR[:,RR]))
                    DATA_case3.update({"xi0xi2zeta202zeta112_mock_%04d" % (RR+1): D_RR})
                
                    # Measurements for case 4
                    D_RR = np.hstack((xi0_RR[:,RR], xi2_RR[:,RR], zeta000_RR[:,RR], zeta202_RR[:,RR]))
                    DATA_case4.update({"xi0xi2zeta000zeta202_mock_%04d" % (RR+1): D_RR})
                
                    # Measurements for case 5
                    D_RR = np.hstack((xi0_RR[:,RR], xi2_RR[:,RR], zeta000_RR[:,RR], zeta202_RR[:,RR], zeta110_RR[:,RR]))
                    DATA_case5.update({"xi0xi2zeta000zeta202zeta110_mock_%04d" % (RR+1): D_RR})
                
                    # Measurements for case 6
                    D_RR = np.hstack((xi0_RR[:,RR], xi2_RR[:,RR], zeta000_RR[:,RR], zeta202_RR[:,RR], zeta112_RR[:,RR]))
                    DATA_case6.update({"xi0xi2zeta000zeta202zeta112_mock_%04d" % (RR+1): D_RR})
                
                    # Measurements for case 7
                    D_RR = np.hstack((xi0_RR[:,RR], xi2_RR[:,RR],\
                                      zeta000_RR[:,RR], zeta202_RR[:,RR],\
                                      zeta110_RR[:,RR], zeta112_RR[:,RR]))
                    DATA_case7.update({"xi0xi2zeta000zeta202zeta110zeta112_mock_%04d" % (RR+1): D_RR})
    
                    # Measurements for case 7
                    D_RR = np.hstack((xi0_RR[:,RR], xi2_RR[:,RR],\
                                      zeta110_RR[:,RR], zeta112_RR[:,RR]))
                    DATA_case8.update({"xi0xi2zeta110zeta112_mock_%04d" % (RR+1): D_RR})
    
                    # Measurements for test 1
                    D_RR = np.hstack((xi0_RR[:,RR], xi2_RR[:,RR], zeta000_RR[:,RR]))
                    DATA_test1.update({"xi0xi2zeta000_mock_%04d" % (RR+1): D_RR})
                
                    # Measurements for test 2
                    D_RR = np.hstack((xi0_RR[:,RR], xi2_RR[:,RR], zeta202_RR[:,RR]))
                    DATA_test2.update({"xi0xi2zeta202_mock_%04d" % (RR+1): D_RR})
                
                    # Measurements for test 3
                    D_RR = np.hstack((xi0_RR[:,RR], xi2_RR[:,RR], zeta110_RR[:,RR]))
                    DATA_test3.update({"xi0xi2zeta110_mock_%04d" % (RR+1): D_RR})
                
                    # Measurements for test 4
                    D_RR = np.hstack((xi0_RR[:,RR], xi2_RR[:,RR], zeta112_RR[:,RR]))
                    DATA_test4.update({"xi0xi2zeta112_mock_%04d" % (RR+1): D_RR})
                
                # Make an output directory #
                OUTPUT = "data_and_cov_for_small_scales"
                try:
                    os.mkdir(OUTPUT)
                except:
                    pass
                
                # Save data #
                
                if recon == False:
                
                    with open("%s/DATA_%s_zbin%d_case1_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case1, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_case2_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case2, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_case3_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case3, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_case4_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case4, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_case5_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case5, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_case6_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case6, pp)               
                    
                    with open("%s/DATA_%s_zbin%d_case7_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case7, pp)               
    
                    with open("%s/DATA_%s_zbin%d_case8_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case8, pp)               
                                                                    
                    with open("%s/DATA_%s_zbin%d_test1_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_test1, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_test2_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_test2, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_test3_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_test3, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_test4_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_test4, pp)
                 
                elif recon == True:
                
                    with open("%s/DATA_%s_zbin%d_recon_case1_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case1, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_recon_case2_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case2, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_recon_case3_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case3, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_recon_case4_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case4, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_recon_case5_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case5, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_recon_case6_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case6, pp)               
                    
                    with open("%s/DATA_%s_zbin%d_recon_case7_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case7, pp)               
    
                    with open("%s/DATA_%s_zbin%d_recon_case8_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_case8, pp)               
                     
                    with open("%s/DATA_%s_zbin%d_recon_test1_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_test1, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_recon_test2_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_test2, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_recon_test3_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_test3, pp)               
                                                                  
                    with open("%s/DATA_%s_zbin%d_recon_test4_r%02d.p" % (OUTPUT, NS, zbin, rbin_min_3pcf), "wb") as pp:
                        pickle.dump(DATA_test4, pp)
                    
                    
