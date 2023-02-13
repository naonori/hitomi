########################
# Purpose of this code #
########################
# This code calculates the 1-sigma errors of parameters through the Fisher analysis. It calculates three types of gravity theories: GR, Horndeski, and DHOST theories. It also calculates the four samples of BOSS, NGC at z=0.38, SGC at z=0.38, NGC at z=0.61, SGC at z=0.61. The minimum scale used, rmin, is calculated for six cases: rmin = 30, 40, 50, 60, 70, 80 Mpc/h. When rmin = 80 Mpc/h, it corresponds to the data used in the MCMC analysis in this work.
######################################################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pickle
import os
from classy import Class

##############################
## Compute the fisher matrix.
##############################
def calc_fisher(Gravity, NS, zbin, case, rmin, b1, sigma8, fsigma8, Omega_m):

    ####################
    # As fiducial values, input the sigma8 and fsigma8 calculated in LCDM.
    # Run 'calc_Omega_M.py' and see the results.
    # Assume b1 = 2, b2 = 0, and bs2 = 0
    ####################
    b1 = b1
    sigma8  = sigma8
    Omega_m = Omega_m
    
    b1sigma8 = b1 * sigma8
    fsigma8 = fsigma8

    FGsigma8 = (1.0 - (4.0/21.0) * Omega_m**(3.0/572.0)) * sigma8
    FSsigma8 = sigma8
    FTsigma8 = (2.0 / 7.0) * Omega_m**(3.0/572.0) * sigma8
    
    GSsigma8 = sigma8
    GTsigma8 = (4.0 / 7.0) * Omega_m**(15.0/1144.0) * sigma8

    if Gravity == "GR":
        fid_param = np.array([b1sigma8, b1sigma8, fsigma8, FGsigma8, FSsigma8, FTsigma8])
    elif Gravity == "Horndeski":
        fid_param = np.array([b1sigma8, b1sigma8, fsigma8, FGsigma8, FSsigma8, FTsigma8, GTsigma8])
    elif Gravity == "DHOST":
        fid_param = np.array([b1sigma8, b1sigma8, fsigma8, FGsigma8, FSsigma8, FTsigma8, GSsigma8, GTsigma8])

    ####################
    ## Parameters 
    ####################
    if Gravity == "GR":

        PARAM_NAME = []
        PARAM_NAME.append("b1sigma8")
        PARAM_NAME.append("fsigma8")
        PARAM_NAME.append("FGsigma8")
        PARAM_NAME.append("FSsigma8")
        PARAM_NAME.append("FTsigma8")

    elif Gravity == "Horndeski":

        PARAM_NAME = []
        PARAM_NAME.append("b1sigma8")
        PARAM_NAME.append("fsigma8")
        PARAM_NAME.append("FGsigma8")
        PARAM_NAME.append("FSsigma8")
        PARAM_NAME.append("FTsigma8")
        PARAM_NAME.append("GTsigma8")

    elif Gravity == "DHOST":

        PARAM_NAME = []
        PARAM_NAME.append("b1sigma8")
        PARAM_NAME.append("fsigma8")
        PARAM_NAME.append("FGsigma8")
        PARAM_NAME.append("FSsigma8")
        PARAM_NAME.append("FTsigma8")
        PARAM_NAME.append("GSsigma8")
        PARAM_NAME.append("GTsigma8")
    
    ####################
    ## Compute partial derivatives of 2 and 3PCFS with respect to parameters
    ####################
    signal = {}
    for param in PARAM_NAME:
        signal.update({param : compute_xi_zeta_diff(NS, zbin, case, rmin, b1sigma8, fsigma8, FGsigma8, FSsigma8, FTsigma8, GSsigma8, GTsigma8, param, Gravity)})
    
    ############
    ## Compute the inverse covariance matrix ##
    ############
    cov_inv = compute_xi_zeta_cov(NS, zbin, case, rmin)

    ############
    ## Compute the inverse Fisher matrix
    ############
    Fisher = np.zeros((len(PARAM_NAME), len(PARAM_NAME)))
    
    for i, name1 in enumerate(PARAM_NAME):
        for j, name2 in enumerate(PARAM_NAME):
            S1 = np.matrix(signal[name1])
            S2 = np.matrix(signal[name2])
    
            Fisher[i,j] = S1 * cov_inv * S2.T
    
    if case == 1:
        Fisher = np.matrix(Fisher[0:2, 0:2])
    else:
        Fisher = np.matrix(Fisher)

    Cov = Fisher.I
    Std = np.sqrt(np.diag(Cov))

    if case == 1:
        return Std
   
    else:
        #######################
        ## Variable transformatino from [fsigma8, FSsigma8, GSsigma8, GTsigma8] to [E_f, FSsigma8, E_s, E_t],
        ## where E_f = fsigma8 / FSsigma8, E_s = GSsigma8 / FSsigma8, E_t = GTsigma8 / FSsigma8
        #######################
        E_f =  fsigma8 / FSsigma8
        E_s = GSsigma8 / FSsigma8
        E_t = (7.0/4.0) * GTsigma8 / FSsigma8
    
        ## Transformation matrix ##
        Trans = np.diag(np.ones(len(PARAM_NAME)))
    
        if Gravity == "GR":
    
            Trans[1,1] = FSsigma8
            Trans[3,1] = E_f
    
        elif Gravity == "Horndeski":
    
            Trans[1,1] = FSsigma8
            Trans[3,1] = E_f
            Trans[3,5] = (4.0/7.0) * E_t
            Trans[5,5] = (4.0/7.0) * FSsigma8
    
        elif Gravity == "DHOST":
    
            Trans[1,1] = FSsigma8
            Trans[3,1] = E_f
            Trans[3,5] = E_s
            Trans[3,6] = (4.0/7.0) * E_t
            Trans[5,5] = FSsigma8
            Trans[6,6] = (4.0/7.0) * FSsigma8
    
        Trans = np.matrix(Trans)
    
        ## New Fisher matrix ##
        Fisher_new = Trans * Fisher * Trans.T
        Cov_new = np.matrix(Fisher_new).I
        Std_new = np.sqrt(np.diag(Cov_new))
    
        #######################
        #######################
        ## Variable transformatino from [E_f, E_s, E_t] to [X_f, X_s, X_t].
        ## E_f = Omega_m**(X_f), E_s = Omega_m**(X_s), and E_t = Omega_m**(X_t)
        #######################
    
        X_f = 6.0 / 11.0
        X_s = 0.0
        X_t = 15.0 / 1144.0
    
        ## Transformation matrix at zbin1 ##
        Trans2 = np.diag(np.ones(len(PARAM_NAME)))
    
        if Gravity == "GR":
            Trans2[1,1] = Omega_m**(X_f) * np.log(Omega_m)
        elif Gravity == "Horndeski":
            Trans2[1,1] = Omega_m**(X_f) * np.log(Omega_m)
            Trans2[5,5] = Omega_m**(X_t) * np.log(Omega_m)
        elif Gravity == "DHOST":
            Trans2[1,1] = Omega_m**(X_f) * np.log(Omega_m)
            Trans2[5,5] = Omega_m**(X_s) * np.log(Omega_m)
            Trans2[6,6] = Omega_m**(X_t) * np.log(Omega_m)

        Trans2 = np.matrix(Trans2)

        ## New Fisher matrix ##
        Fisher_new2 = Trans2 * Fisher_new * Trans2.T
        Cov_new2 = np.matrix(Fisher_new2).I
        Std_new2 = np.sqrt(np.diag(Cov_new2))
    
        return Std, Std_new, Std_new2
 
#
#################################################################################
# Compute partial derivatives of 2 and 3PCFs with respect to parameters   
#####
def compute_xi_zeta_diff(NS, zbin, case, rmin, b1sigma8, fsigma8, FGsigma8, FSsigma8, FTsigma8, GSsigma8, GTsigma8, diff_param, Gravity):

    model_dir  = "/mwork2/sugiymnn/WORK/model/boss/modified_gravity/model_for_small_scales"
    model_file = "MODEL_%s_zbin%d_case%d_r%02d.p" % (NS, zbin, case, rmin)
    
    use_xi  = "xi0xi2"
    data_type  = "galaxy"
    model_type  = "galaxy"

    if case == 2:
        use_zeta  = "zeta000zeta110"
    elif case == 3:
        use_zeta  = "zeta202zeta112"
    elif case == 4:
        use_zeta  = "zeta000zeta202"
    elif case == 5:
        use_zeta  = "zeta000zeta202zeta110"
    elif case == 6:
        use_zeta  = "zeta000zeta202zeta112"
    elif case == 7:
        use_zeta  = "zeta000zeta202zeta110zeta112"
    elif case == 8:
        use_zeta  = "zeta110zeta112"
    elif case == 1:
        pass
    else:
        print("ERROR")
        return exit()
    
    # Reads the 2PCF and 3PCF results calculated using perturbation theory.
    with open("%s/%s" % (model_dir, model_file), "rb") as pp:
        model = pickle.load(pp)
    
    #####################################################################
    PARAM_2pcf = {}
    PARAM_3pcf = {}
    if diff_param == "fiducial":

        PARAM_2pcf.update({"b1_b1": b1sigma8**2 * fsigma8**0})
        PARAM_2pcf.update({"b1_f":  b1sigma8**1 * fsigma8**1})
        PARAM_2pcf.update({"f_f":   b1sigma8**0 * fsigma8**2})

        ################

        PARAM_3pcf.update({"FG_b3_f0" : FGsigma8 * b1sigma8**3})
        PARAM_3pcf.update({"FS_b3_f0" : FSsigma8 * b1sigma8**3})
        PARAM_3pcf.update({"FT_b3_f0" : FTsigma8 * b1sigma8**3})
        PARAM_3pcf.update({"FG_b2_f1" : FGsigma8 * b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"FS_b2_f1" : FSsigma8 * b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"FT_b2_f1" : FTsigma8 * b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"FG_b1_f2" : FGsigma8 * b1sigma8 * fsigma8**2})
        PARAM_3pcf.update({"FS_b1_f2" : FSsigma8 * b1sigma8 * fsigma8**2})
        PARAM_3pcf.update({"FT_b1_f2" : FTsigma8 * b1sigma8 * fsigma8**2})
        
        GGsigma8 = GSsigma8 - (2.0/3.0) * GTsigma8
        PARAM_3pcf.update({"GG_b2_f1" : GGsigma8 * b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"GS_b2_f1" : GSsigma8 * b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"GT_b2_f1" : GTsigma8 * b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"GG_b1_f2" : GGsigma8 * b1sigma8**1 * fsigma8**2})
        PARAM_3pcf.update({"GS_b1_f2" : GSsigma8 * b1sigma8**1 * fsigma8**2})
        PARAM_3pcf.update({"GT_b1_f2" : GTsigma8 * b1sigma8**1 * fsigma8**2})
        PARAM_3pcf.update({"GG_b0_f3" : GGsigma8 * fsigma8**3})
        PARAM_3pcf.update({"GS_b0_f3" : GSsigma8 * fsigma8**3})
        PARAM_3pcf.update({"GT_b0_f3" : GTsigma8 * fsigma8**3})
        
        PARAM_3pcf.update({"b3_f1" : b1sigma8**3 * fsigma8})
        PARAM_3pcf.update({"b2_f2" : b1sigma8**2 * fsigma8**2})
        PARAM_3pcf.update({"b1_f3" : b1sigma8 * fsigma8**3})
        PARAM_3pcf.update({"b0_f4" : fsigma8**4})

    elif diff_param == "b1sigma8":

        PARAM_2pcf.update({"b1_b1": 2.0 * b1sigma8})
        PARAM_2pcf.update({"b1_f":  fsigma8})
        PARAM_2pcf.update({"f_f":   0.0})

        ################

        PARAM_3pcf.update({"FG_b3_f0" : FGsigma8 * 3.0 * b1sigma8**2})
        PARAM_3pcf.update({"FS_b3_f0" : FSsigma8 * 3.0 * b1sigma8**2})
        PARAM_3pcf.update({"FT_b3_f0" : FTsigma8 * 3.0 * b1sigma8**2})
        PARAM_3pcf.update({"FG_b2_f1" : FGsigma8 * 2.0 * b1sigma8**1 * fsigma8})
        PARAM_3pcf.update({"FS_b2_f1" : FSsigma8 * 2.0 * b1sigma8**1 * fsigma8})
        PARAM_3pcf.update({"FT_b2_f1" : FTsigma8 * 2.0 * b1sigma8**1 * fsigma8})
        PARAM_3pcf.update({"FG_b1_f2" : FGsigma8 * fsigma8**2})
        PARAM_3pcf.update({"FS_b1_f2" : FSsigma8 * fsigma8**2})
        PARAM_3pcf.update({"FT_b1_f2" : FTsigma8 * fsigma8**2})
        
        GGsigma8 = GSsigma8 - (2.0/3.0) * GTsigma8
        PARAM_3pcf.update({"GG_b2_f1" : GGsigma8 * 2.0 * b1sigma8**1 * fsigma8})
        PARAM_3pcf.update({"GS_b2_f1" : GSsigma8 * 2.0 * b1sigma8**1 * fsigma8})
        PARAM_3pcf.update({"GT_b2_f1" : GTsigma8 * 2.0 * b1sigma8**1 * fsigma8})
        PARAM_3pcf.update({"GG_b1_f2" : GGsigma8 * fsigma8**2})
        PARAM_3pcf.update({"GS_b1_f2" : GSsigma8 * fsigma8**2})
        PARAM_3pcf.update({"GT_b1_f2" : GTsigma8 * fsigma8**2})
        PARAM_3pcf.update({"GG_b0_f3" : 0.0})
        PARAM_3pcf.update({"GS_b0_f3" : 0.0})
        PARAM_3pcf.update({"GT_b0_f3" : 0.0})
        
        PARAM_3pcf.update({"b3_f1" : 3.0 * b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"b2_f2" : 2.0 * b1sigma8**1 * fsigma8**2})
        PARAM_3pcf.update({"b1_f3" : fsigma8**3})
        PARAM_3pcf.update({"b0_f4" : 0.0})

    elif diff_param == "fsigma8":

        PARAM_2pcf.update({"b1_b1": 0.0})
        PARAM_2pcf.update({"b1_f":  b1sigma8**1})
        PARAM_2pcf.update({"f_f":   2.0 * fsigma8**1})

        ################

        PARAM_3pcf.update({"FG_b3_f0" : 0.0})
        PARAM_3pcf.update({"FS_b3_f0" : 0.0})
        PARAM_3pcf.update({"FT_b3_f0" : 0.0})
        PARAM_3pcf.update({"FG_b2_f1" : FGsigma8 * b1sigma8**2})
        PARAM_3pcf.update({"FS_b2_f1" : FSsigma8 * b1sigma8**2})
        PARAM_3pcf.update({"FT_b2_f1" : FTsigma8 * b1sigma8**2})
        PARAM_3pcf.update({"FG_b1_f2" : FGsigma8 * b1sigma8 * 2.0 * fsigma8})
        PARAM_3pcf.update({"FS_b1_f2" : FSsigma8 * b1sigma8 * 2.0 * fsigma8})
        PARAM_3pcf.update({"FT_b1_f2" : FTsigma8 * b1sigma8 * 2.0 * fsigma8})
        
        GGsigma8 = GSsigma8 - (2.0/3.0) * GTsigma8
        PARAM_3pcf.update({"GG_b2_f1" : GGsigma8 * b1sigma8**2})
        PARAM_3pcf.update({"GS_b2_f1" : GSsigma8 * b1sigma8**2})
        PARAM_3pcf.update({"GT_b2_f1" : GTsigma8 * b1sigma8**2})
        PARAM_3pcf.update({"GG_b1_f2" : GGsigma8 * b1sigma8**1 * 2.0 * fsigma8})
        PARAM_3pcf.update({"GS_b1_f2" : GSsigma8 * b1sigma8**1 * 2.0 * fsigma8})
        PARAM_3pcf.update({"GT_b1_f2" : GTsigma8 * b1sigma8**1 * 2.0 * fsigma8})
        PARAM_3pcf.update({"GG_b0_f3" : GGsigma8 * 3.0 * fsigma8**2})
        PARAM_3pcf.update({"GS_b0_f3" : GSsigma8 * 3.0 * fsigma8**2})
        PARAM_3pcf.update({"GT_b0_f3" : GTsigma8 * 3.0 * fsigma8**2})
        
        PARAM_3pcf.update({"b3_f1" : b1sigma8**3})
        PARAM_3pcf.update({"b2_f2" : b1sigma8**2 * 2.0 * fsigma8})
        PARAM_3pcf.update({"b1_f3" : b1sigma8 * 3.0 * fsigma8**2})
        PARAM_3pcf.update({"b0_f4" : 4.0 * fsigma8**3})

    elif diff_param == "FGsigma8":

        PARAM_2pcf.update({"b1_b1": 0.0})
        PARAM_2pcf.update({"b1_f":  0.0})
        PARAM_2pcf.update({"f_f":   0.0})

        ################

        PARAM_3pcf.update({"FG_b3_f0" : b1sigma8**3})
        PARAM_3pcf.update({"FS_b3_f0" : 0.0})
        PARAM_3pcf.update({"FT_b3_f0" : 0.0})
        PARAM_3pcf.update({"FG_b2_f1" : b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"FS_b2_f1" : 0.0})
        PARAM_3pcf.update({"FT_b2_f1" : 0.0})
        PARAM_3pcf.update({"FG_b1_f2" : b1sigma8 * fsigma8**2})
        PARAM_3pcf.update({"FS_b1_f2" : 0.0})
        PARAM_3pcf.update({"FT_b1_f2" : 0.0})
        
        GGsigma8 = GSsigma8 - (2.0/3.0) * GTsigma8
        PARAM_3pcf.update({"GG_b2_f1" : 0.0})
        PARAM_3pcf.update({"GS_b2_f1" : 0.0})
        PARAM_3pcf.update({"GT_b2_f1" : 0.0})
        PARAM_3pcf.update({"GG_b1_f2" : 0.0})
        PARAM_3pcf.update({"GS_b1_f2" : 0.0})
        PARAM_3pcf.update({"GT_b1_f2" : 0.0})
        PARAM_3pcf.update({"GG_b0_f3" : 0.0})
        PARAM_3pcf.update({"GS_b0_f3" : 0.0})
        PARAM_3pcf.update({"GT_b0_f3" : 0.0})
        
        PARAM_3pcf.update({"b3_f1" : 0.0})
        PARAM_3pcf.update({"b2_f2" : 0.0})
        PARAM_3pcf.update({"b1_f3" : 0.0})
        PARAM_3pcf.update({"b0_f4" : 0.0})

    elif diff_param == "FSsigma8":

        PARAM_2pcf.update({"b1_b1": 0.0})
        PARAM_2pcf.update({"b1_f":  0.0})
        PARAM_2pcf.update({"f_f":   0.0})

        ################

        PARAM_3pcf.update({"FG_b3_f0" : 0.0})
        PARAM_3pcf.update({"FS_b3_f0" : b1sigma8**3})
        PARAM_3pcf.update({"FT_b3_f0" : 0.0})
        PARAM_3pcf.update({"FG_b2_f1" : 0.0})
        PARAM_3pcf.update({"FS_b2_f1" : b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"FT_b2_f1" : 0.0})
        PARAM_3pcf.update({"FG_b1_f2" : 0.0})
        PARAM_3pcf.update({"FS_b1_f2" : b1sigma8 * fsigma8**2})
        PARAM_3pcf.update({"FT_b1_f2" : 0.0})
        
        if Gravity == "DHOST":
            GGsigma8 = GSsigma8 - (2.0/3.0) * GTsigma8
            PARAM_3pcf.update({"GG_b2_f1" : 0.0})
            PARAM_3pcf.update({"GS_b2_f1" : 0.0})
            PARAM_3pcf.update({"GT_b2_f1" : 0.0})
            PARAM_3pcf.update({"GG_b1_f2" : 0.0})
            PARAM_3pcf.update({"GS_b1_f2" : 0.0})
            PARAM_3pcf.update({"GT_b1_f2" : 0.0})
            PARAM_3pcf.update({"GG_b0_f3" : 0.0})
            PARAM_3pcf.update({"GS_b0_f3" : 0.0})
            PARAM_3pcf.update({"GT_b0_f3" : 0.0})
        elif Gravity == "Horndeski":
            GSsigma8 = 1.0
            GGsigma8 = 1.0
            GTsigma8 = 0.0
            PARAM_3pcf.update({"GG_b2_f1" : GGsigma8 * b1sigma8**2 * fsigma8})
            PARAM_3pcf.update({"GS_b2_f1" : GSsigma8 * b1sigma8**2 * fsigma8})
            PARAM_3pcf.update({"GT_b2_f1" : GTsigma8 * b1sigma8**2 * fsigma8})
            PARAM_3pcf.update({"GG_b1_f2" : GGsigma8 * b1sigma8**1 * fsigma8**2})
            PARAM_3pcf.update({"GS_b1_f2" : GSsigma8 * b1sigma8**1 * fsigma8**2})
            PARAM_3pcf.update({"GT_b1_f2" : GTsigma8 * b1sigma8**1 * fsigma8**2})
            PARAM_3pcf.update({"GG_b0_f3" : GGsigma8 * fsigma8**3})
            PARAM_3pcf.update({"GS_b0_f3" : GSsigma8 * fsigma8**3})
            PARAM_3pcf.update({"GT_b0_f3" : GTsigma8 * fsigma8**3})
        elif Gravity == "GR":
            GSsigma8 = 1.0
            GGsigma8 = (1.0 - (8.0/21.0) * Omega_m**(15.0/1144.0))
            GTsigma8 = (4.0/7.0) * Omega_m**(15.0/1144.0)
            PARAM_3pcf.update({"GG_b2_f1" : GGsigma8 * b1sigma8**2 * fsigma8})
            PARAM_3pcf.update({"GS_b2_f1" : GSsigma8 * b1sigma8**2 * fsigma8})
            PARAM_3pcf.update({"GT_b2_f1" : GTsigma8 * b1sigma8**2 * fsigma8})
            PARAM_3pcf.update({"GG_b1_f2" : GGsigma8 * b1sigma8**1 * fsigma8**2})
            PARAM_3pcf.update({"GS_b1_f2" : GSsigma8 * b1sigma8**1 * fsigma8**2})
            PARAM_3pcf.update({"GT_b1_f2" : GTsigma8 * b1sigma8**1 * fsigma8**2})
            PARAM_3pcf.update({"GG_b0_f3" : GGsigma8 * fsigma8**3})
            PARAM_3pcf.update({"GS_b0_f3" : GSsigma8 * fsigma8**3})
            PARAM_3pcf.update({"GT_b0_f3" : GTsigma8 * fsigma8**3})
        
        PARAM_3pcf.update({"b3_f1" : 0.0})
        PARAM_3pcf.update({"b2_f2" : 0.0})
        PARAM_3pcf.update({"b1_f3" : 0.0})
        PARAM_3pcf.update({"b0_f4" : 0.0})

    elif diff_param == "FTsigma8":

        PARAM_2pcf.update({"b1_b1": 0.0})
        PARAM_2pcf.update({"b1_f":  0.0})
        PARAM_2pcf.update({"f_f":   0.0})

        ################

        PARAM_3pcf.update({"FG_b3_f0" : 0.0})
        PARAM_3pcf.update({"FS_b3_f0" : 0.0})
        PARAM_3pcf.update({"FT_b3_f0" : b1sigma8**3})
        PARAM_3pcf.update({"FG_b2_f1" : 0.0})
        PARAM_3pcf.update({"FS_b2_f1" : 0.0})
        PARAM_3pcf.update({"FT_b2_f1" : b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"FG_b1_f2" : 0.0})
        PARAM_3pcf.update({"FS_b1_f2" : 0.0})
        PARAM_3pcf.update({"FT_b1_f2" : b1sigma8 * fsigma8**2})
        
        GGsigma8 = GSsigma8 - (2.0/3.0) * GTsigma8
        PARAM_3pcf.update({"GG_b2_f1" : 0.0})
        PARAM_3pcf.update({"GS_b2_f1" : 0.0})
        PARAM_3pcf.update({"GT_b2_f1" : 0.0})
        PARAM_3pcf.update({"GG_b1_f2" : 0.0})
        PARAM_3pcf.update({"GS_b1_f2" : 0.0})
        PARAM_3pcf.update({"GT_b1_f2" : 0.0})
        PARAM_3pcf.update({"GG_b0_f3" : 0.0})
        PARAM_3pcf.update({"GS_b0_f3" : 0.0})
        PARAM_3pcf.update({"GT_b0_f3" : 0.0})
        
        PARAM_3pcf.update({"b3_f1" : 0.0})
        PARAM_3pcf.update({"b2_f2" : 0.0})
        PARAM_3pcf.update({"b1_f3" : 0.0})
        PARAM_3pcf.update({"b0_f4" : 0.0})

    elif diff_param == "GSsigma8":

        PARAM_2pcf.update({"b1_b1": 0.0})
        PARAM_2pcf.update({"b1_f":  0.0})
        PARAM_2pcf.update({"f_f":   0.0})

        ################

        PARAM_3pcf.update({"FG_b3_f0" : 0.0})
        PARAM_3pcf.update({"FS_b3_f0" : 0.0})
        PARAM_3pcf.update({"FT_b3_f0" : 0.0})
        PARAM_3pcf.update({"FG_b2_f1" : 0.0})
        PARAM_3pcf.update({"FS_b2_f1" : 0.0})
        PARAM_3pcf.update({"FT_b2_f1" : 0.0})
        PARAM_3pcf.update({"FG_b1_f2" : 0.0})
        PARAM_3pcf.update({"FS_b1_f2" : 0.0})
        PARAM_3pcf.update({"FT_b1_f2" : 0.0})
        
        GGsigma8 = 1.0 
        PARAM_3pcf.update({"GG_b2_f1" : GGsigma8 * b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"GS_b2_f1" : b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"GT_b2_f1" : 0.0})
        PARAM_3pcf.update({"GG_b1_f2" : GGsigma8 * b1sigma8**1 * fsigma8**2})
        PARAM_3pcf.update({"GS_b1_f2" : b1sigma8**1 * fsigma8**2})
        PARAM_3pcf.update({"GT_b1_f2" : 0.0})
        PARAM_3pcf.update({"GG_b0_f3" : GGsigma8 * fsigma8**3})
        PARAM_3pcf.update({"GS_b0_f3" : fsigma8**3})
        PARAM_3pcf.update({"GT_b0_f3" : 0.0})
        
        PARAM_3pcf.update({"b3_f1" : 0.0})
        PARAM_3pcf.update({"b2_f2" : 0.0})
        PARAM_3pcf.update({"b1_f3" : 0.0})
        PARAM_3pcf.update({"b0_f4" : 0.0})

    elif diff_param == "GTsigma8":

        PARAM_2pcf.update({"b1_b1": 0.0})
        PARAM_2pcf.update({"b1_f":  0.0})
        PARAM_2pcf.update({"f_f":   0.0})

        ################

        PARAM_3pcf.update({"FG_b3_f0" : 0.0})
        PARAM_3pcf.update({"FS_b3_f0" : 0.0})
        PARAM_3pcf.update({"FT_b3_f0" : 0.0})
        PARAM_3pcf.update({"FG_b2_f1" : 0.0})
        PARAM_3pcf.update({"FS_b2_f1" : 0.0})
        PARAM_3pcf.update({"FT_b2_f1" : 0.0})
        PARAM_3pcf.update({"FG_b1_f2" : 0.0})
        PARAM_3pcf.update({"FS_b1_f2" : 0.0})
        PARAM_3pcf.update({"FT_b1_f2" : 0.0})
        
        GGsigma8 = - (2.0/3.0)
        PARAM_3pcf.update({"GG_b2_f1" : GGsigma8 * b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"GS_b2_f1" : 0.0})
        PARAM_3pcf.update({"GT_b2_f1" : b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"GG_b1_f2" : GGsigma8 * b1sigma8**1 * fsigma8**2})
        PARAM_3pcf.update({"GS_b1_f2" : 0.0})
        PARAM_3pcf.update({"GT_b1_f2" : b1sigma8**1 * fsigma8**2})
        PARAM_3pcf.update({"GG_b0_f3" : GGsigma8 * fsigma8**3})
        PARAM_3pcf.update({"GS_b0_f3" : 0.0})
        PARAM_3pcf.update({"GT_b0_f3" : fsigma8**3})
        
        PARAM_3pcf.update({"b3_f1" : 0.0})
        PARAM_3pcf.update({"b2_f2" : 0.0})
        PARAM_3pcf.update({"b1_f3" : 0.0})
        PARAM_3pcf.update({"b0_f4" : 0.0})

    elif diff_param == "else":

        PARAM_2pcf.update({"b1_b1": 0.0})
        PARAM_2pcf.update({"b1_f":  0.0})
        PARAM_2pcf.update({"f_f":   0.0})

        ################

        PARAM_3pcf.update({"FG_b3_f0" : 0.0})
        PARAM_3pcf.update({"FS_b3_f0" : 0.0})
        PARAM_3pcf.update({"FT_b3_f0" : 0.0})
        PARAM_3pcf.update({"FG_b2_f1" : 0.0})
        PARAM_3pcf.update({"FS_b2_f1" : 0.0})
        PARAM_3pcf.update({"FT_b2_f1" : 0.0})
        PARAM_3pcf.update({"FG_b1_f2" : 0.0})
        PARAM_3pcf.update({"FS_b1_f2" : 0.0})
        PARAM_3pcf.update({"FT_b1_f2" : 0.0})
        
        PARAM_3pcf.update({"GG_b2_f1" : 0.0})
        PARAM_3pcf.update({"GS_b2_f1" : 0.0})
        PARAM_3pcf.update({"GT_b2_f1" : 0.0})
        PARAM_3pcf.update({"GG_b1_f2" : 0.0})
        PARAM_3pcf.update({"GS_b1_f2" : 0.0})
        PARAM_3pcf.update({"GT_b1_f2" : 0.0})
        PARAM_3pcf.update({"GG_b0_f3" : 0.0})
        PARAM_3pcf.update({"GS_b0_f3" : 0.0})
        PARAM_3pcf.update({"GT_b0_f3" : 0.0})
        
        PARAM_3pcf.update({"b3_f1" : 0.0})
        PARAM_3pcf.update({"b2_f2" : 0.0})
        PARAM_3pcf.update({"b1_f3" : 0.0})
        PARAM_3pcf.update({"b0_f4" : 0.0})

    PARAM_NAME_2pcf = list(PARAM_2pcf.keys())
    ####################################################################
    nbin_2pcf = len(model["%s_%s_%s" % (use_xi, PARAM_NAME_2pcf[0], model_type)])
    xi = np.zeros(nbin_2pcf) 
    for param_name in PARAM_NAME_2pcf:
        xi = xi + PARAM_2pcf["%s" % param_name] * model["%s_%s_%s" % (use_xi, param_name, model_type)]
    
    if case == 1:
        signal = xi
        return signal
    else:
        ####################################################################
        PARAM_NAME_3pcf = list(PARAM_3pcf.keys())
        ####################################################################
        nbin_3pcf = len(model["%s_%s_%s" % (use_zeta, PARAM_NAME_3pcf[0], model_type)])
        zeta = np.zeros(nbin_3pcf) 
        for param_name in PARAM_NAME_3pcf:
            zeta = zeta + PARAM_3pcf["%s" % param_name] * model["%s_%s_%s" % (use_zeta, param_name, model_type)]
        
        signal = np.hstack((xi, zeta))
        
        return signal
    
#################################################################################
# Compute covariance matrices of 2 and 3PCFs.
#####
def compute_xi_zeta_cov(NS, zbin, case, rmin):

    data_dir  = "/mwork2/sugiymnn/WORK/covariance/boss/modified_gravity/data_and_cov_for_small_scales"
    data_file = "DATA_%s_zbin%d_case%d_r%02d.p" % (NS, zbin, case, rmin)

    use_xi = "xi0xi2"
    data_type = "galaxy"
    model_type = "galaxy"

    if case == 2:
        use_zeta  = "zeta000zeta110"
    elif case == 3:
        use_zeta  = "zeta202zeta112"
    elif case == 4:
        use_zeta  = "zeta000zeta202"
    elif case == 5:
        use_zeta  = "zeta000zeta202zeta110"
    elif case == 6:
        use_zeta  = "zeta000zeta202zeta112"
    elif case == 7:
        use_zeta  = "zeta000zeta202zeta110zeta112"
    elif case == 8:
        use_zeta  = "zeta110zeta112"
    elif case == 1:
        pass
    else:
        exit()
 
    with open("%s/%s" % (data_dir, data_file), "rb") as pp:
        data = pickle.load(pp)

    if case == 1:
        cov = data["%s_cov" % (use_xi)]
    else:
        cov = data["%s%s_cov" % (use_xi, use_zeta)]

    cov_inv = np.matrix(cov).I

    ###################
    # Hartlap factor
    nsim = 2048.0
    nbin = len(cov)
    Hartlap = (nsim - nbin - 2.0) / (nsim - 1.0)
    cov_inv = Hartlap * cov_inv
    ###################

    return cov_inv

def calc_Omega_M(NS, zbin, data):

    ############
    # Redshift #
    ############
    # Calculates the theoretical models for the first and third BOSS galaxy samples when divided into three redshift ranges:
    # 0.2 < z < 0.5, 0.4 < z < 0.6, and 0.5 < z < 0.75.
    # (see $WORK/data/boss/fits2dat_galaxy.py)
    ##############################################
    if zbin == 1:
        redshift = 0.38 
    elif zbin == 3:
        redshift = 0.61
    
    if data == "galaxy":
        ################################################################################
        # BOSS fiducial parameters #
        # The values of the fiducial parameters should be exactly the same as those used to calculate the distance to the BOSS galaxies.
        # (see $WORK/data/boss/fits2dat_galaxy.py)
        ################################################################################
        params_cosmo = {
            'h': 0.676,
            'omega_b': 0.022,
            'omega_cdm': 0.11966256,
            'n_s': 0.96,
            'ln10^{10}A_s': 3.094,
            'tau_reio': 0.0826026,
            
            'N_ur': 2.0328,
            'N_ncdm': 1,
            'm_ncdm': 0.06,
            'T_ncdm': 0.71611,
        
            'output': 'tCl, mPk',
            'P_k_max_h/Mpc': 50.0,
            'z_max_pk': 3.0,
        }
    
    elif data == "mock":
    
        ################################################################################
        # Patchy mock fiducial parameters #
        # The values of the fiducial parameters should be exactly the same as those used to calculate the distance to the Patchy mock galaxies.
        # (see $WORK/data/boss/fits2dat_galaxy.py)
        ################################################################################
    
        h        = 0.6777
        Omega_b  = 0.0480
        Omega_m  = 0.3071
        params_cosmo = {
            'h': h,
            'omega_b': Omega_b * h**2,
            'omega_cdm': (Omega_m - Omega_b) * h**2,
            'n_s': 0.9645,
            'ln10^{10}A_s': 3.094,
            'tau_reio': 0.0826026,
        
            'output': 'tCl, mPk',
            'P_k_max_h/Mpc': 50.0,
            'z_max_pk': 3.0,
        
        }
    
    else:
        print("ERROR")
        exit()
    
    ###############
    # Call Class() 
    ###############
    cosmo = Class()
    cosmo.set(params_cosmo)
    cosmo.compute()
    
    Omega_m = cosmo.Omega_m() * ( cosmo.Hubble(0.0) / cosmo.Hubble(redshift) )**2 * (1.0 + redshift)**3
    f = cosmo.scale_independent_growth_factor_f(redshift)
    D = cosmo.scale_independent_growth_factor(redshift)
    sigma8 = cosmo.sigma(8.0/cosmo.h(), redshift)
 
    return Omega_m, f * sigma8, sigma8

if __name__ == "__main__":

    # output directory #
    output_dir = "results_fisher"
    try:
        os.mkdir(output_dir)
    except:
        pass

    data = "galaxy"
    # The fiducial value of the linear bias is set to 2.0 #
    # Here, we ignore any non-linear bias parameters: b2 = bs2 = 0 #
    b1 = 2.0

    # The minimum scales for Fisher analysis: rmin #
    rbin = np.array([30.0, 40.0, 50.0, 60.0, 70.0, 80.0])

    for NS in ["North", "South"]:
        for zbin in [1,3]:

            print("%s zbin = %d" % (NS,zbin))
            Omega_m, fsigma8, sigma8 = calc_Omega_M(NS, zbin, data)

            print("Omega_m = %.2f fsigma8 = %.2f sigma8 = %.2f" % (Omega_m, fsigma8, sigma8))

            # Fiducial values of parameters #
            b1sigma8 = b1 * sigma8
            fsigma8 = fsigma8
        
            FGsigma8 = (1.0 - (4.0/21.0) * Omega_m**(3.0/572.0)) * sigma8
            FSsigma8 = sigma8
            FTsigma8 = (2.0 / 7.0) * Omega_m**(3.0/572.0) * sigma8
        
            GSsigma8 = sigma8
            GTsigma8 = (4.0 / 7.0) * Omega_m**(15.0/1144.0) * sigma8

            E_f =  fsigma8 / FSsigma8
            E_s = GSsigma8 / FSsigma8
            E_t = (7.0 / 4.0) * GTsigma8 / FSsigma8
    
            X_f = 6.0 / 11.0
            X_s = 0.0
            X_t = 15.0 / 1144.0

            X = np.array([b1sigma8, fsigma8, FGsigma8, FSsigma8, FTsigma8, GSsigma8, GTsigma8, E_f, E_s, E_t, X_f, X_s, X_t, Omega_m]).T
            np.savetxt("results_fisher/fiducial_values_%s_zbin%d.dat" % (NS, zbin), X, fmt="%.7e", header="b1sigma8 \t fsigma8 \t FGsigma8 \t FSsigma8 \t FTsigma8 \t GSsigma8 \t GTsigma8 \t E_f \t E_s \t E_t \t X_f \t X_s \t X_t \t Omega_m")

            ## Calculate the expected 1 sigma errors in the parameters through the Fisher matrix ##
            for Gravity in ["GR", "Horndeski", "DHOST"]:

                if Gravity == "GR":
                    case = 1
                    Std = np.zeros((6,2))
                    for rmin in [0,1,2,3,4,5]:
                        Std_temp = calc_fisher(Gravity, NS, zbin, case, rmin, b1, sigma8, fsigma8, Omega_m)
                        Std[rmin,:] = Std_temp[:]
                    X = np.array([rbin, Std[:,0], Std[:,1]]).T
                    np.savetxt("results_fisher/fisher_%s_%s_zbin%d_case%d.dat" % (Gravity, NS, zbin, case), X,\
                            fmt="%.7e", header="rmin [Mpc/h] \t b1sigma8 \t fsigma8")

                for case in [2,3,4,5,6,7,8]:
                    if Gravity == "GR":
                        Std = np.zeros((6,5))
                        Std_new = np.zeros((6,5))
                        Std_new2 = np.zeros((6,5))
                    elif Gravity == "Horndeski":
                        Std = np.zeros((6,6))
                        Std_new = np.zeros((6,6))
                        Std_new2 = np.zeros((6,6))
                    elif Gravity == "DHOST":
                        Std = np.zeros((6,7))
                        Std_new = np.zeros((6,7))
                        Std_new2 = np.zeros((6,7))
                    for rmin in [0,1,2,3,4,5]:

                        Std_temp, Std_new_temp, Std_new2_temp =\
                                calc_fisher(Gravity, NS, zbin, case, rmin, b1, sigma8, fsigma8, Omega_m)
                        
                        Std[rmin,:] = Std_temp[:]
                        Std_new[rmin,:] = Std_new_temp[:]
                        Std_new2[rmin,:] = Std_new2_temp[:]

                    if Gravity == "GR":
                        X = np.array([rbin, Std[:,0], Std[:,1], Std[:,2], Std[:,3], Std[:,4], Std_new[:,1], Std_new2[:,1]]).T
                        np.savetxt("results_fisher/fisher_%s_%s_zbin%d_case%d.dat" % (Gravity, NS, zbin, case), X,\
                                fmt="%.7e", header="rmin [Mpc/h] \t b1sigma8 \t fsigma8 \t FGsigma8 \t FSsigma8 \t FTsigma8 \t E_f \t X_f")
                    elif Gravity == "Horndeski":
                        X = np.array([rbin, Std[:,0], Std[:,1], Std[:,2], Std[:,3], Std[:,4], Std[:,5], Std_new[:,1], Std_new[:,5], Std_new2[:,1], Std_new2[:,5]]).T
                        np.savetxt("results_fisher/fisher_%s_%s_zbin%d_case%d.dat" % (Gravity, NS, zbin, case), X,\
                                fmt="%.7e", header="rmin [Mpc/h] \t b1sigma8 \t fsigma8 \t FGsigma8 \t FSsigma8 \t FTsigma8 \t GTsigma8 \t E_f \t E_t \t X_f \t X_t")
                    elif Gravity == "DHOST":
                        X = np.array([rbin, Std[:,0], Std[:,1], Std[:,2], Std[:,3], Std[:,4], Std[:,5], Std[:,6], Std_new[:,1], Std_new[:,5], Std_new[:,6], Std_new2[:,1], Std_new2[:,5], Std_new2[:,6]]).T
                        np.savetxt("results_fisher/fisher_%s_%s_zbin%d_case%d.dat" % (Gravity, NS, zbin, case), X,\
                                fmt="%.7e", header="rmin [Mpc/h] \t b1sigma8 \t fsigma8 \t FGsigma8 \t FSsigma8 \t FTsigma8 \t GSsigma8 \t GTsigma8 \t E_f \t E_s \t E_t \t X_f \t X_s \t X_t")


