######################
# Purpuse of this code
######################
# This code performs a test calculation of cosmological parameter estimation using two- and three-point correlation functions
######################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import montepython.io_mp as io_mp
from montepython.likelihood_class import Likelihood

import pickle

class galaxy_GENERAL_North_zbin1_case7(Likelihood):    # initialization routine

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        # Reads the measured results of the two- and three point correlation functions and the corresponding covariance matrices.
        with open("%s/%s" % (self.data_dir, self.data_file), "rb") as pp:
            data = pickle.load(pp)

        signal = data["%s%s_%s" % (self.use_xi, self.use_zeta, self.data_type)]
        nbin = len(signal)
        cov = data["%s%s_cov" % (self.use_xi, self.use_zeta)].reshape(nbin, nbin)

        ###################
        # Rescaling covariance
        rescale_factor = 1.15
        cov22 = cov[0:30,0:30]
        cov23 = cov[0:30,30:202]
        cov32 = cov[30:202,0:30]
        cov33 = rescale_factor * cov[30:202,30:202]
        cov2 = np.hstack((cov22, cov23))
        cov3 = np.hstack((cov32, cov33))
        cov = np.vstack((cov2, cov3))
        ###################

        self.cov_inv = np.matrix(cov).I
        self.signal_data = signal

        ###################
        # Hartlap factor
        nsim = self.nsim
        Hartlap = (nsim - nbin - 2.0) / (nsim - 1.0)
        self.cov_inv = Hartlap * self.cov_inv
        ###################
       
        # Reads the 2PCF and 3PCF results calculated using perturbation theory.
        with open("%s/%s" % (self.model_dir, self.model_file), "rb") as pp:
            self.model = pickle.load(pp)

    def loglkl(self, cosmo, data):

        # varying parameters 
        b1sigma8 = data.mcmc_parameters['b1sigma8']['current']
        fsigma8  = data.mcmc_parameters['fsigma8']['current']

        FGsigma8 = data.mcmc_parameters['FGsigma8']['current']
        FSsigma8 = data.mcmc_parameters['FSsigma8']['current']
        FTsigma8 = data.mcmc_parameters['FTsigma8']['current']

        GGsigma8 = data.mcmc_parameters['GGsigma8']['current']
        GSsigma8 = data.mcmc_parameters['ES']['current'] * FSsigma8
        GTsigma8 = data.mcmc_parameters['GTsigma8']['current']

        #####################################################################
        PARAM_2pcf = {}
        PARAM_2pcf.update({"b1_b1": b1sigma8**2 * fsigma8**0})
        PARAM_2pcf.update({"b1_f":  b1sigma8**1 * fsigma8**1})
        PARAM_2pcf.update({"f_f":   b1sigma8**0 * fsigma8**2})
        PARAM_NAME_2pcf = list(PARAM_2pcf.keys())
        ####################################################################
        nbin_2pcf = len(self.model["%s_%s_%s" % (self.use_xi, PARAM_NAME_2pcf[0], self.model_type)])
        xi = np.zeros(nbin_2pcf) 
        for param_name in PARAM_NAME_2pcf:
            xi = xi + PARAM_2pcf["%s" % param_name] * self.model["%s_%s_%s" % (self.use_xi, param_name, self.model_type)]

        ####################################################################
        PARAM_3pcf = {}
        PARAM_3pcf.update({"FG_b3_f0" : FGsigma8 * b1sigma8**3})
        PARAM_3pcf.update({"FS_b3_f0" : FSsigma8 * b1sigma8**3})
        PARAM_3pcf.update({"FT_b3_f0" : FTsigma8 * b1sigma8**3})
        PARAM_3pcf.update({"FG_b2_f1" : FGsigma8 * b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"FS_b2_f1" : FSsigma8 * b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"FT_b2_f1" : FTsigma8 * b1sigma8**2 * fsigma8})
        PARAM_3pcf.update({"FG_b1_f2" : FGsigma8 * b1sigma8 * fsigma8**2})
        PARAM_3pcf.update({"FS_b1_f2" : FSsigma8 * b1sigma8 * fsigma8**2})
        PARAM_3pcf.update({"FT_b1_f2" : FTsigma8 * b1sigma8 * fsigma8**2})
        
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
        PARAM_NAME_3pcf = list(PARAM_3pcf.keys())
        ####################################################################
        nbin_3pcf = len(self.model["%s_%s_%s" % (self.use_zeta, PARAM_NAME_3pcf[0], self.model_type)])
        zeta = np.zeros(nbin_3pcf) 
        for param_name in PARAM_NAME_3pcf:
            zeta = zeta + PARAM_3pcf["%s" % param_name] * self.model["%s_%s_%s" % (self.use_zeta, param_name, self.model_type)]
 
        ####################################################################

        signal_theory = np.hstack((xi, zeta))

        S = np.matrix(self.signal_data - signal_theory)

        chi2 = S * self.cov_inv * S.T

        lkl = - 0.5 * chi2

        return lkl

