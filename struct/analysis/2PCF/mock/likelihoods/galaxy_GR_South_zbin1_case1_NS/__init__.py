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

class galaxy_GR_South_zbin1_case1_NS(Likelihood):    # initialization routine

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        # Reads the measured results of the two- and three point correlation functions and the corresponding covariance matrices.
        with open("%s/%s" % (self.data_dir, self.data_file), "rb") as pp:
            data = pickle.load(pp)

        signal = data["%s%s_%s" % (self.use_xi, self.use_zeta, self.data_type)]
        nbin = len(signal)
        cov = data["%s%s_cov" % (self.use_xi, self.use_zeta)].reshape(nbin, nbin)
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
        b1sigma8 = data.mcmc_parameters['b1sigma8_South_zbin1']['current']
        fsigma8 = data.mcmc_parameters['fsigma8_zbin1']['current']

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

        signal_theory = xi

        S = np.matrix(self.signal_data - signal_theory)

        chi2 = S * self.cov_inv * S.T

        lkl = - 0.5 * chi2

        return lkl

