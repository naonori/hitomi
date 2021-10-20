#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import hitomipy

class InputPowerSpectrum():

    def __init__(self, redshift, cosmo):
        self.cosmo = cosmo
        if redshift < 1.0e-10:
            self.redshift = 1.0e-10
        else:
            self.redshift = redshift
        ## matter power spectrum ##
        self.kmax = 50.0
        self.kmin = 2.0e-5
        self.num_kbin = 500
        self.k_pivot = 0.05
        ln_kmin = np.log(self.kmin)
        ln_kmax = np.log(self.kmax)
        self.k = np.logspace(ln_kmin, ln_kmax, self.num_kbin, base=np.e)
        self.pk = np.zeros(len(self.k))
        self.pk_nw = np.zeros(len(self.k))
        self.pk_pri = np.zeros(len(self.k))
        self.Dz = 1.0
        self.fz = 0.0
        self.sigma8_norm = 1.0
        self.sigma8_0 = 0.0
        self.H = 0.0
        self.Da = 0.0
        self._c_ = 2.99792458e5  # [ km/s ] ##

    def calcNoWiggleMatterPowerSpectrum(self):

        h = self.cosmo.h()
        omega0 = self.cosmo.Omega_m()
        omegab = self.cosmo.Omega_b()
        Tcmb = self.cosmo.T_cmb()
        n_s = self.cosmo.n_s()

        for i in range(self.num_kbin):
            self.pk_nw[i] = hitomipy.calcNoWiggleMatterPowerSpectrum_py(self.k[i], h, omega0, omegab, Tcmb, n_s)

    def calcMatterPowerSpectrum(self):

        h = self.cosmo.h()

        ## compute linear matter power spectrum ##
        for i in range(self.num_kbin):
            self.pk[i] = self.cosmo.pk_lin(self.k[i] * h, self.redshift) * h**3

        self.Da = self.cosmo.angular_distance(self.redshift) * h
        self.H = self.cosmo.Hubble(self.redshift) * self._c_
        self.Dz = self.cosmo.scale_independent_growth_factor(self.redshift)
        self.fz = self.cosmo.scale_independent_growth_factor_f(self.redshift)
        self.sigma8_norm = self.cosmo.sigma(8.0 / h, self.redshift)


    def calcPrimordialPowerSpectrum(self, ln10A_s10, n_s):

        h = self.cosmo.h()

        def f_pk_pri(k):
            A_s = 1.0e-10 * np.exp(ln10A_s10)
            P_pri = (3.0 / 5.0)**2 * (2.0 * np.pi**2 / k**3) * \
                    A_s * (k / (self.k_pivot * h))**(n_s - 1.0)
            return P_pri

        for i in range(self.num_kbin):
            self.pk_pri[i] = f_pk_pri(self.k[i] * h) * h**3

    def getMatterPowerSpectrum(self):
        return self.k, self.pk

    def getNoWiggleMatterPowerSpectrum(self):
        return self.k, self.pk_nw

    def getTransferFunctionM(self):
        return self.k, np.sqrt(self.pk / self.pk_pri)

    def getGrowthRate(self):
        return self.fz

    def getGrowthFactor(self):
        return self.Dz

    def getSigma8z(self, sigma8_0=-1.0):
        if sigma8_0 < 0.0:
            return self.sigma8_norm
        else:
            return sigma8_0 * self.Dz

    def getSigma8ForNormalization(self):
        return self.sigma8_norm

    def getAngularDiameterDistance(self):
        return self.Da

    def getHubbleParameter(self):
        return self.H

    def getAlphaPerp(self, Da_fid):
        return self.Da / Da_fid

    def getAlphaParallel(self, H_fid):
        return H_fid / self.H


