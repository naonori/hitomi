#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate
import hitomipy
import pycuba
import os

class ClassPowerSpectrum():

    def __init__(self):
        self.initialize()

    def initialize(self):
        self.k_temp = np.zeros(1)
        self.P_temp = np.zeros(1)
        self.P_nw_temp = np.zeros(1)

        self.sigma8_norm = 1.0

        self.alpha_perp = 1.0
        self.alpha_parallel = 1.0
        self.sigma8 = 0.0
        self.fz = 0.0

        self.b1 = 0.0
        self.b2 = 0.0
        self.b3 = 0.0

        self.bK2 = 0.0
        self.bK3 = 0.0
        self.bDK = 0.0
        self.bO = 0.0

    def set_input_pk(self, k_in, P_in):
        self.k_temp = k_in
        self.P_temp = P_in

    def set_input_pk_nw(self, k_in, P_nw_in):
        self.k_temp = k_in
        self.P_nw_temp = P_nw_in
 
    def set_normalization(self, sigma8_norm):
        self.sigma8_norm = sigma8_norm

    def set_params(self, params):
        try:
            self.alpha_perp = params["alpha_perp"]
        except:
            pass

        try:
            self.alpha_parallel = params["alpha_parallel"]
        except:
            pass

        try:
            self.sigma8 = params["sigma8"]
        except:
            pass
        
        try:
            self.fz = params["fz"]
        except:
            pass
        
        try:
            self.b1 = params["b1"]
        except:
            pass

        try:
            self.b2 = params["b2"]
        except:
            pass

        try:
            self.b3 = params["b3"]
        except:
            pass

        try:
            self.bK2 = params["bK2"]
        except:
            pass

        try:
            self.bK3 = params["bK3"]
        except:
            pass

        try:
            self.bDK = params["bDK"]
        except:
            pass

        try:
            self.bO = params["bO"]
        except:
            pass

    def select_P(self, name):

        n_kbin = len(self.kbin)

        if name == "Tree":
            return hitomipy.integrand_P_Tree_py(
                self.xx_in, self.ndim, self.ff_out, self.ncomp,
                self.kbin, n_kbin, self.ELL,
                self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1)

        elif name == "Tree_NoWiggle":
            return hitomipy.integrand_P_Tree_NoWiggle_py(
                self.xx_in, self.ndim, self.ff_out, self.ncomp,
                self.kbin, n_kbin, self.ELL,
                self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1)

        elif name == "Tree_BAO":
            return hitomipy.integrand_P_Tree_BAO_py(
                self.xx_in, self.ndim, self.ff_out, self.ncomp,
                self.kbin, n_kbin, self.ELL,
                self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                self.sigma2_perp, self.sigma2_para)

        elif name == "Tree_BAO_Template":
            return hitomipy.integrand_P_Tree_BAO_Template_py(
                self.xx_in, self.ndim, self.ff_out, self.ncomp,
                self.kbin, n_kbin, self.ELL,
                self.alpha_perp, self.alpha_parallel,
                self.sigma2_perp, self.sigma2_para,
                self.param_name)

        elif name == "Tree_Reconstructed_Correction":
            return hitomipy.integrand_P_Tree_Reconstructed_Correction_py(
                self.xx_in, self.ndim, self.ff_out, self.ncomp,
                self.kbin, n_kbin, self.ELL,
                self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                self.sigma2_perp_recon_correction, self.sigma2_para_recon_correction,
                self.one_over_b1_fid, self.R)

        else:
            print("select_P: ERROR")

        return 0.0

    def select_ndim(self, name):

        name_dim2 = []
        name_dim2.append("Tree")
        name_dim2.append("Tree_NoWiggle")
        name_dim2.append("Tree_BAO")
        name_dim2.append("Tree_Reconstructed_Correction")
        name_dim2.append("Tree_BAO_Template")

        if name in name_dim2:
            return 2

    def check_flag_BAO(self, name, flag_BAO, sigma8_fid, fz_fid):

        flag = 0
        name_BAO = []
        name_BAO.append("Tree_BAO")
        name_BAO.append("Tree_BAO_Template")
        name_BAO.append("Tree_Reconstructed_Correction")

        name_no_BAO = []
        name_no_BAO.append("Tree")
        name_no_BAO.append("Tree_NoWiggle")

        if name in name_BAO:
            if flag_BAO and sigma8_fid >= 0.0 and fz_fid >= 0.0:
                flag = 0
                self.sigma8_fid = sigma8_fid
                self.fz_fid = fz_fid
            elif flag_BAO and sigma8_fid < 0.0 and fz_fid < 0.0:
                flag = 0
                self.sigma8_fid = self.sigma8
                self.fz_fid = self.fz
            else:
                flag = -1

        elif name in name_no_BAO:
            if flag_BAO:
                flag = -1
            else:
                flag = 0
        else:
            flag = -1

        return flag

    def check_flag_recon(self, name, flag_BAO, flag_recon, one_over_b1_fid, R):

        flag = 0
        name_recon = []
        name_recon.append("Tree_BAO")
        name_recon.append("Tree_BAO_Template")
        name_recon.append("Tree_Reconstructed_Correction")
 
        name_no_recon = []
        name_no_recon.append("Tree")
        name_no_recon.append("Tree_NoWiggle")

        if name in name_recon:
            if flag_BAO and flag_recon and one_over_b1_fid >= 0.0 and R >= 0.0:
                flag = 0
                self.one_over_b1_fid = one_over_b1_fid
                self.R = R 
            elif flag_BAO and flag_recon == False:
                flag = 0
                self.one_over_b1_fid = 1.0e-10
                self.R = 1.0e-10
            else:
                flag = - 1
         
        elif name in name_no_recon:
            if flag_recon:
                flag = -1
            else:
                flag = 0
        else:
            flag = -1

        return flag

    def check_flag_damping(self, name, flag_BAO, flag_damping, sigma2_perp, sigma2_para):

        flag = 0
        name_damping = []
        name_damping.append("Tree_BAO")
        name_damping.append("Tree_BAO_Template")
        name_damping.append("Tree_Reconstructed_Correction")

        name_no_damping = []
        name_no_damping.append("Tree")
        name_no_damping.append("Tree_NoWiggle")
        
        if name in name_damping:
            if flag_BAO and flag_damping and sigma2_perp >= 0.0 and sigma2_para >= 0.0:
                flag = 0
                self.sigma2_perp = sigma2_perp
                self.sigma2_para = sigma2_para
            elif flag_BAO and flag_damping == False:
                flag = 0
                self.sigma2_perp = 1.0e-10
                self.sigma2_para = 1.0e-10
            else:
                flag = - 1
        elif name in name_no_damping:
            if flag_damping:
                flag = -1
            else:
                flag = 0
        else:
            flag = -1

        return flag


    def Integrand_P(self, ndim, xx, ncomp, ff, userdata):

        self.xx_in = np.zeros(ndim[0])
        self.ff_out = np.zeros(ncomp[0])

        for i in range(ndim[0]):
            self.xx_in[i] = xx[i]
        
        self.ndim = ndim[0]
        self.ncomp = ncomp[0]

        self.select_P(self.name)

        for i in range(ncomp[0]):
            ff[i] = self.ff_out[i]

        return 0

    def calc_P( 
            self,  name,
            kbin=np.linspace(0.01, 0.2, 20), ELL=0,
            flag_2pcf=False,
            flag_BAO=False, sigma8_fid=-1.0, fz_fid=- 1.0,
            flag_recon = False, one_over_b1_fid = - 1.0, R = - 1.0,
            flag_damping = False, sigma2_perp = -1.0, sigma2_para = -1.0,
            param_name = None):

        ## flags ##
        output_dict_ini = {
                "kbin": kbin,
                "P": np.zeros(len(kbin)),
                "kbin_fft": kbin,
                "P_fft": np.zeros(len(kbin)),
                "ELL": ELL,
                "flag_2pcf": flag_2pcf,
                "flag_BAO": flag_BAO,
                "flag_recon": flag_recon,
                "sigma2_perp" : -1.0,
                "sigma2_para" : -1.0
                }

        check_bao = self.check_flag_BAO(name, flag_BAO, sigma8_fid, fz_fid)
        check_recon = self.check_flag_recon(name, flag_BAO, flag_recon, one_over_b1_fid, R)
        check_damping = self.check_flag_damping(name, flag_BAO, flag_damping, sigma2_perp, sigma2_para)

        if check_bao < 0:
            print("FLAG_BAO: ERROR")
            return output_dict_ini

        if check_recon < 0:
            print("FLAG_RECON: ERROR")
            return output_dict_ini

        if check_damping < 0:
            print("FLAG_DAMPING: ERROR")
            return output_dict_ini

#        ## type of power spectra, e.g,, "Tree", "NoWiggle" and etc. ##
        self.name = name
        if param_name != None and name == "Tree_BAO_Template":
            self.param_name = param_name
        elif param_name == None and name != "Tree_BAO_Template":
            pass
        else:
            print("FLAG_TEMPLATE: ERROR")
            return output_dict_ini

        ## set kbin ##
        if not flag_2pcf:
            self.kbin = kbin
        elif flag_2pcf:
            kbin0 = np.logspace(np.log(3.0e-4), np.log(0.2), 100, base=np.e)
            kbin1 = np.logspace(np.log(0.201), np.log(10.0), 50, base=np.e)
            self.kbin = np.hstack([kbin0, kbin1])

        ## set multipole indices ##
        self.ELL = ELL

        ## initialization ##
        hitomipy.initializeInputPowerSpectrum_py()

        ## read linear power spectrum ##
        hitomipy.readInputPowerSpectrum_py(
            self.k_temp, self.P_temp, len(self.k_temp))

        hitomipy.readInputNoWigglePowerSpectrum_py(
            self.k_temp, self.P_nw_temp, len(self.k_temp))

        ## normalization ##
        hitomipy.calcNormalizationUsingSigma8_py(self.sigma8_norm)
        hitomipy.calcNormalizationNoWiggle_py(1.0)

        ## sigma2_perp and sigma2_para ##
        if flag_BAO == True and flag_damping == False:

            if flag_recon == False:

                self.sigma2_perp = hitomipy.calcSigma_dd_py(self.sigma8_fid)
                self.sigma2_para = (1.0 + self.fz_fid) * (1.0 + self.fz_fid) * self.sigma2_perp

            elif flag_recon == True:

                self.sigma2_perp = pycuba.Cuhre(
                        self.Integrand_P_sigma2_perp_Reconstructed, 
                        2, 
                        ncomp=1, 
                        key=0, verbose=0 | 4)["results"][0]['integral']


                self.sigma2_para = pycuba.Cuhre(
                        self.Integrand_P_sigma2_para_Reconstructed, 
                        2, 
                        ncomp=1, 
                        key=0, verbose=0 | 4)["results"][0]['integral']

                self.sigma2_perp_recon_correction = hitomipy.calcSigma_dd_py(self.sigma8_fid)
                self.sigma2_para_recon_correction = (1.0 + self.fz_fid) * (1.0 + self.fz_fid) * self.sigma2_perp

        elif flag_BAO == False:
            self.sigma2_perp = sigma2_perp
            self.sigma2_para = sigma2_para

        ## compute power spectra ##
        NDIM = self.select_ndim(self.name)
        NCOMP = len(self.kbin)
        if NCOMP > 1024:
            print("# of NCOMP should be <= 1024, otherwise results become zero.")
            return output_dict_ini

        if NDIM == 2:
            AA = pycuba.Cuhre(self.Integrand_P, NDIM, ncomp=NCOMP, key=0, verbose=0 | 4)
        elif NDIM > 2:   
            NNEW = 2000
            NMIN = 2
            FLATNESS = 50
            MAXEVAL = 2000
            print("Suave")
            AA = pycuba.Suave(
                    self.Integrand_P, NDIM, NNEW, NMIN, FLATNESS,\
                    ncomp=NCOMP, maxeval = MAXEVAL, verbose=0 | 4)

        pk_temp = np.zeros(NCOMP)
        for i in range(NCOMP):
            pk_temp[i] = AA["results"][i]['integral']

        ## finalize parameters ##
        hitomipy.finalizeInputPowerSpectrum_py()

        if flag_2pcf:
            f_pk = interpolate.interp1d(
                self.kbin, pk_temp, fill_value="extrapolate", kind="cubic")
            kbin_out = kbin
            pk_out = f_pk(kbin)
            kbin_fft = self.kbin
            pk_fft = pk_temp
        else:
            kbin_out = self.kbin
            pk_out = pk_temp
            kbin_fft = self.kbin
            pk_fft = pk_temp

        ## output dict ##
        output_dict = {
                "kbin": kbin_out,
                "P": pk_out,
                "kbin_fft": kbin_fft,
                "P_fft": pk_fft,
                "ELL": self.ELL,
                "flag_2pcf": flag_2pcf,
                "flag_BAO": flag_BAO,
                "flag_recon": flag_recon,
                "sigma2_perp": self.sigma2_perp,
                "sigma2_para": self.sigma2_para,
        }

        return output_dict

    def calc_P_to_2PCF(self, pk_in, rbin = np.linspace(0.0, 200, 41), N_fftlog = 1000):

        if not pk_in["flag_2pcf"]:

            output_dict_ini = {
                "rbin": rbin, 
                "2pcf": np.zeros(len(rbin)), 
                "rbin_fft": rbin, 
                "2pcf_fft": np.zeros(len(rbin)), 
                "ELL": pk_in["ELL"],
                "flag_2pcf": pk_in["flag_2pcf"], 
                "flag_BAO": pk_in["flag_BAO"], 
                "flag_recon": pk_in["flag_recon"], 
                "N_fftlog": N_fftlog}
            
            return output_dict_ini

        ## set bispec. ##
        BB = pk_in["P_fft"]

        ## set kbin ##
        self.kbin = pk_in["kbin_fft"]
        
        ## set multipole indices ##
        self.ELL  = pk_in["ELL"]
        
        ## set rbin ##
        self.rbin = rbin
        
        ## input parameter for fftlog ##
        NNN = N_fftlog
        
        ## compute 2PCF ##
        kbin_for_xi = np.logspace(np.log(self.kbin[0]), np.log(self.kbin[-1]), NNN, base=np.e)
        
        f_pk = interpolate.interp1d(self.kbin, BB[:], fill_value = "extrapolate", kind="cubic")
        pk_for_xi = f_pk(kbin_for_xi)
        
        r_temp = np.zeros(NNN)
        xi_temp = np.zeros(NNN)
        hitomipy.hankel_py(self.ELL, 2, NNN, kbin_for_xi, pk_for_xi, r_temp, xi_temp)
        
        f_xi = interpolate.interp1d(r_temp, xi_temp, fill_value = "extrapolate", kind="cubic")
        
        CC = np.zeros(len(self.rbin))
        CC[:] = f_xi(self.rbin[:])
        
        CC_for_pk = np.zeros(len(r_temp))
        CC_for_pk[:] = f_xi(r_temp[:])
        
        sign = np.real(1.0j**(self.ELL))
        xi_out = sign * CC
        xi_fft = sign * CC_for_pk
        
        ## rbin_out ##
        rbin_out = self.rbin

        ## output dict ##
        output_dict = {"rbin": rbin_out, 
                       "2pcf": xi_out, 
                       "rbin_fft": r_temp, 
                       "2pcf_fft": xi_fft, 
                       "ELL": self.ELL, 
                       "flag_2pcf": pk_in["flag_2pcf"], 
                       "flag_BAO": pk_in["flag_BAO"], 
                       "flag_recon": pk_in["flag_recon"], 
                       "N_fftlog": N_fftlog}

        return output_dict

    def calc_2PCF_to_P(self, xi_in, kbin = np.linspace(0.01, 0.2, 20)):

        if not xi_in["flag_2pcf"]:
            
            output_dict_ini = {
                "kbin": kbin,
                "P": np.zeros(len(kbin)),
                "kbin_fft": kbin,
                "P_fft": np.zeros(len(kbin)),
                "ELL": xi_in["ELL"],
                "flag_2pcf": xi_in["flag_2pcf"],
                "flag_BAO": xi_in["flag_BAO"],
                "flag_recon": xi_in["flag_recon"],
                "N_fftlog": xi_in["N_fftlog"]}
    
            return output_dict_ini

        ## set powerspec. ##
        BB = xi_in["2pcf_fft"]
        
        ## set rbin ##
        self.rbin = xi_in["rbin_fft"]
        
        ## set multipole indices ##
        self.ELL  = xi_in["ELL"]
        
        ## set kbin ##
        self.kbin = kbin
        
        ## input parameter for fftlog ##
        NNN = xi_in["N_fftlog"]
        
        ## compute 3PCF ##
        rbin_for_pk = np.logspace(np.log(self.rbin[0]), np.log(self.rbin[-1]), NNN, base=np.e)
        CC = np.zeros(len(self.kbin))
        
        f_xi = interpolate.interp1d(self.rbin, BB[:], fill_value = "extrapolate", kind="cubic")
        xi_for_pk = f_xi(rbin_for_pk)
        
        k_temp = np.zeros(NNN)
        pk_temp = np.zeros(NNN)
        hitomipy.hankel_py(self.ELL, 2, NNN, rbin_for_pk, xi_for_pk, k_temp, pk_temp)
        
        f_pk = interpolate.interp1d(k_temp, pk_temp, fill_value = "extrapolate", kind="cubic")
        CC[:] = f_pk(self.kbin[:])
        
        sign = np.real((-1.0j)**(self.ELL)) * (2.0*np.pi)**3
        pk_out = sign * CC
        pk_temp = sign * pk_temp
        
        ## rbin_out ##
        kbin_out = self.kbin
        
        ## output dict ##
        output_dict = {
            "kbin": kbin_out, 
            "P": pk_out, 
            "kbin_fft": k_temp, 
            "P_fft": pk_temp, 
            "ELL": self.ELL,
            "flag_2pcf": xi_in["flag_2pcf"],
            "flag_BAO": xi_in["flag_BAO"],
            "flag_recon": xi_in["flag_recon"],
            "N_fftlog": xi_in["N_fftlog"]}
        
        return output_dict

    def Integrand_P_sigma2_perp_Reconstructed(self, ndim, xx, ncomp, ff, userdata):
        self.xx_in = np.zeros(ndim[0])
        self.ff_out = np.zeros(ncomp[0])
        for i in range(ndim[0]):
            self.xx_in[i] = xx[i]

        self.ndim = ndim[0]
        self.ncomp = ncomp[0]

        hitomipy.integrand_P_sigma2_perp_Reconstructed_py(
            self.xx_in, self.ndim, self.ff_out, self.ncomp, 
            self.kbin, len(self.kbin), 
            self.sigma8_fid, self.fz_fid, self.b1,
            self.one_over_b1_fid, self.R)

        for i in range(ncomp[0]):
            ff[i] = self.ff_out[i]
        
        return 0

    def Integrand_P_sigma2_para_Reconstructed(self, ndim, xx, ncomp, ff, userdata):
        self.xx_in = np.zeros(ndim[0])
        self.ff_out = np.zeros(ncomp[0])
        for i in range(ndim[0]):
            self.xx_in[i] = xx[i]

        self.ndim = ndim[0]
        self.ncomp = ncomp[0]

        hitomipy.integrand_P_sigma2_para_Reconstructed_py(
            self.xx_in, self.ndim, self.ff_out, self.ncomp, 
            self.kbin, len(self.kbin), 
            self.sigma8_fid, self.fz_fid, self.b1,
            self.one_over_b1_fid, self.R)

        for i in range(ncomp[0]):
            ff[i] = self.ff_out[i]
        
        return 0




