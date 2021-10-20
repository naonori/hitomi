#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate
import hitomipy
import pycuba
import os

class ClassBiSpectrum():

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


    def select_B(self, name):

        n_kbin = len(self.kbin)

        if name == "Tree":
            return hitomipy.integrand_B_Tree_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ell1, self.ell2, self.ELL, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, 
                    self.sigma8, self.fz, self.b1,
                    self.b2, self.bK2)

        elif name == "Tree_NoWiggle":
            return hitomipy.integrand_B_Tree_NoWiggle_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ell1, self.ell2, self.ELL, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, 
                    self.sigma8, self.fz, self.b1,
                    self.b2, self.bK2)

        elif name == "Tree_BAO":
            return hitomipy.integrand_B_Tree_BAO_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ell1, self.ell2, self.ELL, self.kmag1,
                    self.alpha_perp, self.alpha_parallel,
                    self.sigma8, self.fz, self.b1,
                    self.b2, self.bK2,
                    self.sigma2_perp, self.sigma2_para)

        elif name == "Tree_BAO_Reconstructed":
            return hitomipy.integrand_B_Tree_BAO_Reconstructed_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ell1, self.ell2, self.ELL, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, 
                    self.sigma8, self.fz, self.b1,
                    self.b2, self.bK2,
                    self.one_over_b1_fid, self.R,
                    self.sigma2_perp, self.sigma2_para)

        elif name == "Tree_BAO_Template":
            return hitomipy.integrand_B_Tree_BAO_Template_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ell1, self.ell2, self.ELL, self.kmag1,
                    self.alpha_perp, self.alpha_parallel,
                    self.sigma2_perp, self.sigma2_para,
                    self.param_name)

        elif name == "Tree_NonGaussian_Local":
            return hitomipy.integrand_B_Tree_NonGaussian_Local_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ell1, self.ell2, self.ELL, self.kmag1,
                    self.alpha_perp, self.alpha_parallel,
                    self.sigma8, self.fz, self.b1)

        elif name == "Tree_NonGaussian_Equilateral":
            return hitomipy.integrand_B_Tree_NonGaussian_Equilateral_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ell1, self.ell2, self.ELL, self.kmag1,
                    self.alpha_perp, self.alpha_parallel,
                    self.sigma8, self.fz, self.b1)

        elif name == "Tree_NonGaussian_Orthogonal":
            return hitomipy.integrand_B_Tree_NonGaussian_Orthogonal_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ell1, self.ell2, self.ELL, self.kmag1,
                    self.alpha_perp, self.alpha_parallel,
                    self.sigma8, self.fz, self.b1)

        else:
            print("select_B: ERROR")

        return 0.0
     
    
    def select_ndim(self, name):

        name_dim3 = []
        name_dim3.append("Tree")
        name_dim3.append("Tree_NoWiggle")
        name_dim3.append("Tree_BAO")
        name_dim3.append("Tree_BAO_Reconstructed")
        name_dim3.append("Tree_BAO_Template")
        name_dim3.append("Tree_NonGaussian_Local")
        name_dim3.append("Tree_NonGaussian_Equilateral")
        name_dim3.append("Tree_NonGaussian_Orthogonal")

        if name in name_dim3:
            return 3

    def check_flag_BAO(self, name, flag_BAO, sigma8_fid, fz_fid):

        flag = 0
        name_BAO = []
        name_BAO.append("Tree_BAO")
        name_BAO.append("Tree_BAO_Reconstructed")
        name_BAO.append("Tree_BAO_Template")

        name_no_BAO = []
        name_no_BAO.append("Tree")
        name_no_BAO.append("Tree_NoWiggle")
        name_no_BAO.append("Tree_NonGaussian_Local")
        name_no_BAO.append("Tree_NonGaussian_Equilateral")
        name_no_BAO.append("Tree_NonGaussian_Orthogonal")

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
        name_recon.append("Tree_BAO_Reconstructed")

        name_no_recon = []
        name_no_recon.append("Tree")
        name_no_recon.append("Tree_NoWiggle")
        name_no_recon.append("Tree_BAO")
        name_no_recon.append("Tree_BAO_Template")
        name_no_recon.append("Tree_NonGaussian_Local")
        name_no_recon.append("Tree_NonGaussian_Equilateral")
        name_no_recon.append("Tree_NonGaussian_Orthogonal")
        
        if name in name_recon:
            if flag_recon and one_over_b1_fid >= 0.0 and R >= 0.0:
                flag = 0
                self.one_over_b1_fid = one_over_b1_fid
                self.R = R 
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
        name_damping.append("Tree_BAO_Reconstructed")
        name_damping.append("Tree_BAO_Template")

        name_no_damping = []
        name_no_damping.append("Tree")
        name_no_damping.append("Tree_NoWiggle")
        name_no_damping.append("Tree_NonGaussian_Local")
        name_no_damping.append("Tree_NonGaussian_Equilateral")
        name_no_damping.append("Tree_NonGaussian_Orthogonal")
        
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

    def check_flag_PNG(self, name, flag_PNG, k_pri, mk_pri):
        
        flag = 0
        name_PNG = []
        name_PNG.append("Tree_NonGaussian_Local")
        name_PNG.append("Tree_NonGaussian_Equilateral")
        name_PNG.append("Tree_NonGaussian_Orthogonal")

        if name in name_PNG:
            if flag_PNG and len(k_pri) > 1 and np.sum(mk_pri) > 0.0:
                flag = 0
            else:
                flag = -1
        return flag
 

    def Integrand_B(self, ndim, xx, ncomp, ff, userdata):

        self.xx_in = np.zeros(ndim[0])
        self.ff_out = np.zeros(ncomp[0])

        for i in range(ndim[0]):
            self.xx_in[i] = xx[i]
        
        self.ndim = ndim[0]
        self.ncomp = ncomp[0]

        self.select_B(self.name)

        for i in range(ncomp[0]):
            ff[i] = self.ff_out[i]

        return 0

    def calc_B(
            self, name,
            kbin=np.linspace(0.01, 0.2, 20), ell1=0, ell2=0, ELL=0,
            flag_3pcf=False,
            flag_BAO=False, sigma8_fid = - 1.0, fz_fid = - 1.0,
            flag_recon = False, one_over_b1_fid = - 1.0, R = - 1.0,
            flag_damping = False, sigma2_perp = -1.0, sigma2_para = -1.0,
            flag_PNG = False, k_pri=np.zeros(1), mk_pri=np.zeros(1),
            param_name = None):

        (kbin2_out, kbin1_out) = np.meshgrid(kbin, kbin)
        ## flags ##
        output_dict_ini = {
                "kbin1": kbin1_out,
                "kbin2": kbin2_out,
                "B": np.zeros((len(kbin), len(kbin))),
                "kbin1_fft": kbin1_out,
                "kbin2_fft": kbin2_out,
                "B_fft": np.zeros((len(kbin), len(kbin))),
                "ell1": ell1,
                "ell2": ell2,
                "ELL" : ELL,
                "flag_3pcf": flag_3pcf,
                "flag_BAO": flag_BAO,
                "flag_recon": flag_recon,
                "sigma2_perp" : -1.0,
                "sigma2_para" : -1.0
                }

        check_bao = self.check_flag_BAO(name, flag_BAO, sigma8_fid, fz_fid)
        check_recon = self.check_flag_recon(name, flag_BAO, flag_recon, one_over_b1_fid, R)
        check_damping = self.check_flag_damping(name, flag_BAO, flag_damping, sigma2_perp, sigma2_para)
        check_png = self.check_flag_PNG(name, flag_PNG, k_pri, mk_pri)

        if check_bao < 0:
            print("FLAG_BAO: ERROR")
            return output_dict_ini

        if check_recon < 0:
            print("FLAG_RECON: ERROR")
            return output_dict_ini

        if check_damping < 0:
            print("FLAG_DAMPING: ERROR")
            return output_dict_ini

        if check_png < 0:
            print("FLAG_PNG: ERROR")
            return output_dict_ini
 
        ##################
        self.name = name
        self.param_name = param_name
        
        ## set kbin ##
        if not flag_3pcf:
            self.kbin = kbin
        elif flag_3pcf:
            kbin0 = np.logspace(np.log(3.0e-4), np.log(0.2), 100, base=np.e)
            kbin1 = np.logspace(np.log(0.201), np.log(10.0), 50, base=np.e)
            self.kbin = np.hstack([kbin0, kbin1])
        
        ## set multipole indices ##
        self.ell1 = ell1
        self.ell2 = ell2
        self.ELL = ELL
        
        ## initialization ##
        hitomipy.initializeInputPowerSpectrum_py()
        
        ## calc. wigner 3j symbols ##
        hitomipy.setWigner3j_py()
        
        ## read linear power spectrum ##
        hitomipy.readInputPowerSpectrum_py(
            self.k_temp, self.P_temp, len(self.k_temp))

        hitomipy.readInputNoWigglePowerSpectrum_py(
            self.k_temp, self.P_nw_temp, len(self.k_temp))

        if flag_PNG:
            hitomipy.readInputTransferFunctionM_py(k_pri, mk_pri, len(k_pri))
        
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

        elif flag_BAO == False:
            self.sigma2_perp = sigma2_perp
            self.sigma2_para = sigma2_para

        ## compute bispectra ##
        NDIM = self.select_ndim(self.name)
        NCOMP = len(self.kbin)
        if NCOMP > 1024:
            print("# of NCOMP should be <= 1024, otherwise results become zero.")
            return output_dict_ini

        AA = []
        for i in range(NCOMP):
            print("k1 = %.6f [h/Mpc]" % self.kbin[i])
            self.kmag1 = self.kbin[i]
            if NDIM > 3:
                NNEW = 5000
                NMIN = 2
                FLATNESS = 50
                MAXEVAL = 5000
                AA.append(pycuba.Suave(self.Integrand_B, NDIM, NNEW, NMIN, FLATNESS, ncomp=NCOMP, maxeval = MAXEVAL, verbose=0 | 4)["results"])
            else:
                AA.append(pycuba.Cuhre(self.Integrand_B, NDIM, ncomp=NCOMP, key=0, verbose=0 | 4)["results"])

        bk_temp = np.zeros((NCOMP,NCOMP))
        for i in range(NCOMP):
            for j in range(NCOMP):
                bk_temp[i,j] = AA[i][j]["integral"]

        ## finalize parameters ##
        hitomipy.finalizeInputPowerSpectrum_py()

        #########################
        if flag_3pcf:
            f_bk = interpolate.interp2d(self.kbin, self.kbin, bk_temp, kind="cubic")
            (kbin2_out, kbin1_out) = np.meshgrid(kbin, kbin)
            bk_out = f_bk(kbin, kbin)
            (kbin2_fft, kbin1_fft) = np.meshgrid(self.kbin, self.kbin)
            bk_fft = bk_temp
        else:
            (kbin2_out, kbin1_out) = np.meshgrid(self.kbin, self.kbin)
            bk_out = bk_temp
            (kbin2_fft, kbin1_fft) = np.meshgrid(self.kbin, self.kbin)
            bk_fft = bk_temp

        ## output dict ##
        output_dict = {
                "kbin1": kbin1_out,
                "kbin2": kbin2_out,
                "B": bk_out,
                "kbin1_fft": kbin1_fft,
                "kbin2_fft": kbin2_fft,
                "B_fft": bk_fft,
                "ell1": ell1,
                "ell2": ell2,
                "ELL" : ELL,
                "flag_3pcf": flag_3pcf,
                "flag_BAO": flag_BAO,
                "flag_recon": flag_recon,
                "sigma2_perp": self.sigma2_perp,
                "sigma2_para": self.sigma2_para,
        }

        return output_dict

    def calc_B_to_3PCF(self, bk_in, rbin = np.linspace(0.0, 200, 41), N_fftlog = 1000):

        if not bk_in["flag_3pcf"]:
            (rbin2_out, rbin1_out) = np.meshgrid(rbin, rbin)
            output_dict_ini = {
                    "rbin1": rbin1_out, 
                    "rbin2": rbin2_out, 
                    "3pcf": np.zeros((len(rbin),len(rbin))), 
                    "rbin1_fft": rbin1_out, 
                    "rbin2_fft": rbin2_out, 
                    "3pcf_fft": np.zeros((len(rbin),len(rbin))), 
                    "ell1": bk_in["ell1"],
                    "ell2": bk_in["ell2"],
                    "ELL":  bk_in["ELL"],
                    "flag_3pcf": bk_in["flag_3pcf"],
                    "flag_BAO": bk_in["flag_BAO"],
                    "flag_recon": bk_in["flag_recon"],
                    "N_fftlog": N_fftlog}
           
            return output_dict_ini

        ## set bispec. ##
        BB = bk_in["B_fft"]

        ## set kbin ##
        self.kbin = bk_in["kbin1_fft"][:,0]
        
        ## set multipole indices ##
        self.ell1  = bk_in["ell1"]
        self.ell2  = bk_in["ell2"]
        self.ELL   = bk_in["ELL"]
        
        ## set rbin ##
        self.rbin = rbin
        
        ## input parameter for fftlog ##
        NNN = N_fftlog
        
        ## compute 3PCF ##
        kbin_for_zeta = np.logspace(np.log(self.kbin[0]), np.log(self.kbin[-1]), NNN, base=np.e)
        CC = np.zeros((len(self.kbin), NNN))
        for i in range(len(self.kbin)):

            f_bk = interpolate.interp1d(self.kbin, BB[i,:], fill_value = "extrapolate", kind="cubic")
            bk_for_zeta = f_bk(kbin_for_zeta)
            
            r_temp = np.zeros(NNN)
            zeta_temp = np.zeros(NNN)
            hitomipy.hankel_py(self.ell2, 2, NNN, kbin_for_zeta, bk_for_zeta, r_temp, zeta_temp)
            
            f_zeta = interpolate.interp1d(r_temp, zeta_temp, fill_value = "extrapolate", kind="cubic")
            CC[i,:] = f_zeta(r_temp[:])

        DD = np.zeros((NNN, NNN))   
        for j in range(NNN):

            f_bk = interpolate.interp1d(self.kbin, CC[:,j], fill_value = "extrapolate", kind="cubic")
            bk_for_zeta = f_bk(kbin_for_zeta)
            
            r_temp = np.zeros(NNN)
            zeta_temp = np.zeros(NNN)
            hitomipy.hankel_py(self.ell1, 2, NNN, kbin_for_zeta, bk_for_zeta, r_temp, zeta_temp)
            
            f_zeta = interpolate.interp1d(r_temp, zeta_temp, fill_value = "extrapolate", kind="cubic")
            DD[:,j] = f_zeta(r_temp[:])

        sign = np.real(1.0j**(self.ell1+self.ell2))
        zeta_fft = sign * DD
        
        f_zeta = interpolate.interp2d(r_temp, r_temp, zeta_fft, kind="cubic")
        zeta_out = np.zeros((len(self.rbin),len(self.rbin)))
        zeta_out[:,:] = f_zeta(self.rbin[:], self.rbin[:])

        ## rbin_out ##
        (rbin2_out, rbin1_out) = np.meshgrid(self.rbin, self.rbin)
        (rbin2_fft, rbin1_fft) = np.meshgrid(r_temp, r_temp)

        ## output dict ##
        output_dict = {
                "rbin1": rbin1_out, 
                "rbin2": rbin2_out, 
                "3pcf": zeta_out, 
                "rbin1_fft": rbin1_fft, 
                "rbin2_fft": rbin2_fft, 
                "3pcf_fft": zeta_fft, 
                "ell1": bk_in["ell1"],
                "ell2": bk_in["ell2"],
                "ELL":  bk_in["ELL"],
                "flag_3pcf": bk_in["flag_3pcf"],
                "flag_BAO": bk_in["flag_BAO"],
                "flag_recon": bk_in["flag_recon"],
                "N_fftlog": N_fftlog}
 
        return output_dict

    def calc_3PCF_to_B(self, zeta_in, kbin = np.linspace(0.01, 0.2, 20)):

        if not zeta_in["flag_3pcf"]:
            (kbin2_out, kbin1_out) = np.meshgrid(kbin, kbin)
            output_dict_ini = {
                    "kbin1": kbin1_out,
                    "kbin2": kbin2_out,
                    "B": np.zeros((len(kbin),len(kbin))),
                    "ell1": zeta_in["ell1"],
                    "ell2": zeta_in["ell2"],
                    "ELL": zeta_in["ELL"],
                    "flag_3pcf": zeta_in["flag_3pcf"],
                    "flag_BAO": zeta_in["flag_BAO"],
                    "flag_recon": zeta_in["flag_recon"],
                    "N_fftlog": zeta_in["N_fftlog"]}
    
            return output_dict_ini

        ## set 3pcf ##
        BB = zeta_in["3pcf_fft"]

        ## set rbin ##
        self.rbin = zeta_in["rbin1_fft"][:,0]
        
        ## set multipole indices ##
        self.ell1 = zeta_in["ell1"]
        self.ell2 = zeta_in["ell2"]
        self.ELL  = zeta_in["ELL"]
        
        ## set kbin ##
        self.kbin = kbin

        ## input parameter for fftlog ##
        NNN = zeta_in["N_fftlog"]
        
        ## compute bispec ##
        rbin_for_bk = np.logspace(np.log(self.rbin[0]), np.log(self.rbin[-1]), NNN, base=np.e)
        CC = np.zeros((len(self.rbin), len(self.kbin)))
        for i in range(len(self.rbin)):

            f_zeta = interpolate.interp1d(self.rbin, BB[i,:], fill_value = "extrapolate", kind="cubic")
            zeta_for_bk = f_zeta(rbin_for_bk)
            
            k_temp = np.zeros(NNN)
            bk_temp = np.zeros(NNN)
            hitomipy.hankel_py(self.ell2, 2, NNN, rbin_for_bk, zeta_for_bk, k_temp, bk_temp)
            
            f_bk = interpolate.interp1d(k_temp, bk_temp, fill_value = "extrapolate", kind="cubic")
            CC[i,:] = f_bk(self.kbin[:])

        DD = np.zeros((len(self.kbin), len(self.kbin)))
        for j in range(len(self.kbin)):

            f_zeta = interpolate.interp1d(self.rbin, CC[:,j], fill_value = "extrapolate", kind="cubic")
            zeta_for_bk = f_zeta(rbin_for_bk)
            
            k_temp = np.zeros(NNN)
            bk_temp = np.zeros(NNN)
            hitomipy.hankel_py(self.ell1, 2, NNN, rbin_for_bk, zeta_for_bk, k_temp, bk_temp)
            
            f_bk = interpolate.interp1d(k_temp, bk_temp, fill_value = "extrapolate", kind="cubic")
            DD[:,j] = f_bk(self.kbin[:])

        sign = np.real((-1.0j)**(self.ell1+self.ell2)) * (2.0*np.pi)**6
        bk_out = sign * DD
    	
        ## rbin_out ##
        (kbin2_out, kbin1_out) = np.meshgrid(self.kbin, self.kbin)
        
        ## output dict ##
        output_dict = {
                "kbin1": kbin1_out, 
                "kbin2": kbin2_out, 
                "B": bk_out, 
                "ell1": self.ell1,
                "ell2": self.ell2,
                "ELL": self.ELL,
                "flag_3pcf": zeta_in["flag_3pcf"],
                "flag_BAO": zeta_in["flag_BAO"],
                "flag_recon": zeta_in["flag_recon"],
                "N_fftlog": zeta_in["N_fftlog"]}
                        
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


