# cython: c_string_type=unicode, c_string_encoding=utf8

import numpy as np
cimport numpy as np

# distutils: language=["c++"]
# distutils: sources=["cpp/pyfftlog.hpp"]
cdef extern from "cpp/pyfftlog.hpp":
    void fftlog_ComputeXiLM(int l, int m, int N, double * k, double * pk, double * r, double * xi)

def hankel_py(
        int l, int m, int N,
        np.ndarray[double, ndim=1, mode="c"] k not None,
        np.ndarray[double, ndim=1, mode="c"] pk not None,
        np.ndarray[double, ndim=1, mode="c"] r not None,
        np.ndarray[double, ndim=1, mode="c"] xi not None
        ):
    return fftlog_ComputeXiLM(l, m, N, & k[0], & pk[0], & r[0], & xi[0])

# distutils: sources=["cpp/wigner.hpp"]
cdef extern from "cpp/wigner.hpp":
    int setWigner3j()

def setWigner3j_py():
    return setWigner3j()

# distutils: sources=["cpp/pk_lin.hpp"]
cdef extern from "cpp/pk_lin.hpp":
    int readInputPowerSpectrum(double * kbin_in, double * pk_in, int pk_num_in)
    int readInputNoWigglePowerSpectrum(double * kbin_in, double * pk_in, int pk_num_in)
    int readInputTransferFunctionM(double * kbin_in, double * pk_in, int pk_num_in)

    void initializeInputPowerSpectrum()
    void finalizeInputPowerSpectrum()
    void calcNormalizationUsingSigma8(double sigma8)
    void calcNormalizationNoWiggle(double sigma8)
    double set_kmin(double kmin)
    double set_kmax(double kmax)
    double calcSigma_dd(double sigma8)
    double calcNoWiggleMatterPowerSpectrum(double kmag, double h, double omega0, double omegab, double Tcmb, double n_s)

def calcSigma_dd_py(double sigma8):
    return calcSigma_dd(sigma8)

def calcNoWiggleMatterPowerSpectrum_py(
        double kmag, double h, double omega0, double omegab, double Tcmb, double n_s):
    return calcNoWiggleMatterPowerSpectrum(kmag, h, omega0, omegab, Tcmb, n_s)

def readInputPowerSpectrum_py(
        np.ndarray[double, ndim=1, mode="c"] kbin_in not None,
        np.ndarray[double, ndim=1, mode="c"] pk_in not None,
        int pk_num_in
        ):
    return readInputPowerSpectrum( & kbin_in[0], & pk_in[0], pk_num_in)

def readInputNoWigglePowerSpectrum_py(
        np.ndarray[double, ndim=1, mode="c"] kbin_in not None,
        np.ndarray[double, ndim=1, mode="c"] pk_in not None,
        int pk_num_in
        ):
    return readInputNoWigglePowerSpectrum( & kbin_in[0], & pk_in[0], pk_num_in)

def readInputTransferFunctionM_py(
        np.ndarray[double, ndim=1, mode="c"] kbin_in not None,
        np.ndarray[double, ndim=1, mode="c"] pk_in not None,
        int pk_num_in
        ):
    return readInputTransferFunctionM( & kbin_in[0], & pk_in[0], pk_num_in)

def initializeInputPowerSpectrum_py():
    return initializeInputPowerSpectrum()

def finalizeInputPowerSpectrum_py():
    return finalizeInputPowerSpectrum()

def calcNormalizationUsingSigma8_py(double sigma8):
    return calcNormalizationUsingSigma8(sigma8)

def calcNormalizationNoWiggle_py(double sigma8):
    return calcNormalizationNoWiggle(sigma8)

def set_kmin_py(double kmin):
    return set_kmin(kmin)

def set_kmax_py(double kmax):
    return set_kmax(kmax)

# distutils: sources=["cpp/calc_P.hpp"]
cdef extern from "cpp/calc_P.hpp":

    int integrand_P_Tree(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel,
            double sigma8, double fz, double b1)

    int integrand_P_Tree_NoWiggle(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, 
            double sigma8, double fz, double b1)

    int integrand_P_Tree_BAO(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel,
            double sigma8, double fz, double b1,
            double sigma2_perp, double sigma2_para)

    int integrand_P_Tree_BAO_Template(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para,
            char * param_name)

    int integrand_P_sigma2_perp_Reconstructed(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_kbin, 
            double sigma8, double fz, double b1,
            double one_over_b1_fid, double R)

    int integrand_P_sigma2_para_Reconstructed(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_kbin, 
            double sigma8, double fz, double b1,
            double one_over_b1_fid, double R)

def integrand_P_Tree_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel,
        double sigma8, double fz, double b1):
    return integrand_P_Tree(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            alpha_perp, alpha_parallel, sigma8, fz, b1)

def integrand_P_Tree_NoWiggle_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, 
        double alpha_perp, double alpha_parallel,
        double sigma8, double fz, double b1):
    return integrand_P_Tree_NoWiggle(
            &xx_in[0], ndim, & ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, 
            alpha_perp, alpha_parallel, sigma8, fz, b1)

def integrand_P_Tree_BAO_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, 
        double alpha_perp, double alpha_parallel, 
        double sigma8, double fz, double b1,
        double sigma2_perp, double sigma2_para):
    return integrand_P_Tree_BAO( 
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            sigma2_perp, sigma2_para)

def integrand_P_Tree_BAO_Template_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel,
        double sigma2_perp, double sigma2_para,
        char * param_name):
    return integrand_P_Tree_BAO_Template( 
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, 
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para,
            &param_name[0])

def integrand_P_sigma2_perp_Reconstructed_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, 
        double sigma8, double fz, double b1,
        double one_over_b1_fid, double R):
    return integrand_P_sigma2_perp_Reconstructed(
            &xx_in[0], ndim, & ff_out[0], ncomp, 
            &kbin[0], num_k_bin, 
            sigma8, fz, b1,
            one_over_b1_fid, R)

def integrand_P_sigma2_para_Reconstructed_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, 
        double sigma8, double fz, double b1,
        double one_over_b1_fid, double R):
    return integrand_P_sigma2_para_Reconstructed(
            &xx_in[0], ndim, & ff_out[0], ncomp, 
            &kbin[0], num_k_bin, 
            sigma8, fz, b1,
            one_over_b1_fid, R)


# distutils: sources=["cpp/calc_B.hpp"]
cdef extern from "cpp/calc_B.hpp":

    int integrand_B_Tree(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, 
            double sigma8, double fz, double b1, 
            double b2, double bK2)

    int integrand_B_Tree_NoWiggle(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, 
            double sigma8, double fz, double b1,
            double b2, double bK2)

    int integrand_B_Tree_BAO(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel,
            double sigma8, double fz, double b1, 
            double b2, double bK2,
            double sigma2_perp, double sigma2_para)

    int integrand_B_Tree_BAO_Reconstructed(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel,
            double sigma8, double fz, double b1, 
            double b2, double bK2,
            double one_over_b1_fid, double R, 
            double sigma2_perp, double sigma2_para)

    int integrand_B_Tree_BAO_Template(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para,
            char * param_name)

    int integrand_B_Tree_NonGaussian_Local(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel,
            double sigma8, double fz, double b1)

    int integrand_B_Tree_NonGaussian_Equilateral(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, 
            double sigma8, double fz, double b1)

    int integrand_B_Tree_NonGaussian_Orthogonal(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel,
            double sigma8, double fz, double b1)


def integrand_B_Tree_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma8, double fz, double b1, 
        double b2, double bK2):
    return integrand_B_Tree(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma8, fz, b1, 
            b2, bK2)

def integrand_B_Tree_NoWiggle_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel,
        double sigma8, double fz, double b1,
        double b2, double bK2):
    return integrand_B_Tree_NoWiggle(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel, 
            sigma8, fz, b1, 
            b2, bK2)

def integrand_B_Tree_BAO_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma8, double fz, double b1, 
        double b2, double bK2, 
        double sigma2_perp, double sigma2_para):
    return integrand_B_Tree_BAO(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma8, fz, b1,
            b2, bK2,
            sigma2_perp, sigma2_para)


def integrand_B_Tree_BAO_Reconstructed_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel,
        double sigma8, double fz, double b1, 
        double b2, double bK2, 
        double one_over_b1_fid, double R,
        double sigma2_perp, double sigma2_para):
    return integrand_B_Tree_BAO_Reconstructed(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma8, fz, b1,
            b2, bK2,
            one_over_b1_fid, R,
            sigma2_perp, sigma2_para)

def integrand_B_Tree_BAO_Template_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para,
        char * param_name):
    return integrand_B_Tree_BAO_Template(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para,
            &param_name[0])

def integrand_B_Tree_NonGaussian_Local_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma8, double fz, double b1):
    return integrand_B_Tree_NonGaussian_Local(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma8, fz, b1)

def integrand_B_Tree_NonGaussian_Equilateral_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma8, double fz, double b1):
    return integrand_B_Tree_NonGaussian_Equilateral(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma8, fz, b1)


def integrand_B_Tree_NonGaussian_Orthogonal_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel,
        double sigma8, double fz, double b1):
    return integrand_B_Tree_NonGaussian_Orthogonal(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma8, fz, b1)

#
#def integrand_SS_py(
#        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#        int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
#        int n, int m,
#        np.ndarray[double, ndim=1, mode="c"] epsilon not None, int num_epsilon
#        ):
#    return integrand_SS(
#            &xx_in[0], ndim, &ff_out[0], ncomp, 
#            ell1, ell2, ELL, ell1_dash, ell2_dash, ELL_dash,
#            n, m,
#            &epsilon[0], num_epsilon)
#
#
