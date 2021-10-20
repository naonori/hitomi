#ifndef __calc_B__
#define __calc_B__

#ifndef __common__
#include "common.hpp"
#endif

#ifndef __pk_lin__
#include "pk_lin.hpp"
#endif

#ifndef __kernel__
#include "kernel.hpp"
#endif

#ifndef __sph__
#include "sph.hpp"
#endif

#ifndef __B__
#include "B.hpp"
#endif

int integrand_B_Tree(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
        double alpha_perp, double alpha_parallel,
        double sigma8, double fz, double b1, 
        double b2, double bK2) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];

	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_Tree_NoWiggle(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel,
        double sigma8, double fz, double b1, 
        double b2, double bK2) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_Tree_BAO(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel,
        double sigma8, double fz, double b1, 
        double b2, double bK2,
        double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_Tree_BAO_Reconstructed(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
        double alpha_perp, double alpha_parallel,
        double sigma8, double fz, double b1, 
        double b2, double bK2,
        double one_over_b1_fid, double R, 
        double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_Reconstructed(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, one_over_b1_fid, R, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_Tree_BAO_Template(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para,
        char * param_name) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    double bispec = Bispectrum_Tree_BAO_Template(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para, param_name);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}

	return 0;

}

int integrand_B_Tree_NonGaussian_Local(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
        double alpha_perp, double alpha_parallel,
        double sigma8, double fz, double b1) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NonGaussian_Local(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_Tree_NonGaussian_Equilateral(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
        double alpha_perp, double alpha_parallel, 
        double sigma8, double fz, double b1) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NonGaussian_Equilateral(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}

	return 0;

}

int integrand_B_Tree_NonGaussian_Orthogonal(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
        double alpha_perp, double alpha_parallel,
        double sigma8, double fz, double b1) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NonGaussian_Orthogonal(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}




//
//
//double get_nm(int n) {
//
//    if(0) {
//    } else if (n == 0) {
//	return 1.0;
//    } else if (n == 1) {
//	return 1.0;
//    } else if (n == 2) {
//	return 2.0;
//    } else if (n == 3) {
//	return 6.0;
//    } else if (n == 4) {
//	return 24.0;
//    } else {
//	return 1.0e-10;
//    }
//
//}
//
//int integrand_SS(
//        double * xx_in, int ndim, double * ff_out, int ncomp,  
//        int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
//        int n, int m,
//        double * epsilon, int num_epsilon) {
//
//	/********/
//    double mu1  = 1.0;
//    double phi1 = 0.0;
//    double mu2  = - 1.0 + 2.0 * xx_in[0];
//    double phi2 =   0.0;
//    double mu   = - 1.0 + 2.0 * xx_in[1];      
//    double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
//	/********/
//
//    for(int i = 0; i < num_epsilon; i++) {
//
//        /********/
//    	double kvec1_hat[3] = {sqrt(1.0 - mu1 * mu1) * cos(phi1), sqrt(1.0 - mu1 * mu1) * sin(phi1), mu1};
//    	double kvec2_hat[3] = {sqrt(1.0 - mu2 * mu2) * cos(phi2), sqrt(1.0 - mu2 * mu2) * sin(phi2), mu2};
//    	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
//     
//    	double kvec1_hat_dash[3] = { 0.0, 0.0, 1.0 };
//    	double kvec2_hat_dash[3] = { 0.0, 0.0, 1.0 };
//    	calcTrueWavevectorHat(kvec1_hat, los, epsilon[i], kvec1_hat_dash);
//    	calcTrueWavevectorHat(kvec2_hat, los, epsilon[i], kvec2_hat_dash);
//    	/********/
//    
//    	/********/
//    	double mu1_n = MU(kvec1_hat, los); 
//    	double mu2_m = MU(kvec2_hat, los); 
//    	double L1 = (3.0 * mu1_n * mu1_n - 1.0) / 2.0;
//    	double L2 = (3.0 * mu2_m * mu2_m - 1.0) / 2.0;
//    
//    	double fac1_a =  (2.0/3.0) * L1 * ( pow(1.0 + epsilon[i], - 6.0) - 1.0 );
//    	double fac1_b = 1.0 + (1.0/3.0) * ( pow(1.0 + epsilon[i], - 6.0) - 1.0 );
//    	double Delta1 = pow( 1.0 + fac1_a / fac1_b, 0.5 ) - 1.0;
//    	double Delta1_n = pow( Delta1, double(n));
//    
//    	double fac2_a =  (2.0/3.0) * L2 * ( pow(1.0 + epsilon[i], - 6.0) - 1.0 );
//    	double fac2_b = 1.0 + (1.0/3.0) * ( pow(1.0 + epsilon[i], - 6.0) - 1.0 );
//    	double Delta2 = pow( 1.0 + fac2_a / fac2_b, 0.5 ) - 1.0;
//    	double Delta2_m = pow( Delta2, double(m));
//    	/********/
//    
//    	/********/
//    	double YYY_YYY = calcYYY_YYY(kvec1_hat, kvec2_hat, kvec1_hat_dash, kvec2_hat_dash, los, ell1, ell2, ELL, ell1_dash, ell2_dash, ELL_dash);
//    	/********/
//    	
//    	/********/
//    	double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
//    	double NN = get_nm(n);
//    	double MM = get_nm(m);
//    	double result = Nlll * YYY_YYY * Delta1_n * Delta2_m / (NN * MM);
//    	/********/
//    	
//    	double jacobian = 1.0;
//    	ff_out[i] = result * jacobian;
//
//    }
//
//	return 0;
//
//}
//
#endif

