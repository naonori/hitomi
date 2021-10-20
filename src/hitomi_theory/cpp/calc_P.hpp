#ifndef __calc_P__
#define __calc_P__

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

#ifndef __P__
#include "P.hpp"
#endif

int integrand_P_Tree(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_kbin, int ELL, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Tree_NoWiggle(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_kbin, int ELL, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
         
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_NoWiggle(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Tree_BAO(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_kbin, int ELL, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double sigma2_perp, double sigma2_para) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_BAO(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Tree_BAO_Template(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_kbin, int ELL, 
        double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para,
        char * parameters) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_BAO_Template(kvec, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para, parameters);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}
	return 0;
}

int integrand_P_sigma2_perp_Reconstructed(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_kbin, 
        double sigma8, double fz, double b1,
        double one_over_b1_fid, double R) {

    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double kmag = exp( log(kmin) + dlnk * xx_in[0] );
    double mu = -1.0 + 2.0 * xx_in[1];
    mu = mu;
    double kmag2 = kmag * kmag;
    
    double W = exp( - pow(kmag * R, 2) / 2.0 );
    double result = (1.0/3.0) * (pow(sigma8, 2) * f_pk(kmag) / kmag2) * 
                   ( ( 1.0 - (W*b1 * one_over_b1_fid) ) * ( 1.0 - (W*b1 * one_over_b1_fid) )
                   + ( 1.0 - (W*b1 * one_over_b1_fid) ) * (-W*fz * one_over_b1_fid) * (2.0/5.0)
                   + (-W*fz * one_over_b1_fid) * (-W*fz * one_over_b1_fid) * (3.0/35.0) );
    
    double jacobian = dlnk * pow(kmag,3) / (2.0 * M_PI * M_PI);
    ff_out[0] = result * jacobian;
    
    return 0;
}

int integrand_P_sigma2_para_Reconstructed(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_kbin, 
        double sigma8, double fz, double b1,
        double one_over_b1_fid, double R) {

    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double kmag = exp( log(kmin) + dlnk * xx_in[0] );
    double mu = -1.0 + 2.0 * xx_in[1];
    mu = mu;
    double kmag2 = kmag * kmag;
    
    double W = exp( - pow(kmag * R, 2) / 2.0 );
    double result = (1.0/3.0) * (pow(sigma8, 2) * f_pk(kmag) / kmag2) * 
                   ( ( ( 1.0 - (W*b1 * one_over_b1_fid) ) + fz ) * ( ( 1.0 - (W*b1 * one_over_b1_fid) ) + fz )
                   + ( 1.0 - (W*b1 * one_over_b1_fid) ) * (-W*fz * one_over_b1_fid) * (6.0/5.0)
                   + (-W*fz * one_over_b1_fid) * (6.0 * fz / 5.0)
                   + (-W*fz * one_over_b1_fid) * (-W*fz * one_over_b1_fid) * (3.0/7.0)
                   );

    double jacobian = dlnk * pow(kmag,3) / (2.0 * M_PI * M_PI);
    ff_out[0] = result * jacobian;
    
    return 0;

}

#endif

