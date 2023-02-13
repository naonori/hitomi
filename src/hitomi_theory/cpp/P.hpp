#ifndef __P__
#define __P__

#ifndef __common__
#include "common.hpp"
#endif

#ifndef __pk_lin__
#include "pk_lin.hpp"
#endif

#ifndef __sph__
#include "sph.hpp"
#endif

#ifndef __kernel__
#include "kernel.hpp"
#endif

double Powerspectrum_Tree(
        double * kvec_in, double * los, 
        double alpha_perp, double alpha_parallel,
        double sigma8, double f, double b1) {
	
    double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
	double P = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2) * f_pk(k);
	return P / alpha3;

}

double Powerspectrum_Tree_BAO(
        double * kvec_in, double * los,
        double alpha_perp, double alpha_parallel,
        double sigma8, double f, double b1,
	    double sigma2_perp, double sigma2_para) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
    double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);
    double BAO = f_pk(k) - f_pk_no_wiggle(k);
	double Kaiser = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2);
    double G = D * D * Kaiser * BAO;
	double MC = Kaiser * f_pk_no_wiggle(k);
	return (G+MC) / alpha3;

}

double Powerspectrum_Tree_NoWiggle(
        double * kvec_in, double * los,
        double alpha_perp, double alpha_parallel,
        double sigma8, double f, double b1) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
	double Kaiser = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2);
	double P = Kaiser * f_pk_no_wiggle(k);
	return P / alpha3;

}

double PowerspectrumKernelFunctionTemplate(
        double * kvec_in, double * los, char * param_name_in) {

    std::string param_name = param_name_in;

    if(0) {
    } else if (param_name == "b1_b1") {
        return D1() * D1();
    } else if (param_name == "b1_f") {
        return 2.0 * D1() * V1(kvec_in, los);
    } else if (param_name == "f_f") {
        return V1(kvec_in, los) * V1(kvec_in, los);
    } else {
        return 0.0;
    }

}

double Powerspectrum_Tree_BAO_Template(
        double * kvec_in, double * los, double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para,
        char * param_name) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
    double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);
     	
	double BAO = f_pk(k) - f_pk_no_wiggle(k);

	double Kaiser = PowerspectrumKernelFunctionTemplate(kvec, los, param_name);

    double G = D * D * Kaiser * BAO;
	double MC = Kaiser * f_pk_no_wiggle(k);
	return (G+MC) / alpha3;

}

double Powerspectrum_Tree_Reconstructed_Correction(
        double * kvec_in, double * los,
        double alpha_perp, double alpha_parallel,
        double sigma8, double f, double b1,
	    double sigma2_perp, double sigma2_para,
        double one_over_b1_fid, double R
        ) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
    double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);
    double W = exp( - pow(k * R, 2) / 2.0 );
	double Kaiser = 2.0 * (1.0 - D) * (1.0 - one_over_b1_fid * W) * (- one_over_b1_fid * W) * Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2);
	double MC = Kaiser * f_pk_no_wiggle(k);
	return (MC) / alpha3;

}

#endif

