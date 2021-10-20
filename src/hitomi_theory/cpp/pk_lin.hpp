#ifndef __pk_lin__
#define __pk_lin__

#ifndef __common__
#include "common.hpp"
#endif

gsl_interp_accel * acc_pk = NULL;
gsl_spline       * spline_pk = NULL;

gsl_interp_accel * acc_pk_nw = NULL;
gsl_spline       * spline_pk_nw = NULL;


gsl_interp_accel * acc_pk_pri = NULL;
gsl_spline       * spline_pk_pri = NULL;


gsl_interp_accel * acc_pk_nl_0 = NULL;
gsl_spline       * spline_pk_nl_0 = NULL;
gsl_interp_accel * acc_pk_nl_2 = NULL;
gsl_spline       * spline_pk_nl_2 = NULL;

double kmax_nl = 0.0;
double kmin_nl = 0.0;

double pk_kmin = 0.0;
double pk_kmax = 0.0;
double norm = 1.0;
double norm_no_wiggle = 1.0;
double norm_pri = 1.0;

double get_norm() {
    return norm;
}

int set_kmin(double kmin) {
    if(kmin > pk_kmin) {
	pk_kmin = kmin;
    }
    return 0;
}

int set_kmax(double kmax) {
    if(kmax < pk_kmax) {
	pk_kmax = kmax;
    }
    return 0;
}

double f_pk(double kmag) {

	if( (kmag < pk_kmin) || (kmag > pk_kmax) ) { 
		return 0.0;
	} else if (spline_pk == NULL || acc_pk == NULL) {
		return 0.0;
	} else {
		return norm * gsl_spline_eval(spline_pk, kmag, acc_pk);	  
	}
}

double f_pk_no_wiggle(double kmag) {

	if( (kmag < pk_kmin) || (kmag > pk_kmax) ) { 
		return 0.0;
	} else if (spline_pk_nw == NULL || acc_pk_nw == NULL) {
		return 0.0;
	} else {
		return norm_no_wiggle * gsl_spline_eval(spline_pk_nw, kmag, acc_pk_nw);	  
	}

}


double f_pk_nl_0(double kmag) {

	if( (kmag < kmin_nl) || (kmag > kmax_nl) ) { 
		return 0.0;
	} else if (spline_pk_nl_0 == NULL || acc_pk_nl_0 == NULL) {
		return 0.0;
	} else {
		return gsl_spline_eval(spline_pk_nl_0, kmag, acc_pk_nl_0);	  
	}
}

double f_pk_nl_2(double kmag) {

	if( (kmag < kmin_nl) || (kmag > kmax_nl) ) { 
		return 0.0;
	} else if (spline_pk_nl_2 == NULL || acc_pk_nl_2 == NULL) {
		return 0.0;
	} else {
		return gsl_spline_eval(spline_pk_nl_2, kmag, acc_pk_nl_2);	  
	}
}

double f_Mk_pri(double kmag) {

	if( (kmag < pk_kmin) || (kmag > pk_kmax) ) { 
		return 0.0;
	} else if (spline_pk_pri == NULL || acc_pk_pri == NULL) {
		return 0.0;
	} else {
		return norm_pri * gsl_spline_eval(spline_pk_pri, kmag, acc_pk_pri);	  
	}
}

int readInputPowerSpectrum(double * kmag_in, double * pk_in, int pk_num_in) {

	if(kmag_in[0] <= 0.0) {
		return -1;
	}

//	/* set kmin and kmax */
	pk_kmin = kmag_in[0];
	pk_kmax = kmag_in[pk_num_in-1];

	/* spline interpolation */
	acc_pk = gsl_interp_accel_alloc();
	spline_pk = gsl_spline_alloc(gsl_interp_cspline, pk_num_in);
	gsl_spline_init(spline_pk, kmag_in, pk_in, pk_num_in);

	return 0;
}

int readInputNoWigglePowerSpectrum(double * kmag_in, double * pk_in, int pk_num_in) {

	if(kmag_in[0] <= 0.0) {
		return -1;
	}

//	/* set kmin and kmax */
	pk_kmin = kmag_in[0];
	pk_kmax = kmag_in[pk_num_in-1];

	/* spline interpolation */
	acc_pk_nw = gsl_interp_accel_alloc();
	spline_pk_nw = gsl_spline_alloc(gsl_interp_cspline, pk_num_in);
	gsl_spline_init(spline_pk_nw, kmag_in, pk_in, pk_num_in);

	return 0;
}


int readNonLinearPowerSpectrum(double * kmag_in, double * pk_0_in, double * pk_2_in, int pk_num_in) {

	if(kmag_in[0] <= 0.0) {
		return -1;
	}

//	/* set kmin and kmax */
	kmin_nl = kmag_in[0];
	kmax_nl = kmag_in[pk_num_in-1];

	/* spline interpolation */
	acc_pk_nl_0 = gsl_interp_accel_alloc();
	spline_pk_nl_0 = gsl_spline_alloc(gsl_interp_cspline, pk_num_in);
	gsl_spline_init(spline_pk_nl_0, kmag_in, pk_0_in, pk_num_in);

	acc_pk_nl_2 = gsl_interp_accel_alloc();
	spline_pk_nl_2 = gsl_spline_alloc(gsl_interp_cspline, pk_num_in);
	gsl_spline_init(spline_pk_nl_2, kmag_in, pk_2_in, pk_num_in);

	return 0;
}


void calcNormalizationUsingSigma8(double sigma8) {
	norm = 1.0 / (sigma8 * sigma8);
}

int readInputTransferFunctionM(double * kmag_in, double * pk_in, int pk_num_in) {

	if(kmag_in[0] <= 0.0) {
		return -1;
	}

//	/* set kmin and kmax */
	pk_kmin = kmag_in[0];
	pk_kmax = kmag_in[pk_num_in-1];

	/* spline interpolation */
	acc_pk_pri = gsl_interp_accel_alloc();
	spline_pk_pri = gsl_spline_alloc(gsl_interp_cspline, pk_num_in);
	gsl_spline_init(spline_pk_pri, kmag_in, pk_in, pk_num_in);

	return 0;
}

double calcNoWiggleMatterPowerSpectrum(double kmag, double h, double omega0, double omegab, double Tcmb, double n_s) {

	double ss = 44.5 * h * log( 9.83 / (omega0*h*h) ) / sqrt( 1.0 + 10.0 * pow(omegab*h*h,0.75) );
	double alpha_gam = 1.0 - 0.328 * log( 431.0 * omega0*h*h ) * omegab/omega0 + 0.38 * log( 22.3 * omega0*h*h ) * pow(omegab/omega0,2);
	double theta_cmb = Tcmb / 2.70;
	double gamma_eff = omega0 * h * ( alpha_gam + (1.0 - alpha_gam) / (1.0 + pow(0.43*kmag*ss, 4)) );

	double q = kmag * pow(theta_cmb,2) / gamma_eff;
	double L0 = log( 2.0 * exp(1.0) + 1.8 * q );
	double C0 = 14.2 + 731.0 / ( 1.0 + 62.5 * q );

	double T_EH = L0 / (L0 + C0*q*q );
	
	double Pk_lin_EH = pow(kmag,n_s) * pow(T_EH,2);

	return Pk_lin_EH;

}

void calcNormalizationNoWiggle(double sigma8) {

	const int N = 100000;
	const double R = 8.0;
	double kmag, kmag2, kmag3, dlnk;
	double kr, kr2, kr3, w = 0.0, w2 = 0.0;

	double integ = 0.0;
	for(int n = 0; n <= N; n++) {
		dlnk = log(pk_kmax/pk_kmin) / double(N);
		kmag = pk_kmin * exp(dlnk * double(n));
		kmag2 = kmag * kmag;
		kmag3 = kmag * kmag2;

		kr = R * kmag;
		kr2 = kr * kr;
		kr3 = kr2 * kr;
  		w = 3 * (sin(kr) / kr3 - cos(kr) / kr2);
		w2 = w * w;
		
		double power =  kmag3 * w2 * f_pk_no_wiggle(kmag);

		if(0) {
		} else if ( (n == 0) || (n == N) ) {
			integ += power;
		} else if ( (n % 2) == 0 ) {
			integ += 2.0 * power;
		} else if ( (n % 2) == 1 ) {
			integ += 4.0 * power;
		}
	}

	integ *= dlnk / 3.0;
	integ /= (2.0 * M_PI * M_PI);

    if(fabs(norm_no_wiggle - 1.0) < 1.0e-10) {
	    norm_no_wiggle = sigma8 * sigma8 / integ;
    }

}

double calcSigma_dd(double sigma8) {

	const int N = 100000;
	double kmag, dlnk;

	double integ = 0.0;
	for(int n = 0; n <= N; n++) {
		dlnk = log(pk_kmax/pk_kmin) / double(N);
		kmag = pk_kmin * exp(dlnk * double(n));
		
		double power =  kmag * f_pk(kmag);

		if(0) {
		} else if ( (n == 0) || (n == N) ) {
			integ += power;
		} else if ( (n % 2) == 0 ) {
			integ += 2.0 * power;
		} else if ( (n % 2) == 1 ) {
			integ += 4.0 * power;
		}
	}

	integ *= dlnk / 3.0;
	integ /= (2.0 * M_PI * M_PI);
	
	return sigma8 * sigma8 * integ / 3.0;

}



void initializeInputPowerSpectrum() {

	acc_pk = NULL;
	spline_pk= NULL;
	
    acc_pk_nw = NULL;
	spline_pk_nw = NULL;
	
	acc_pk_pri = NULL;
	spline_pk_pri= NULL;

	acc_pk_nl_0 = NULL;
	spline_pk_nl_0 = NULL;
	acc_pk_nl_2 = NULL;
	spline_pk_nl_2 = NULL;
	
	pk_kmin = 0.0;
	pk_kmax = 0.0;
	norm = 1.0;
	norm_pri = 1.0;
	norm_no_wiggle = 1.0;

	kmin_nl = 0.0;
	kmax_nl = 0.0;

}


void finalizeInputPowerSpectrum() {
	gsl_interp_accel_free(acc_pk);
	acc_pk = NULL;
	gsl_spline_free(spline_pk);
	spline_pk = NULL;

    gsl_interp_accel_free(acc_pk_nw);
	acc_pk_nw = NULL;
	gsl_spline_free(spline_pk_nw);
	spline_pk_nw = NULL;

	gsl_interp_accel_free(acc_pk_pri);
	acc_pk_pri = NULL;
	gsl_spline_free(spline_pk_pri);
	spline_pk_pri = NULL;

	norm = 1.0;
	norm_pri = 1.0;
	norm_no_wiggle = 1.0;

	gsl_interp_accel_free(acc_pk_nl_0);
	acc_pk_nl_0 = NULL;
	gsl_spline_free(spline_pk_nl_0);
	spline_pk_nl_0 = NULL;

	gsl_interp_accel_free(acc_pk_nl_2);
	acc_pk_nl_2 = NULL;
	gsl_spline_free(spline_pk_nl_2);
	spline_pk_nl_2 = NULL;

}



#endif

