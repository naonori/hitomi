#ifndef __sph__
#define __sph__

#ifndef __common__
#include "common.hpp"
#endif

#ifndef __wigner__
#include "wigner.hpp"
#endif

#define wigner_3j(j1,j2,j3,m1,m2,m3) (w3j[j1][j2][j3][m1+j1][m2+j2][m3+j3])

double LegendrePolynomials(int ell, double mu) {

    if(0) {
    } else if (ell == 0) {
	return 1.0;
    } else if (ell == 1) {
	return mu;
    } else if (ell == 2) {
	return (3.0 * mu * mu - 1.0) / 2.0;
    } else if (ell == 3) {
	return (5.0 * mu * mu * mu - 3.0 * mu) / 2.0;
    } else if (ell == 4) {
	double mu2 = mu * mu;
	return (35.0 * mu2 * mu2 - 30.0 * mu2 + 3.0) / 8.0;
    } else {
	return gsl_sf_legendre_Pl(ell, mu);
    }

}

/* spherical harmonics */
std::complex<double> calcSphericalHarmonics(int ell, int M, double * kvec) {
	
	if((ell ==0) && (M==0)) {
		return 1.0;
	}

	/*************************************/
	/* calc. mu1 and phi1 */
	/*************************************/
	double kmag = NORM(kvec);
        double mu = 0.0;
        double phi = 0.0;
        if( fabs(kmag) < 1.0e-10 ) {
                return 0.0;
        } else {
                mu = kvec[2] / kmag; 
        }
        double kmag_xy = sqrt(kvec[0] * kvec[0] + kvec[1] * kvec[1] );
        if( fabs(kmag_xy) < 1.0e-15 ) {
                phi = 0.0;
        } else {
                double cosphi = kvec[0] / kmag_xy;
                phi = acos(cosphi);
        }
        /* IMPORTANT! */
	if( kvec[1] < 0.0 ) { phi = - phi + 2.0 * M_PI; }
	/*************************************/

	std::complex<double> _I_(0.0, 1.0);
	std::complex<double> Ylm = 0.0;
	int S = abs(M);
	Ylm = pow(-1.0, (M-S)/2.0 ) * gsl_sf_legendre_sphPlm(ell, S, mu) * exp( _I_ * double(M) * phi );
	return Ylm;

}

/* spherical harmonics */
std::complex<double> calcSphericalHarmonics(int ell, int M, double mu, double phi) {
	
	if((ell ==0) && (M==0)) {
		return 1.0;
	}

	std::complex<double> _I_(0.0, 1.0);
	std::complex<double> Ylm = 0.0;
	int S = abs(M);
	Ylm = pow(-1.0, (M-S)/2.0 ) * gsl_sf_legendre_sphPlm(ell, S, mu) * exp( _I_ * double(M) * phi );
	return Ylm;

}


/* spherical harmonics */
std::complex<double> calcReducedSphericalHarmonics(int ell, int M, double mu, double phi) {
	
	if((ell ==0) && (M==0)) {
		return 1.0;
	}

	std::complex<double> _I_(0.0, 1.0);
	std::complex<double> Ylm = 0.0;
	int S = abs(M);
	Ylm = pow(-1.0, (M-S)/2.0 ) * gsl_sf_legendre_sphPlm(ell, S, mu) * exp( _I_ * double(M) * phi );
	Ylm = std::conj(Ylm);
	Ylm *= sqrt(4.0 * M_PI/(2.0 * double(ell) + 1.0));
	return Ylm;

}

double calcYYY(double mu1, double phi1, double mu2, double phi2, double mu, double phi, int ell1, int ell2, int ELL) {

	if(0) {
	} else if(ell1 == 0 && ell2 == 0 && ELL == 0) {
	    return 1.0;
	} else if((ell1 == ell2) && (ELL == 0)) {
	    double fac = 2.0 * double(ell1) + 1.0;
            double n1[3] = {sqrt(1.0 - mu1 * mu1) * cos(phi1), sqrt(1.0 - mu1 * mu1) * sin(phi1), mu1};
            double n2[3] = {sqrt(1.0 - mu2 * mu2) * cos(phi2), sqrt(1.0 - mu2 * mu2) * sin(phi2), mu2};
	    double mu12 = DOT(n1, n2);
	    return LegendrePolynomials(ell1, mu12) / fac;
	} else if((ell2 == ELL) && (ell1 == 0)) {
	    double fac = 2.0 * double(ELL) + 1.0;
            double  n[3] = {sqrt(1.0 - mu  * mu ) * cos(phi),  sqrt(1.0 - mu  * mu ) * sin(phi ), mu };
            double n2[3] = {sqrt(1.0 - mu2 * mu2) * cos(phi2), sqrt(1.0 - mu2 * mu2) * sin(phi2), mu2};
	    double mu02 = DOT(n, n2);
	    return LegendrePolynomials(ELL, mu02) / fac;
	} else if((ell1 == ELL) && (ell2 == 0)) {
	    double fac = 2.0 * double(ELL) + 1.0;
            double  n[3] = {sqrt(1.0 - mu  * mu ) * cos(phi),  sqrt(1.0 - mu  * mu ) * sin(phi ), mu };
            double n1[3] = {sqrt(1.0 - mu1 * mu1) * cos(phi1), sqrt(1.0 - mu1 * mu1) * sin(phi1), mu1};
	    double mu01 = DOT(n, n1);
	    return LegendrePolynomials(ELL, mu01) / fac;
	} else {
        	std::complex<double> S = 0.0;
        	for(int m1 = - ell1; m1 <= ell1; m1++) {
        	for(int m2 = - ell2; m2 <= ell2; m2++) {
        	for(int M  = -  ELL; M  <=  ELL;  M++) {
        
        		double w = wigner_3j(ell1, ell2, ELL, 0, 0, 0) * wigner_3j(ell1, ell2, ELL, m1, m2, M);
        	
        		if ( fabs(w) < 1.0e-10 ) {
        			continue;
        		}
        	
        		std::complex<double> Ylm1 = calcReducedSphericalHarmonics(ell1, m1, mu1, phi1);
        		std::complex<double> Ylm2 = calcReducedSphericalHarmonics(ell2, m2, mu2, phi2);
        		std::complex<double> Ylm  = calcReducedSphericalHarmonics( ELL,  M,  mu,  phi);
        
        		S += (w) * Ylm1 * Ylm2 * Ylm;
        	}}}
		return S.real();
	}
}

double calcYYY(double * kvec1, double * kvec2, double * los, int ell1, int ell2, int ELL) {


	/*************************************/
	/* calc. mu1 and phi1 */
	/*************************************/
	double kmag1 = NORM(kvec1);
    double mu1 = 0.0;
    double phi1 = 0.0;

    if( fabs(kmag1) < 1.0e-10 ) {
        return 0.0;
    } else {
            mu1 = kvec1[2] / kmag1; 
    }
    double kmag1_xy = sqrt(kvec1[0] * kvec1[0] + kvec1[1] * kvec1[1] );
    if( fabs(kmag1_xy) < 1.0e-15 ) {
            phi1 = 0.0;
    } else {
            double cosphi = kvec1[0] / kmag1_xy;
            phi1 = acos(cosphi);
    }
    /* IMPORTANT! */
	if( kvec1[1] < 0.0 ) { phi1 = - phi1 + 2.0 * M_PI; }
	/*************************************/

	/*************************************/
	/* calc. mu2 and phi2 */
	/*************************************/
	double kmag2 = NORM(kvec2);
        double mu2 = 0.0;
        double phi2 = 0.0;
        if( fabs(kmag2) < 1.0e-10 ) {
                return 0.0;
        } else {
                mu2 = kvec2[2] / kmag2; 
        }
        double kmag2_xy = sqrt(kvec2[0] * kvec2[0] + kvec2[1] * kvec2[1] );
        if( fabs(kmag2_xy) < 1.0e-15 ) {
                phi2 = 0.0;
        } else {
                double cosphi = kvec2[0] / kmag2_xy;
                phi2 = acos(cosphi);
        }
        /* IMPORTANT! */
	if( kvec2[1] < 0.0 ) { phi2 = - phi2 + 2.0 * M_PI; }
	/*************************************/

	/*************************************/
	/* calc. mu and phi */
	/*************************************/
	double lmag = NORM(los);
        double mu = 0.0;
        double phi = 0.0;
        if( fabs(lmag) < 1.0e-10 ) {
                return 0.0;
        } else {
                mu = los[2] / lmag; 
        }
        double lmag_xy = sqrt(los[0] * los[0] + los[1] * los[1] );
        if( fabs(lmag_xy) < 1.0e-15 ) {
                phi = 0.0;
        } else {
                double cosphi = los[0] / lmag_xy;
                phi = acos(cosphi);
        }
        /* IMPORTANT! */
	if( los[1] < 0.0 ) { phi = - phi + 2.0 * M_PI; }
	/*************************************/

	std::complex<double> S = 0.0;
	for(int m1 = - ell1; m1 <= ell1; m1++) {
	for(int m2 = - ell2; m2 <= ell2; m2++) {
	for(int M  = -  ELL; M  <=  ELL;  M++) {

		double w = wigner_3j(ell1, ell2, ELL, 0, 0, 0) * wigner_3j(ell1, ell2, ELL, m1, m2, M);
	
		if ( fabs(w) < 1.0e-10 ) {
			continue;
		}
	
		std::complex<double> Ylm1 = calcReducedSphericalHarmonics(ell1, m1, mu1, phi1);
		std::complex<double> Ylm2 = calcReducedSphericalHarmonics(ell2, m2, mu2, phi2);
		std::complex<double> Ylm  = calcReducedSphericalHarmonics( ELL,  M,  mu,  phi);

		S += (w) * Ylm1 * Ylm2 * Ylm;
    

	}}}

	return S.real();
}
 

int calc_mu_phi_from_vec(double * kvec_in, double & mu_out, double & phi_out) {

	double kmag = NORM(kvec_in);
        double mu = 0.0;
        double phi = 0.0;
        if( fabs(kmag) < 1.0e-10 ) {
                return 0.0;
        } else {
                mu = kvec_in[2] / kmag; 
        }

        double kmag_xy = sqrt(kvec_in[0] * kvec_in[0] + kvec_in[1] * kvec_in[1] );
        if( fabs(kmag_xy) < 1.0e-15 ) {
                phi = 0.0;
        } else {
                double cosphi = kvec_in[0] / kmag_xy;
                phi = acos(cosphi);
        }
        /* IMPORTANT! */
	if( kvec_in[1] < 0.0 ) { phi = - phi + 2.0 * M_PI; }

	mu_out = mu;
	phi_out = phi;

	return 0;
}

double calcYYY_YYY(double * kvec1, double * kvec2, double * kvec1_dash, double * kvec2_dash, double * los, int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash) {

	/*************************************/
	/* calc. mu1 and phi1 */
	/*************************************/
	double mu1 = 0.0;
	double phi1 = 0.0;
	calc_mu_phi_from_vec(kvec1, mu1, phi1);
	/*************************************/

	/*************************************/
	/* calc. mu2 and phi2 */
	/*************************************/
	double mu2 = 0.0;
	double phi2 = 0.0;
	calc_mu_phi_from_vec(kvec2, mu2, phi2);
	/*************************************/

	/*************************************/
	/* calc. mu1_dash and phi1_dash */
	/*************************************/
	double mu1_dash = 0.0;
	double phi1_dash = 0.0;
	calc_mu_phi_from_vec(kvec1_dash, mu1_dash, phi1_dash);
	/*************************************/

	/*************************************/
	/* calc. mu2_dash and phi2_dash */
	/*************************************/
	double mu2_dash = 0.0;
	double phi2_dash = 0.0;
	calc_mu_phi_from_vec(kvec2_dash, mu2_dash, phi2_dash);
	/*************************************/

	/*************************************/
	/* calc. mu and phi */
	/*************************************/
	double mu = 0.0;
	double phi = 0.0;
	calc_mu_phi_from_vec(los, mu, phi);
	/*************************************/

	std::complex<double> S = 0.0;
	for(int m1 = - ell1; m1 <= ell1; m1++) {
	for(int m2 = - ell2; m2 <= ell2; m2++) {
	for(int M  = -  ELL; M  <=  ELL;  M++) {

		double w000 = wigner_3j(ell1, ell2, ELL, 0, 0, 0);
		if ( fabs(w000) < 1.0e-10 ) {
			continue;
		}
		
        double w = wigner_3j(ell1, ell2, ELL, m1, m2, M) * w000;
	
		std::complex<double> Ylm1 = calcReducedSphericalHarmonics(ell1, m1, mu1, phi1);
		std::complex<double> Ylm2 = calcReducedSphericalHarmonics(ell2, m2, mu2, phi2);
		std::complex<double> Ylm  = calcReducedSphericalHarmonics( ELL,  M,  mu,  phi);

		S += (w) * Ylm1 * Ylm2 * Ylm;
	
	}}}

	std::complex<double> S_dash = 0.0;
	for(int m1 = - ell1_dash; m1 <= ell1_dash; m1++) {
	for(int m2 = - ell2_dash; m2 <= ell2_dash; m2++) {
	for(int M  = -  ELL_dash; M  <=  ELL_dash;  M++) {

	    double w000 =  wigner_3j(ell1_dash, ell2_dash, ELL_dash, 0, 0, 0);
		if ( fabs(w000) < 1.0e-10 ) {
			continue;
		}

		double w = wigner_3j(ell1_dash, ell2_dash, ELL_dash, m1, m2, M) / w000;
	
		std::complex<double> Ylm1 = calcReducedSphericalHarmonics(ell1_dash, m1, mu1_dash, phi1_dash);
		std::complex<double> Ylm2 = calcReducedSphericalHarmonics(ell2_dash, m2, mu2_dash, phi2_dash);
		std::complex<double> Ylm  = calcReducedSphericalHarmonics( ELL_dash,  M,  mu,  phi);

		S_dash += (w) * std::conj(Ylm1 * Ylm2 * Ylm);
	}}}

	return (S * S_dash).real();

}
 
double ClebschGordan(int ell1, int ell2, int ell3, int m1, int m2, int m3) {
    return pow(-1.0, ell1-ell2+m3) * sqrt(2.0*double(ell3)+1.0) 
          * wigner_3j(ell1, ell2, ell3, m1, m2, -m3);
}

std::complex<double> calcBipoSH(double mu1, double phi1, double mu2, double phi2, int ell1, int ell2, int _L_, int _M_) {

    std::complex<double> S = 0.0;
    for(int m1 = - ell1; m1 <= ell1; m1++) {
    for(int m2 = - ell2; m2 <= ell2; m2++) {
    
    	double w = wigner_3j(ell1, ell2, _L_, 0, 0, 0) * wigner_3j(ell1, ell2, _L_, m1, m2, - _M_);
    	if ( fabs(w) < 1.0e-10 ) {
    		continue;
    	}
    
    	std::complex<double> Ylm1 = calcReducedSphericalHarmonics(ell1, m1, mu1, phi1);
    	std::complex<double> Ylm2 = calcReducedSphericalHarmonics(ell2, m2, mu2, phi2);
    
    	S += (2.0 * double(_L_) + 1.0) * pow(-1.0, _M_) * (w) * Ylm1 * Ylm2;

    }}
    
    return S;

}

std::complex<double> calcTripoSH(double mu1, double phi1, double mu2, double phi2, double mu,
                	         double phi, int ell1, int ell2, int ell12, int ELL, int _J_, int _MJ_) {

	std::complex<double> S = 0.0;
	for(int m1 = - ell1; m1 <= ell1; m1++) {
	for(int m2 = - ell2; m2 <= ell2; m2++) {
	for(int m12 = - ell12; m12 <= ell12; m12++) {
	for(int M  = -  ELL; M  <=  ELL;  M++) {

		double C1 = ClebschGordan(ell1, ell2, ell12, m1, m2, m12);
		double C2 = ClebschGordan(ell12, ELL, _J_, m12, M, _MJ_);
		double C = C1 * C2;
		if ( fabs(C) < 1.0e-10 ) {
			continue;
		}
	
		std::complex<double> Ylm1 = calcReducedSphericalHarmonics(ell1, m1, mu1, phi1);
		std::complex<double> Ylm2 = calcReducedSphericalHarmonics(ell2, m2, mu2, phi2);
		std::complex<double> Ylm  = calcReducedSphericalHarmonics( ELL,  M,  mu,  phi);

		S += (C) * Ylm1 * Ylm2 * Ylm;
	}}}}

	return S;
}

std::complex<double> calcTripoSH(double * kvec1, double * kvec2, double *los, int ell1, int ell2, int ell12, int ELL, int _J_, int _MJ_) {

	/*************************************/
	/* calc. mu1 and phi1 */
	/*************************************/
	double kmag1 = NORM(kvec1);
        double mu1 = 0.0;
        double phi1 = 0.0;
        if( fabs(kmag1) < 1.0e-10 ) {
                return 0.0;
        } else {
                mu1 = kvec1[2] / kmag1; 
        }
        double kmag1_xy = sqrt(kvec1[0] * kvec1[0] + kvec1[1] * kvec1[1] );
        if( fabs(kmag1_xy) < 1.0e-15 ) {
                phi1 = 0.0;
        } else {
                double cosphi = kvec1[0] / kmag1_xy;
                phi1 = acos(cosphi);
        }
        /* IMPORTANT! */
	if( kvec1[1] < 0.0 ) { phi1 = - phi1 + 2.0 * M_PI; }
	/*************************************/

	/*************************************/
	/* calc. mu2 and phi2 */
	/*************************************/
	double kmag2 = NORM(kvec2);
        double mu2 = 0.0;
        double phi2 = 0.0;
        if( fabs(kmag2) < 1.0e-10 ) {
                return 0.0;
        } else {
                mu2 = kvec2[2] / kmag2; 
        }
        double kmag2_xy = sqrt(kvec2[0] * kvec2[0] + kvec2[1] * kvec2[1] );
        if( fabs(kmag2_xy) < 1.0e-15 ) {
                phi2 = 0.0;
        } else {
                double cosphi = kvec2[0] / kmag2_xy;
                phi2 = acos(cosphi);
        }
        /* IMPORTANT! */
	if( kvec2[1] < 0.0 ) { phi2 = - phi2 + 2.0 * M_PI; }
	/*************************************/

	/*************************************/
	/* calc. mu and phi */
	/*************************************/
	double lmag = NORM(los);
        double mu = 0.0;
        double phi = 0.0;
        if( fabs(lmag) < 1.0e-10 ) {
                return 0.0;
        } else {
                mu = los[2] / lmag; 
        }
        double lmag_xy = sqrt(los[0] * los[0] + los[1] * los[1] );
        if( fabs(lmag_xy) < 1.0e-15 ) {
                phi = 0.0;
        } else {
                double cosphi = los[0] / lmag_xy;
                phi = acos(cosphi);
        }
        /* IMPORTANT! */
	if( los[1] < 0.0 ) { phi = - phi + 2.0 * M_PI; }
	/*************************************/




	std::complex<double> S = 0.0;
	for(int m1 = - ell1; m1 <= ell1; m1++) {
	for(int m2 = - ell2; m2 <= ell2; m2++) {
	for(int m12 = - ell12; m12 <= ell12; m12++) {
	for(int M  = -  ELL; M  <=  ELL;  M++) {

		double C1 = ClebschGordan(ell1, ell2, ell12, m1, m2, m12);
		double C2 = ClebschGordan(ell12, ELL, _J_, m12, M, _MJ_);
		double C = C1 * C2;
		if ( fabs(C) < 1.0e-10 ) {
			continue;
		}
	
		std::complex<double> Ylm1 = calcReducedSphericalHarmonics(ell1, m1, mu1, phi1);
		std::complex<double> Ylm2 = calcReducedSphericalHarmonics(ell2, m2, mu2, phi2);
		std::complex<double> Ylm  = calcReducedSphericalHarmonics( ELL,  M,  mu,  phi);

		S += (C) * Ylm1 * Ylm2 * Ylm;
	}}}}

	return S;
}

double calcDeltaFunction(double k1, double k2, double DeltaK) {

	if( (k1 < 1.0e-10) || (k2 < 1.0e-10) ) {
		return 0.0;
	}

	double x = fabs(k1-k2);
	double D = pow(2.0*M_PI,3) / (4.0 * M_PI * k1 * k2 * DeltaK);

	if( x < DeltaK / 2.0) {
	    return D;
	} else {
	    return 0.0;
	}

}


#endif
