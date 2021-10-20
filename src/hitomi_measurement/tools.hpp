#ifndef __tools__
#define __tools__

#ifndef __parameter__
#include "parameter.hpp"
#endif

class ToolClass {
private:
public:

	static std::complex<double> calcReducedSphericalHarmonics(int ell, int M, double * pos) {

        /*************************************************************************************************/
        /* This function calculates the complex conjugate of the normalized spherical harmonic function. */
        /*************************************************************************************************/

		if((ell ==0) && (M==0)) {
			return 1.0;
		}

		const int dim = 3;
		double xmag_xyz2 = 0.0;
		for(int axes = 0; axes < dim; axes++) {
			xmag_xyz2 += pos[axes] * pos[axes];
		}
		double xmag_xyz = sqrt(xmag_xyz2);
		double xmag_xy = sqrt(pos[0] * pos[0] + pos[1] * pos[1] );

		/* Calculate mu. */
		double mu = 0.0;
		if( fabs(xmag_xyz) < 1.0e-15 ) {
			return 0.0;
		} else {
			mu = pos[2] / xmag_xyz;
		}
		
		/* Calculate phi. */
		double phi = 0.0;
		if( fabs(xmag_xy) < 1.0e-15 ) {
			phi = 0.0;
		} else {
			double cosphi = pos[0] / xmag_xy;
			phi = acos(cosphi);
		}

		if(pos[1] < 0.0) { phi = - phi + 2.0 * M_PI; }

		/* Define an imaginary number, _I_*/
		std::complex<double> _I_(0.0, 1.0);

        /* gsl_sf_legendre_sphPlm() computes the normalized Legendre function, 
         * sqrt((2l+1)/4pi) * sqrt(l-m/l+m) P_lm */
        /* Since 'gsl_sf_legendre_sphPlm()' can be computed only for M>=0, 
         * if M<0, the absolute value of M is computed and substituted into 'gsl_sf_legendre_sphPlm()' */
		/* Ylm = (-1)**((m-|s|)/2) * sqrt(l-|m|/(l+|m|)) * P_l^|m| * exp(i m phi) */
		std::complex<double> Ylm = 0.0;
		int S = abs(M);
		Ylm = pow(-1.0, (M-S)/2.0 ) * gsl_sf_legendre_sphPlm(ell, S, mu) * exp( _I_ * double(M) * phi );
		
        /* Compute the complex conjecture of Ylm */
        Ylm = std::conj(Ylm);
        /* Multiply sqrt(4pi/(2ell+1)) by Ylm to normalize */
		Ylm *= sqrt(4.0 * M_PI/(2.0 * double(ell) + 1.0));
		
        return Ylm;
	
    }

	static int storeReducedSphericalHarmonicsInFourierSpace(int _ell_, int _m_, ParameterClass & param, std::complex<double> * Ylm_out ) {
	
		if(Ylm_out == NULL) {
			return -1;
		}
		double kvec[3];
		double dk[3];
		dk[0] = 2.0 * M_PI / param.boxsize[0];
		dk[1] = 2.0 * M_PI / param.boxsize[1];
		dk[2] = 2.0 * M_PI / param.boxsize[2];
		for(int i = 0; i < param.n_mesh[0]; i++) {
		for(int j = 0; j < param.n_mesh[1]; j++) {
		for(int k = 0; k < param.n_mesh[2]; k++) {
			long long coord = ( i * param.n_mesh[1] + j ) * param.n_mesh[2] + k;
			kvec[0] = (i < param.n_mesh[0]/2) ? (double) i * dk[0] : (double) (i - param.n_mesh[0]) * dk[0];
			kvec[1] = (j < param.n_mesh[1]/2) ? (double) j * dk[1] : (double) (j - param.n_mesh[1]) * dk[1];
			kvec[2] = (k < param.n_mesh[2]/2) ? (double) k * dk[2] : (double) (k - param.n_mesh[2]) * dk[2];
			Ylm_out[coord] = calcReducedSphericalHarmonics(_ell_, _m_, kvec);
		}}}

		return 0;
	}

	static int storeReducedSphericalHarmonicsInConfigurationSpace(int _ell_, int _m_, ParameterClass & param, std::complex<double> * Ylm_out ) {
	
		if(Ylm_out == NULL) {
			return -1;
		}

		double rvec[3];
		double dr[3];
		dr[0] = param.boxsize[0] / double(param.n_mesh[0]);
		dr[1] = param.boxsize[1] / double(param.n_mesh[1]);
		dr[2] = param.boxsize[2] / double(param.n_mesh[2]);
		for(int i = 0; i < param.n_mesh[0]; i++) {
		for(int j = 0; j < param.n_mesh[1]; j++) {
		for(int k = 0; k < param.n_mesh[2]; k++) {
			long long coord = ( i * param.n_mesh[1] + j ) * param.n_mesh[2] + k;
			rvec[0] = (i < param.n_mesh[0]/2) ? (double) i * dr[0] : (double) (i - param.n_mesh[0]) * dr[0];
			rvec[1] = (j < param.n_mesh[1]/2) ? (double) j * dr[1] : (double) (j - param.n_mesh[1]) * dr[1];
			rvec[2] = (k < param.n_mesh[2]/2) ? (double) k * dr[2] : (double) (k - param.n_mesh[2]) * dr[2];
			Ylm_out[coord] = calcReducedSphericalHarmonics(_ell_, _m_, rvec);
		}}}

		return 0;
	}


	static int setKbin(ParameterClass & param, double * kbin_out) {
		double dk = (param.kmax - param.kmin) / double(param.n_kbin-1);
		for(int i = 0; i < param.n_kbin; i++) {
		        kbin_out[i] = param.kmin + dk * double(i);
		}
		return 0;
	}

	static int setRbin(ParameterClass & param, double * rbin_out) {
		double dr = (param.rmax - param.rmin) / double(param.n_rbin-1);
		for(int i = 0; i < param.n_rbin; i++) {
		        rbin_out[i] = param.rmin + dr * double(i);
		}
		return 0;
	}
};

#endif

