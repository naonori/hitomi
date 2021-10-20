#ifndef __wigner__
#define __wigner__
#include <gsl/gsl_sf_coupling.h>

double w3j[10][10][10][19][19][19];

int setWigner3j(void) {
	const int ell_wigner = 9;
	const int m_wigner = 2 * ell_wigner + 1;

	for(int ell1 = 0; ell1 <= ell_wigner; ell1++) {
	for(int ell2 = 0; ell2 <= ell_wigner; ell2++) {
	for(int ELL = 0; ELL <= ell_wigner; ELL++) {
	
		for(int m1 = 0; m1 < m_wigner; m1++) {
		for(int m2 = 0; m2 < m_wigner; m2++) {
		for(int M =  0; M < m_wigner; M++) {
			w3j[ell1][ell2][ELL][m1][m2][M] = 0.0;
		}}}
	}}}


	for(int ell1 = 0; ell1 <= ell_wigner; ell1++) {
	for(int ell2 = 0; ell2 <= ell_wigner; ell2++) {
	for(int ELL = 0; ELL <= ell_wigner; ELL++) {
	
		for(int m1 = -ell1; m1 <= ell1; m1++) {
		for(int m2 = -ell2; m2 <= ell2; m2++) {
		for(int M = -ELL; M <= ELL; M++) {
			double w =  gsl_sf_coupling_3j(2*ell1,2*ell2,2*ELL,2*m1,2*m2,2*M);
			if(fabs(w) < 1.0e-10) {
				continue;
			}
			w3j[ell1][ell2][ELL][m1+ell1][m2+ell2][M+ELL] = w;
		
		}}}

	}}}



	return 0;
}

#endif

