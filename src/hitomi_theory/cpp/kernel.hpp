#ifndef __kernel__
#define __kernel__

#ifndef __common__
#include "common.hpp"
#endif

#ifndef __pk_lin__
#include "pk_lin.hpp"
#endif

#ifndef __sph__
#include "sph.hpp"
#endif

double MU(double * kvec1, double * kvec2) {
 	double k1 = NORM(kvec1);
 	double k2 = NORM(kvec2);
 	double result = 0.0;
 	if( (k1 > pk_kmin) && (k2 > pk_kmin)) {
 		double m = DOT(kvec1, kvec2) / k1 / k2;
 		if(m > 1.0) {
 			m = 1.0;
 		} else if (m < -1.0) {
 			m = - 1.0;
 		}
 		result += m;
 	}
 	return result;
}

double W1(double kmag, double R) {

	double kr = R * kmag;
	double kr2 = kr * kr;
	double kr3 = kr2 * kr;
  	double w = 3.0 * (sin(kr) / kr3 - cos(kr) / kr2);
	return w;
}


double W2(double kmag, double R) {

	double kr = R * kmag;
	double kr2 = kr * kr;
	double kr3 = kr2 * kr;
  	double w = 3.0 * (sin(kr) / kr3 - cos(kr) / kr2);
	double w2 = w * w;

	return w2;
}

int calcTrueWavevector(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double * kvec_out) {

	double kn = DOT(kvec_in, los);

	double kvec_para[3] = { kn * los[0] / alpha_parallel, kn * los[1] / alpha_parallel, kn * los[2] / alpha_parallel };
	double kvec_perp[3] = { (kvec_in[0] - kn * los[0]) / alpha_perp, (kvec_in[1] - kn * los[1]) / alpha_perp, (kvec_in[2] - kn * los[2]) / alpha_perp };

	double kvec[3] = PLUS(kvec_para, kvec_perp);
	
	kvec_out[0] = kvec[0];
	kvec_out[1] = kvec[1];
	kvec_out[2] = kvec[2];

	return 0;

}


/***********************************************/
/* Kernel functions for dark matter: Fn and Gn */
/***********************************************/
double alpha(double * kvec1, double * kvec2) {
	double k1 = NORM(kvec1);
	double kvec12[3] = PLUS(kvec1, kvec2);
	double result = 0.0;
	if( (k1 > pk_kmin) ) {
		result += DOT(kvec12, kvec1) / k1 / k1;
	}
	
	return result;
}

double beta(double * kvec1, double * kvec2) {
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double kvec12[3] = PLUS(kvec1, kvec2);
	double result = 0.0;
	if( (k1 > pk_kmin) && (k2 > pk_kmin)) {
		result += DOT(kvec12, kvec12) * DOT(kvec1, kvec2) / 2.0 / pow(k1,2) / pow(k2,2);
	}
	return result;
}

/***********************************************/

double F2(double * kvec1, double * kvec2) {
	return (5.0/14.0) * alpha(kvec1, kvec2) + (5.0/14.0) * alpha(kvec2, kvec1) + (2.0/7.0) * beta(kvec1, kvec2);
} 

double G2(double * kvec1, double * kvec2) {
	return (3.0/14.0) * alpha(kvec1, kvec2) + (3.0/14.0) * alpha(kvec2, kvec1) + (4.0/7.0) * beta(kvec1, kvec2);
}

double F2_Growth(double * kvec1, double * kvec2) {
	return (17.0/21.0);
} 

double F2_Shift(double * kvec1, double * kvec2) {
    double mu = MU(kvec1, kvec2);
    double k1 = NORM(kvec1);
    double k2 = NORM(kvec2);
	return (1.0/2.0) * mu * (k1/k2 + k2/k1);
} 

double F2_Tidal(double * kvec1, double * kvec2) {
    double mu = MU(kvec1, kvec2);
	return (2.0/7.0) * (mu * mu - (1.0/3.0));
} 


/***********************************************/

double F3_temp(double * kvec1, double * kvec2, double * kvec3) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec23[3] = PLUS(kvec2, kvec3);
	double F = (7.0/18.0) * alpha(kvec1, kvec23) * F2(kvec2, kvec3) + (1.0/9.0) * beta(kvec1, kvec23) * G2(kvec2, kvec3)
             + (7.0/18.0) * alpha(kvec12, kvec3) * G2(kvec1, kvec2) + (1.0/9.0) * beta(kvec12, kvec3) * G2(kvec1, kvec2); 
	return F;
} 

double F3(double * kvec1, double * kvec2, double * kvec3) {
	return ( F3_temp(kvec1, kvec2, kvec3) + F3_temp(kvec2, kvec3, kvec1) + F3_temp(kvec3, kvec1, kvec2) ) / 3.0;
} 

double G3_temp(double * kvec1, double * kvec2, double * kvec3) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec23[3] = PLUS(kvec2, kvec3);
	double G = (1.0/6.0) * alpha(kvec1, kvec23) * F2(kvec2, kvec3) + (1.0/3.0) * beta(kvec1, kvec23) * G2(kvec2, kvec3)
                 + (1.0/6.0) * alpha(kvec12, kvec3) * G2(kvec1, kvec2) + (1.0/3.0) * beta(kvec12, kvec3) * G2(kvec1, kvec2); 
	return G;
} 

double G3(double * kvec1, double * kvec2, double * kvec3) {
	return ( G3_temp(kvec1, kvec2, kvec3) + G3_temp(kvec2, kvec3, kvec1) + G3_temp(kvec3, kvec1, kvec2) ) / 3.0;
}


double F4_temp_1(double * kvec1, double * kvec2, double * kvec3, double * kvec4) {

	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double kvec234[3] = PLUS3(kvec2, kvec3, kvec4);

	double F = (9.0) * alpha(kvec1, kvec234) * F3(kvec2, kvec3, kvec4) 
		 + (2.0) *  beta(kvec1, kvec234) * G3(kvec2, kvec3, kvec4)
		 + (9.0) * G3(kvec1, kvec2, kvec3) * alpha(kvec123, kvec4)
		 + (2.0) * G3(kvec1, kvec2, kvec3) * beta(kvec123, kvec4);

	return F / 33.0;
} 

double F4_temp_2(double * kvec1, double * kvec2, double * kvec3, double * kvec4) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec34[3] = PLUS(kvec3, kvec4);

	double F = (9.0) * G2(kvec1, kvec2) * alpha(kvec12, kvec34) * F2(kvec3, kvec4) 
		 + (2.0) * G2(kvec1, kvec2) * beta(kvec12, kvec34) * G2(kvec3, kvec4);

	return F / 33.0;
} 

double F4(double * kvec1, double * kvec2, double * kvec3, double * kvec4) {

	double Ftemp1 = F4_temp_1(kvec1, kvec2, kvec3, kvec4) 
	              + F4_temp_1(kvec2, kvec3, kvec4, kvec1) 
	              + F4_temp_1(kvec3, kvec4, kvec1, kvec2) 
	              + F4_temp_1(kvec4, kvec1, kvec2, kvec3);

	double Ftemp2 = F4_temp_2(kvec1, kvec2, kvec3, kvec4)
		      + F4_temp_2(kvec1, kvec3, kvec2, kvec4)
		      + F4_temp_2(kvec1, kvec4, kvec2, kvec3)
		      + F4_temp_2(kvec2, kvec3, kvec1, kvec4)
		      + F4_temp_2(kvec2, kvec4, kvec1, kvec3)
		      + F4_temp_2(kvec3, kvec4, kvec1, kvec2);

	return Ftemp1 / 4.0 + Ftemp2 / 6.0;
}


double G4_temp_1(double * kvec1, double * kvec2, double * kvec3, double * kvec4) {
	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double kvec234[3] = PLUS3(kvec2, kvec3, kvec4);

	double G = (3.0) * alpha(kvec1, kvec234) * F3(kvec2, kvec3, kvec4) 
                 + (8.0) *  beta(kvec1, kvec234) * G3(kvec2, kvec3, kvec4)
		 + (3.0) * G3(kvec1, kvec2, kvec3) * alpha(kvec123, kvec4)
		 + (8.0) * G3(kvec1, kvec2, kvec3) * beta(kvec123, kvec4);

	return G / 33.0;
} 

double G4_temp_2(double * kvec1, double * kvec2, double * kvec3, double * kvec4) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec34[3] = PLUS(kvec3, kvec4);

	double G = (3.0) * G2(kvec1, kvec2) * alpha(kvec12, kvec34) * F2(kvec3, kvec4) 
		 + (8.0) * G2(kvec1, kvec2) * beta(kvec12, kvec34) * G2(kvec3, kvec4);

	return G / 33.0;
} 


double G4(double * kvec1, double * kvec2, double * kvec3, double * kvec4) {

	double Gtemp1 = G4_temp_1(kvec1, kvec2, kvec3, kvec4) 
		      + G4_temp_1(kvec2, kvec3, kvec4, kvec1) 
		      + G4_temp_1(kvec3, kvec4, kvec1, kvec2) 
		      + G4_temp_1(kvec4, kvec1, kvec2, kvec3);
	double Gtemp2 = G4_temp_2(kvec1, kvec2, kvec3, kvec4)
		      + G4_temp_2(kvec1, kvec3, kvec2, kvec4)
		      + G4_temp_2(kvec1, kvec4, kvec2, kvec3)
		      + G4_temp_2(kvec2, kvec3, kvec1, kvec4)
		      + G4_temp_2(kvec2, kvec4, kvec1, kvec3)
		      + G4_temp_2(kvec3, kvec4, kvec1, kvec2);

	return Gtemp1 / 4.0 + Gtemp2 / 6.0;
}

double F5_temp_1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {
	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double kvec2345[3] = PLUS4(kvec2, kvec3, kvec4, kvec5);

	double F = (3.0) * alpha(kvec1, kvec2345) * F4(kvec2, kvec3, kvec4, kvec5)
		 + (10.0) * beta(kvec1, kvec2345) * G4(kvec2, kvec3, kvec4, kvec5)
		 + (3.0) * alpha(kvec1234, kvec5) * G4(kvec1, kvec2, kvec3, kvec4)
		 + (10.0) * beta(kvec1234, kvec5) * G4(kvec1, kvec2, kvec3, kvec4);

	return F / 52.0;
} 

double F5_temp_2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec345[3] = PLUS3(kvec3, kvec4, kvec5);

	double F = (3.0) * G2(kvec1, kvec2) * alpha(kvec12, kvec345) * F3(kvec3, kvec4, kvec5)
		 + (10.0) * G2(kvec1, kvec2) * beta(kvec12, kvec345) * G3(kvec3, kvec4, kvec5);
	return F / 52.0;
} 

double F5_temp_3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {
	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double kvec45[3] = PLUS(kvec4, kvec5);

	double F = (3.0) * G3(kvec1, kvec2, kvec3) * alpha(kvec123, kvec45) * F2(kvec4, kvec5)
		 + (10.0) * G3(kvec1, kvec2, kvec3) * beta(kvec123, kvec45) * G2(kvec4, kvec5);
	return F / 52.0;
} 

double F5(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {

	double Ftemp1 = F5_temp_1(kvec1, kvec2, kvec3, kvec4, kvec5)
		      + F5_temp_1(kvec2, kvec3, kvec4, kvec5, kvec1)
		      + F5_temp_1(kvec3, kvec4, kvec5, kvec1, kvec2)
		      + F5_temp_1(kvec4, kvec5, kvec1, kvec2, kvec3)
		      + F5_temp_1(kvec5, kvec1, kvec2, kvec3, kvec4);

	double Ftemp2 = F5_temp_2(kvec1, kvec2, kvec3, kvec4, kvec5)
	              + F5_temp_2(kvec1, kvec3, kvec2, kvec4, kvec5)
	              + F5_temp_2(kvec1, kvec4, kvec2, kvec3, kvec5)
	              + F5_temp_2(kvec1, kvec5, kvec2, kvec3, kvec4)
	              + F5_temp_2(kvec2, kvec3, kvec1, kvec4, kvec5)
	              + F5_temp_2(kvec2, kvec4, kvec1, kvec3, kvec5)
	              + F5_temp_2(kvec2, kvec5, kvec1, kvec3, kvec4)
	              + F5_temp_2(kvec3, kvec4, kvec1, kvec2, kvec5)
	              + F5_temp_2(kvec3, kvec5, kvec1, kvec2, kvec4)
	              + F5_temp_2(kvec4, kvec5, kvec1, kvec2, kvec3);


	double Ftemp3 = F5_temp_3(kvec3, kvec4, kvec5, kvec1, kvec2)
	              + F5_temp_3(kvec2, kvec4, kvec5, kvec1, kvec3)
	              + F5_temp_3(kvec2, kvec3, kvec5, kvec1, kvec4)
	              + F5_temp_3(kvec2, kvec3, kvec4, kvec1, kvec5)
	              + F5_temp_3(kvec1, kvec4, kvec5, kvec2, kvec3)
	              + F5_temp_3(kvec1, kvec3, kvec5, kvec2, kvec4)
	              + F5_temp_3(kvec1, kvec3, kvec4, kvec2, kvec5)
	              + F5_temp_3(kvec1, kvec2, kvec5, kvec3, kvec4)
	              + F5_temp_3(kvec1, kvec2, kvec4, kvec3, kvec5)
	              + F5_temp_3(kvec1, kvec2, kvec3, kvec4, kvec5);

	return Ftemp1 / 5.0 + Ftemp2 / 10.0 + Ftemp3 / 10.0;
}

double G5_temp_1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {
	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double kvec2345[3] = PLUS4(kvec2, kvec3, kvec4, kvec5);

	double G = (3.0) * alpha(kvec1, kvec2345) * F4(kvec2, kvec3, kvec4, kvec5)
		 + (10.0) * beta(kvec1, kvec2345) * G4(kvec2, kvec3, kvec4, kvec5)
		 + (3.0) * alpha(kvec1234, kvec5) * G4(kvec1, kvec2, kvec3, kvec4)
		 + (10.0) * beta(kvec1234, kvec5) * G4(kvec1, kvec2, kvec3, kvec4);

	return G / 52.0;
} 

double G5_temp_2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec345[3] = PLUS3(kvec3, kvec4, kvec5);

	double G = (3.0) * G2(kvec1, kvec2) * alpha(kvec12, kvec345) * F3(kvec3, kvec4, kvec5)
		 + (10.0) * G2(kvec1, kvec2) * beta(kvec12, kvec345) * G3(kvec3, kvec4, kvec5);
	return G / 52.0;
} 

double G5_temp_3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {
	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double kvec45[3] = PLUS(kvec4, kvec5);

	double G = (3.0) * G3(kvec1, kvec2, kvec3) * alpha(kvec123, kvec45) * F2(kvec4, kvec5)
		 + (10.0) * G3(kvec1, kvec2, kvec3) * beta(kvec123, kvec45) * G2(kvec4, kvec5);
	return G / 52.0;
} 

double G5(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {

	double Gtemp1 = G5_temp_1(kvec1, kvec2, kvec3, kvec4, kvec5)
		      + G5_temp_1(kvec2, kvec3, kvec4, kvec5, kvec1)
		      + G5_temp_1(kvec3, kvec4, kvec5, kvec1, kvec2)
		      + G5_temp_1(kvec4, kvec5, kvec1, kvec2, kvec3)
		      + G5_temp_1(kvec5, kvec1, kvec2, kvec3, kvec4);

	double Gtemp2 = G5_temp_2(kvec1, kvec2, kvec3, kvec4, kvec5)
	              + G5_temp_2(kvec1, kvec3, kvec2, kvec4, kvec5)
	              + G5_temp_2(kvec1, kvec4, kvec2, kvec3, kvec5)
	              + G5_temp_2(kvec1, kvec5, kvec2, kvec3, kvec4)
	              + G5_temp_2(kvec2, kvec3, kvec1, kvec4, kvec5)
	              + G5_temp_2(kvec2, kvec4, kvec1, kvec3, kvec5)
	              + G5_temp_2(kvec2, kvec5, kvec1, kvec3, kvec4)
	              + G5_temp_2(kvec3, kvec4, kvec1, kvec2, kvec5)
	              + G5_temp_2(kvec3, kvec5, kvec1, kvec2, kvec4)
	              + G5_temp_2(kvec4, kvec5, kvec1, kvec2, kvec3);


	double Gtemp3 = G5_temp_3(kvec3, kvec4, kvec5, kvec1, kvec2)
	              + G5_temp_3(kvec2, kvec4, kvec5, kvec1, kvec3)
	              + G5_temp_3(kvec2, kvec3, kvec5, kvec1, kvec4)
	              + G5_temp_3(kvec2, kvec3, kvec4, kvec1, kvec5)
	              + G5_temp_3(kvec1, kvec4, kvec5, kvec2, kvec3)
	              + G5_temp_3(kvec1, kvec3, kvec5, kvec2, kvec4)
	              + G5_temp_3(kvec1, kvec3, kvec4, kvec2, kvec5)
	              + G5_temp_3(kvec1, kvec2, kvec5, kvec3, kvec4)
	              + G5_temp_3(kvec1, kvec2, kvec4, kvec3, kvec5)
	              + G5_temp_3(kvec1, kvec2, kvec3, kvec4, kvec5);

	return Gtemp1 / 5.0 + Gtemp2 / 10.0 + Gtemp3 / 10.0;
}

/***********************************************/

/*****************************/
/* kernel functions for bias */
/*****************************/

/* second order */
double D1D1(double * kvec1, double * kvec2) {
	return 1.0;
}

double K1K1(double * kvec1, double * kvec2) {
	double mu = MU(kvec1, kvec2);
	return mu * mu - 1.0/3.0;
}

double F2_Bias(double * kvec1, double * kvec2, double b1, double b2, double bK2) {
	return   b1 * F2(kvec1, kvec2) 
	      +  (b2/2.0) * D1D1(kvec1, kvec2) 
	      +  bK2 * K1K1(kvec1, kvec2);
}


/* third order */

double D1D2(double * kvec1, double * kvec2, double * kvec3) {
	return ( F2(kvec1, kvec2) + F2(kvec1,kvec3) + F2(kvec2,kvec3) ) / 3.0;
}

double K1K2(double * kvec1, double * kvec2, double * kvec3) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec23[3] = PLUS(kvec2, kvec3);

	return ( K1K1(kvec1, kvec23) * F2(kvec2, kvec3) 
	       + K1K1(kvec2, kvec13) * F2(kvec1, kvec3)
	       + K1K1(kvec3, kvec12) * F2(kvec1, kvec2) ) / 3.0; 
}

double D1D1D1(double * kvec1, double * kvec2, double * kvec3) {
	return 1.0;
}

double K1K1K1(double * kvec1, double * kvec2, double * kvec3) {
	double mu12 = MU(kvec1, kvec2);
	double mu13 = MU(kvec1, kvec3);
	double mu23 = MU(kvec2, kvec3);

	return mu12 * mu13 * mu23 - (mu12*mu12 + mu13*mu13 + mu23*mu23) / 3.0 + (2.0/9.0);

}

double D1K1K1(double * kvec1, double * kvec2, double * kvec3) {
	double mu12 = MU(kvec1, kvec2);
	double mu13 = MU(kvec1, kvec3);
	double mu23 = MU(kvec2, kvec3);

	return (mu12 * mu12 + mu13 * mu13 + mu23 * mu23) / 3.0  - (1.0/3.0);
}

double O3(double * kvec1, double * kvec2, double * kvec3) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec23[3] = PLUS(kvec2, kvec3);

	double mu12 = MU(kvec1, kvec2);
	double mu13 = MU(kvec1, kvec3);
	double mu23 = MU(kvec2, kvec3);

	double K23 = K1K1(kvec1, kvec23);
	double K13 = K1K1(kvec2, kvec13);
	double K12 = K1K1(kvec3, kvec12);

	double KDS2 = (4.0/7.0) * ( K23 * (1.0 - mu23 * mu23) + K13 * (1.0 - mu13 * mu13) + K12 * (1.0 - mu12 * mu12) ) / 3.0; 

	return KDS2;

}

double F3_Bias(double * kvec1, double * kvec2, double * kvec3, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO) {

	return b1 * F3(kvec1, kvec2, kvec3)
	    +  b2 * D1D2(kvec1, kvec2, kvec3)
	    +  2.0 * bK2 * K1K2(kvec1, kvec2, kvec3)
	    +  (b3 / 6.0) * D1D1D1(kvec1, kvec2, kvec3)
	    +  bK3 * K1K1K1(kvec1, kvec2, kvec3)
	    +  bDK * D1K1K1(kvec1, kvec2, kvec3)
	    +  bO * O3(kvec1, kvec2, kvec3);

}

/* forth order */

double F4_Bias(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double b1) {
	return b1 * F4(kvec1, kvec2, kvec3, kvec4);
}

/* fifth order */

double F5_Bias(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double b1) {
	return b1 * F5(kvec1, kvec2, kvec3, kvec4, kvec5);
}

/******************************/
/* kernel functions with RSDs */
/******************************/

double LV1(double * kvec1, double * los) {

	double k1 = NORM(kvec1);
	double result = 0.0;
	if( k1 > pk_kmin ) {
		double k1n = DOT(kvec1, los);
		result += k1n / pow(k1,2);
	}
	return result;
}

double LV2(double * kvec1, double * kvec2, double * los) {

	double kvec12[3] = PLUS(kvec1, kvec2);
	double k12 = NORM(kvec12);

	double result = 0.0;
	if( k12 > pk_kmin ) {
		double k12n = DOT(kvec12, los);
		result += k12n * G2(kvec1, kvec2) / pow(k12,2);
	}
	return result;
}

double LV3(double * kvec1, double * kvec2, double * kvec3, double * los) {

	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double k123 = NORM(kvec123);

	double result = 0.0;
	if( k123 > pk_kmin ) {
		double k123n = DOT(kvec123, los);
		result += k123n * G3(kvec1, kvec2, kvec3) / pow(k123,2);
	}
	return result;
}

double LV4(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234 = NORM(kvec1234);

	double result = 0.0;
	if( k1234 > pk_kmin ) {
		double k1234n = DOT(kvec1234, los);
		result += k1234n * G4(kvec1, kvec2, kvec3, kvec4) / pow(k1234,2);
	}
	return result;
}


/* first order */

double Z1_Bias(double * kvec1, double * los, double f, double b1) {
	double mu = MU(kvec1, los);
	return b1 + f * mu * mu;
}


/* second order */
double V2(double * kvec1, double * kvec2, double * los, double f) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double mu = MU(kvec12, los);
	
	return f * mu * mu * G2(kvec1, kvec2);
}

double V1V1(double * kvec1, double * kvec2, double * los, double f) {

	double kvec12[3] = PLUS(kvec1, kvec2);
	double kn  = DOT(kvec12, los); 
	double result = (f * f) * (kn * kn) * LV1(kvec1,los) * LV1(kvec2,los);
	return result;
}

double D1V1(double * kvec1, double * kvec2, double * los, double f) {
	
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kn  = DOT(kvec12, los); 
	
	double result = (f) * (kn) * LV1(kvec1,los) / 2.0
		          + (f) * (kn) * LV1(kvec2,los) / 2.0;
	return result;
}

double Z2_Bias(double * kvec1, double * kvec2, double * los, double f, double b1, double b2, double bK2) {
	return  V2(kvec1, kvec2, los, f)
	      + V1V1(kvec1, kvec2, los, f) / 2.0
	      + b1 * D1V1(kvec1, kvec2, los, f) 
	      + F2_Bias(kvec1, kvec2, b1, b2, bK2);
}

/********************************/
/*  decomposed kernel functions */
/********************************/

double D1() {
	return 1.0;
}

double V1(double * kvec1, double * los) {
	double mu = MU(kvec1, los);
	return mu * mu;
}
    
double V2(double * kvec1, double * kvec2, double * los) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double mu = MU(kvec12, los);
	return mu * mu * G2(kvec1, kvec2);
}

double V1V1(double * kvec1, double * kvec2, double * los) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kn  = DOT(kvec12, los); 
	double result = (kn * kn) * LV1(kvec1,los) * LV1(kvec2,los);
	return result;
}

double D1V1(double * kvec1, double * kvec2, double * los) {
	
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kn  = DOT(kvec12, los); 
	
	double result = (kn) * LV1(kvec1,los) / 2.0
		          + (kn) * LV1(kvec2,los) / 2.0;
	return result;
}

/****************************/

/* third order */

double V3(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {
	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double mu = MU(kvec123, los);
	return f * mu * mu * G3(kvec1, kvec2, kvec3);

}


double V1V2(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double k123n = DOT(kvec123, los);

	double result = (f * f) * (k123n * k123n) * LV1(kvec1, los) * LV2(kvec2, kvec3, los) / 3.0
		      + (f * f) * (k123n * k123n) * LV1(kvec2, los) * LV2(kvec1, kvec3, los) / 3.0
		      + (f * f) * (k123n * k123n) * LV1(kvec3, los) * LV2(kvec1, kvec2, los) / 3.0;

	return result;

}

double V1V1V1(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double k123n = DOT(kvec123, los);

	double result = (f * f * f) * (k123n * k123n * k123n) * LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los);

	return  result; 

}

double D1V1V1(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double k123n = DOT(kvec123, los);

	double result = (f * f) * (k123n * k123n) * LV1(kvec1, los) * LV1(kvec2, los) / 3.0
		      + (f * f) * (k123n * k123n) * LV1(kvec1, los) * LV1(kvec3, los) / 3.0
		      + (f * f) * (k123n * k123n) * LV1(kvec2, los) * LV1(kvec3, los) / 3.0;

	return result;

}


double D1V2(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {
	
	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double k123n = DOT(kvec123, los);

	double result = (f) * (k123n) * LV2(kvec1, kvec2, los)  / 3.0
		      + (f) * (k123n) * LV2(kvec1, kvec3, los)  / 3.0
		      + (f) * (k123n) * LV2(kvec2, kvec3, los)  / 3.0;

	return result;

}

double D2V1(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double b1, double b2, double bK2) {

	double kvec123[3] = PLUS3(kvec1,kvec2, kvec3);
	double k123n = DOT(kvec123, los);

	double result = (f) * (k123n) * F2_Bias(kvec2,kvec3, b1, b2, bK2) * LV1(kvec1, los) / 3.0
		      + (f) * (k123n) * F2_Bias(kvec1,kvec3, b1, b2, bK2) * LV1(kvec2, los) / 3.0
		      + (f) * (k123n) * F2_Bias(kvec1,kvec2, b1, b2, bK2) * LV1(kvec3, los) / 3.0;

	return result;

}

double Z3_Bias(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO) {

	return V3(kvec1, kvec2, kvec3, los, f)
	     + V1V2(kvec1, kvec2, kvec3, los, f)
	     + V1V1V1(kvec1, kvec2, kvec3, los, f) / 6.0
	     + b1 * D1V1V1(kvec1, kvec2, kvec3, los, f) / 2.0
	     + b1 * D1V2(kvec1, kvec2, kvec3, los, f)
	     + D2V1(kvec1, kvec2, kvec3, los, f, b1, b2, bK2)
	     + F3_Bias(kvec1, kvec2, kvec3, b1, b2, b3, bK2, bK3, bDK, bO);

}

/********************************/
/********************************/
/********************************/
/* fourth order */
/********************************/
/********************************/
/********************************/

double V4(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {
	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double mu = MU(kvec1234, los);
	return f * mu * mu * G4(kvec1, kvec2, kvec3, kvec4);
}

double V2V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n = DOT(kvec1234, los);

	double result = (f * f) * (k1234n * k1234n) * LV2(kvec1, kvec2, los) * LV2(kvec3, kvec4, los)  / 3.0
		      + (f * f) * (k1234n * k1234n) * LV2(kvec1, kvec3, los) * LV2(kvec2, kvec4, los)  / 3.0
		      + (f * f) * (k1234n * k1234n) * LV2(kvec1, kvec4, los) * LV2(kvec2, kvec3, los)  / 3.0;

	return result;

}


double V1V3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n= DOT(kvec1234, los);

	double result = (f * f) * (k1234n * k1234n) * LV1(kvec1, los) * LV3(kvec2, kvec3, kvec4, los) / 4.0
		          + (f * f) * (k1234n * k1234n) * LV1(kvec2, los) * LV3(kvec1, kvec3, kvec4, los) / 4.0
		          + (f * f) * (k1234n * k1234n) * LV1(kvec3, los) * LV3(kvec1, kvec2, kvec4, los) / 4.0
		          + (f * f) * (k1234n * k1234n) * LV1(kvec4, los) * LV3(kvec1, kvec2, kvec3, los) / 4.0;

	return result;

}

double V1V1V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n= DOT(kvec1234, los);

	double result = (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec2, los) * LV2(kvec3, kvec4, los) / 6.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec3, los) * LV2(kvec2, kvec4, los) / 6.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec4, los) * LV2(kvec2, kvec3, los) / 6.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec2, los) * LV1(kvec3, los) * LV2(kvec1, kvec4, los) / 6.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec2, los) * LV1(kvec4, los) * LV2(kvec1, kvec3, los) / 6.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec3, los) * LV1(kvec4, los) * LV2(kvec1, kvec2, los) / 6.0;
	return result;

}

double V1V1V1V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n = DOT(kvec1234, los);

	double result = (f * f * f * f) * (k1234n * k1234n * k1234n * k1234n)
		      * LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec4, los); 

	return result;
}

double D1V1V1V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n = DOT(kvec1234, los);

	double result = (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los) / 4.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec4, los) / 4.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec3, los) * LV1(kvec4, los) / 4.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec4, los) / 4.0;

	return result;

}

double D2V1V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f, double b1, double b2, double bK2) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n = DOT(kvec1234, los);

	double result = (f * f) * (k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec2, los)  * F2_Bias(kvec3, kvec4, b1, b2, bK2) / 6.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec3, los)  * F2_Bias(kvec2, kvec4, b1, b2, bK2) / 6.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec4, los)  * F2_Bias(kvec2, kvec3, b1, b2, bK2) / 6.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec2, los) * LV1(kvec3, los)  * F2_Bias(kvec1, kvec4, b1, b2, bK2) / 6.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec2, los) * LV1(kvec4, los)  * F2_Bias(kvec1, kvec3, b1, b2, bK2) / 6.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec3, los) * LV1(kvec4, los)  * F2_Bias(kvec1, kvec2, b1, b2, bK2) / 6.0;

	return result;

}

double D1V1V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n= DOT(kvec1234, los);

	double result = (f * f) * (k1234n * k1234n) * LV1(kvec1, los) * LV2(kvec2, kvec3, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec1, los) * LV2(kvec2, kvec4, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec1, los) * LV2(kvec3, kvec4, los) / 12.0
	
		      + (f * f) * (k1234n * k1234n) * LV1(kvec2, los) * LV2(kvec1, kvec3, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec2, los) * LV2(kvec1, kvec4, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec2, los) * LV2(kvec3, kvec4, los) / 12.0

		      + (f * f) * (k1234n * k1234n) * LV1(kvec3, los) * LV2(kvec1, kvec2, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec3, los) * LV2(kvec1, kvec4, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec3, los) * LV2(kvec2, kvec4, los) / 12.0

		      + (f * f) * (k1234n * k1234n) * LV1(kvec4, los) * LV2(kvec1, kvec2, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec4, los) * LV2(kvec1, kvec3, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec4, los) * LV2(kvec2, kvec3, los) / 12.0;

	return result;

}

double D3V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n = DOT(kvec1234, los);

	double result = (f) * (k1234n) * LV1(kvec1, los) * F3_Bias(kvec2, kvec3, kvec4, b1, b2, b3, bK2, bK3, bDK, bO) / 4.0
		      + (f) * (k1234n) * LV1(kvec2, los) * F3_Bias(kvec1, kvec3, kvec4, b1, b2, b3, bK2, bK3, bDK, bO) / 4.0
		      + (f) * (k1234n) * LV1(kvec3, los) * F3_Bias(kvec1, kvec2, kvec4, b1, b2, b3, bK2, bK3, bDK, bO) / 4.0
		      + (f) * (k1234n) * LV1(kvec4, los) * F3_Bias(kvec1, kvec2, kvec3, b1, b2, b3, bK2, bK3, bDK, bO) / 4.0;

	return result;

}

double D2V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f, double b1, double b2, double bK2) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n = DOT(kvec1234, los);

	double result = (f) * (k1234n) * LV2(kvec1, kvec2, los) * F2_Bias(kvec3, kvec4, b1, b2, bK2) / 6.0
		      + (f) * (k1234n) * LV2(kvec1, kvec3, los) * F2_Bias(kvec2, kvec4, b1, b2, bK2) / 6.0
		      + (f) * (k1234n) * LV2(kvec1, kvec4, los) * F2_Bias(kvec2, kvec3, b1, b2, bK2) / 6.0
		      + (f) * (k1234n) * LV2(kvec2, kvec3, los) * F2_Bias(kvec1, kvec4, b1, b2, bK2) / 6.0
		      + (f) * (k1234n) * LV2(kvec2, kvec4, los) * F2_Bias(kvec1, kvec3, b1, b2, bK2) / 6.0
		      + (f) * (k1234n) * LV2(kvec3, kvec4, los) * F2_Bias(kvec1, kvec2, b1, b2, bK2) / 6.0;

	return result;

}

double D1V3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n = DOT(kvec1234, los);

	double result = (f) * (k1234n) * LV3(kvec1, kvec2, kvec3, los) / 4.0
	              + (f) * (k1234n) * LV3(kvec1, kvec2, kvec4, los) / 4.0
	              + (f) * (k1234n) * LV3(kvec1, kvec3, kvec4, los) / 4.0
	              + (f) * (k1234n) * LV3(kvec2, kvec3, kvec4, los) / 4.0;

	return result;

}

double Z4_Bias(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO) {

	return V4(kvec1, kvec2, kvec3, kvec4, los, f)
	     + V2V2(kvec1, kvec2, kvec3, kvec4, los, f) / 2.0
	     + V1V3(kvec1, kvec2, kvec3, kvec4, los, f)
	     + V1V1V2(kvec1, kvec2, kvec3, kvec4, los, f) / 2.0
	     + V1V1V1V1(kvec1, kvec2, kvec3, kvec4, los, f) / 24.0

	     + b1 * D1V1V1V1(kvec1, kvec2, kvec3, kvec4, los, f) / 6.0
	     + b1 * D1V1V2(kvec1, kvec2, kvec3, kvec4, los, f) 
	     + b1 * D1V3(kvec1, kvec2, kvec3, kvec4, los, f) 

	     + D2V1V1(kvec1, kvec2, kvec3, kvec4, los, f, b1, b2, bK2) / 2.0
	     + D2V2(kvec1, kvec2, kvec3, kvec4, los, f, b1, b2, bK2) 
	     + D3V1(kvec1, kvec2, kvec3, kvec4, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 

	     + F4_Bias(kvec1, kvec2, kvec3, kvec4, b1);

}
/********************************/
/********************************/
/********************************/
/* fifth order */
/********************************/
/********************************/
/********************************/

double V5(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {
	double kvec12345[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double mu = MU(kvec12345, los);
	
	return (f) * (mu * mu) * G5(kvec1, kvec2, kvec3, kvec4, kvec5);
}

double V1V4(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {
	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV4(kvec2, kvec3, kvec4, kvec5, los)
		      + LV1(kvec2, los) * LV4(kvec1, kvec3, kvec4, kvec5, los)
		      + LV1(kvec3, los) * LV4(kvec1, kvec2, kvec4, kvec5, los)
		      + LV1(kvec4, los) * LV4(kvec1, kvec2, kvec3, kvec5, los)
		      + LV1(kvec5, los) * LV4(kvec1, kvec2, kvec3, kvec4, los);

	return pow(f, 2) * pow(kn, 2) * result / 5.0;

}

double V2V3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {
	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV2(kvec1, kvec2, los) * LV3(kvec3, kvec4, kvec5, los)
	              + LV2(kvec1, kvec3, los) * LV3(kvec2, kvec4, kvec5, los)
	              + LV2(kvec1, kvec4, los) * LV3(kvec2, kvec3, kvec5, los)
	              + LV2(kvec1, kvec5, los) * LV3(kvec2, kvec3, kvec4, los)
	              + LV2(kvec2, kvec3, los) * LV3(kvec1, kvec4, kvec5, los)
	              + LV2(kvec2, kvec4, los) * LV3(kvec1, kvec3, kvec5, los)
	              + LV2(kvec2, kvec5, los) * LV3(kvec1, kvec3, kvec4, los)
	              + LV2(kvec3, kvec4, los) * LV3(kvec1, kvec2, kvec5, los)
	              + LV2(kvec3, kvec5, los) * LV3(kvec1, kvec2, kvec4, los)
	              + LV2(kvec4, kvec5, los) * LV3(kvec1, kvec2, kvec3, los);

	return pow(f, 2) * pow(kn, 2) * result / 10.0;

}

double V1V2V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {
	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV2(kvec2, kvec3, los)  * LV2(kvec4, kvec5, los)
		      + LV1(kvec1, los) * LV2(kvec2, kvec4, los)  * LV2(kvec3, kvec5, los)
		      + LV1(kvec1, los) * LV2(kvec2, kvec5, los)  * LV2(kvec3, kvec4, los)

		      + LV1(kvec2, los) * LV2(kvec1, kvec3, los)  * LV2(kvec4, kvec5, los)
		      + LV1(kvec2, los) * LV2(kvec1, kvec4, los)  * LV2(kvec3, kvec5, los)
		      + LV1(kvec2, los) * LV2(kvec1, kvec5, los)  * LV2(kvec3, kvec4, los)

		      + LV1(kvec3, los) * LV2(kvec1, kvec2, los)  * LV2(kvec4, kvec5, los)
		      + LV1(kvec3, los) * LV2(kvec1, kvec4, los)  * LV2(kvec2, kvec5, los)
		      + LV1(kvec3, los) * LV2(kvec1, kvec5, los)  * LV2(kvec2, kvec4, los)

		      + LV1(kvec4, los) * LV2(kvec1, kvec2, los)  * LV2(kvec3, kvec5, los)
		      + LV1(kvec4, los) * LV2(kvec1, kvec3, los)  * LV2(kvec2, kvec5, los)
		      + LV1(kvec4, los) * LV2(kvec1, kvec5, los)  * LV2(kvec2, kvec3, los)

		      + LV1(kvec5, los) * LV2(kvec1, kvec2, los)  * LV2(kvec3, kvec4, los)
		      + LV1(kvec5, los) * LV2(kvec1, kvec3, los)  * LV2(kvec2, kvec4, los)
		      + LV1(kvec5, los) * LV2(kvec1, kvec4, los)  * LV2(kvec2, kvec3, los);

	return pow(f, 3) * pow(kn, 3) * result / 15.0;

}

double V1V1V3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV1(kvec2, los) * LV3(kvec3, kvec4, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV3(kvec2, kvec4, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec4, los) * LV3(kvec2, kvec3, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec5, los) * LV3(kvec2, kvec3, kvec4, los)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV3(kvec1, kvec4, kvec5, los)
		      + LV1(kvec2, los) * LV1(kvec4, los) * LV3(kvec1, kvec3, kvec5, los)
		      + LV1(kvec2, los) * LV1(kvec5, los) * LV3(kvec1, kvec3, kvec4, los)
		      + LV1(kvec3, los) * LV1(kvec4, los) * LV3(kvec1, kvec2, kvec5, los)
		      + LV1(kvec3, los) * LV1(kvec5, los) * LV3(kvec1, kvec2, kvec4, los)
		      + LV1(kvec4, los) * LV1(kvec5, los) * LV3(kvec1, kvec2, kvec3, los);

	return pow(f, 3) * pow(kn, 3) * result / 10.0;

}

double V1V1V1V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec3, los) * LV1(kvec4, los) * LV1(kvec5, los) * LV2(kvec1, kvec2, los)
		      + LV1(kvec2, los) * LV1(kvec4, los) * LV1(kvec5, los) * LV2(kvec1, kvec3, los)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec5, los) * LV2(kvec1, kvec4, los)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec4, los) * LV2(kvec1, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec4, los) * LV1(kvec5, los) * LV2(kvec2, kvec3, los)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV1(kvec5, los) * LV2(kvec2, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV1(kvec4, los) * LV2(kvec2, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec5, los) * LV2(kvec3, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec4, los) * LV2(kvec3, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los) * LV2(kvec4, kvec5, los);

	return pow(f, 4) * pow(kn, 4) * result / 10.0;

}

double V1V1V1V1V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec4, los) * LV1(kvec5, los);

	return pow(f, 5) * pow(kn, 5) * result;

}

double D1V1V1V1V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec4, los)
	              + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec5, los)
	              + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec4, los) * LV1(kvec5, los)
	              + LV1(kvec1, los) * LV1(kvec3, los) * LV1(kvec4, los) * LV1(kvec5, los)
	              + LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec4, los) * LV1(kvec5, los);

	return pow(f, 4) * pow(kn, 4) * result / 5.0;
}

double D1V1V1V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV1(kvec2, los) * LV2(kvec3, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV2(kvec2, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec4, los) * LV2(kvec2, kvec3, los)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV2(kvec1, kvec4, los)
		      + LV1(kvec2, los) * LV1(kvec4, los) * LV2(kvec1, kvec3, los)
		      + LV1(kvec3, los) * LV1(kvec4, los) * LV2(kvec1, kvec2, los)

		      + LV1(kvec1, los) * LV1(kvec2, los) * LV2(kvec3, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV2(kvec2, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec5, los) * LV2(kvec2, kvec3, los)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV2(kvec1, kvec5, los)
		      + LV1(kvec2, los) * LV1(kvec5, los) * LV2(kvec1, kvec3, los)
		      + LV1(kvec3, los) * LV1(kvec5, los) * LV2(kvec1, kvec2, los)

		      + LV1(kvec1, los) * LV1(kvec2, los) * LV2(kvec5, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec5, los) * LV2(kvec2, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec4, los) * LV2(kvec2, kvec5, los)
		      + LV1(kvec2, los) * LV1(kvec5, los) * LV2(kvec1, kvec4, los)
		      + LV1(kvec2, los) * LV1(kvec4, los) * LV2(kvec1, kvec5, los)
		      + LV1(kvec5, los) * LV1(kvec4, los) * LV2(kvec1, kvec2, los)

		      + LV1(kvec1, los) * LV1(kvec5, los) * LV2(kvec3, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV2(kvec5, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec4, los) * LV2(kvec5, kvec3, los)
		      + LV1(kvec5, los) * LV1(kvec3, los) * LV2(kvec1, kvec4, los)
		      + LV1(kvec5, los) * LV1(kvec4, los) * LV2(kvec1, kvec3, los)
		      + LV1(kvec3, los) * LV1(kvec4, los) * LV2(kvec1, kvec5, los)

		      + LV1(kvec5, los) * LV1(kvec2, los) * LV2(kvec3, kvec4, los)
		      + LV1(kvec5, los) * LV1(kvec3, los) * LV2(kvec2, kvec4, los)
		      + LV1(kvec5, los) * LV1(kvec4, los) * LV2(kvec2, kvec3, los)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV2(kvec5, kvec4, los)
		      + LV1(kvec2, los) * LV1(kvec4, los) * LV2(kvec5, kvec3, los)
		      + LV1(kvec3, los) * LV1(kvec4, los) * LV2(kvec5, kvec2, los);

	return pow(f, 3) * pow(kn, 3) * result / 30.0;

}

double D1V1V3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV3(kvec2, kvec3, kvec4, los)
		      + LV1(kvec2, los) * LV3(kvec1, kvec3, kvec4, los)
		      + LV1(kvec3, los) * LV3(kvec1, kvec2, kvec4, los)
		      + LV1(kvec4, los) * LV3(kvec1, kvec2, kvec3, los)

		      + LV1(kvec1, los) * LV3(kvec2, kvec3, kvec5, los)
		      + LV1(kvec2, los) * LV3(kvec1, kvec3, kvec5, los)
		      + LV1(kvec3, los) * LV3(kvec1, kvec2, kvec5, los)
		      + LV1(kvec5, los) * LV3(kvec1, kvec2, kvec3, los)

		      + LV1(kvec1, los) * LV3(kvec2, kvec5, kvec4, los)
		      + LV1(kvec2, los) * LV3(kvec1, kvec5, kvec4, los)
		      + LV1(kvec5, los) * LV3(kvec1, kvec2, kvec4, los)
		      + LV1(kvec4, los) * LV3(kvec1, kvec2, kvec5, los)

		      + LV1(kvec1, los) * LV3(kvec5, kvec3, kvec4, los)
		      + LV1(kvec5, los) * LV3(kvec1, kvec3, kvec4, los)
		      + LV1(kvec3, los) * LV3(kvec1, kvec5, kvec4, los)
		      + LV1(kvec4, los) * LV3(kvec1, kvec5, kvec3, los)

		      + LV1(kvec5, los) * LV3(kvec2, kvec3, kvec4, los)
		      + LV1(kvec2, los) * LV3(kvec5, kvec3, kvec4, los)
		      + LV1(kvec3, los) * LV3(kvec5, kvec2, kvec4, los)
		      + LV1(kvec4, los) * LV3(kvec5, kvec2, kvec3, los);

	return pow(f, 2) * pow(kn, 2) * result / 20.0;
}


double D1V2V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV2(kvec1, kvec2, los) * LV2(kvec3, kvec4, los)
		      + LV2(kvec1, kvec3, los) * LV2(kvec2, kvec4, los)
		      + LV2(kvec1, kvec4, los) * LV2(kvec2, kvec3, los)

		      + LV2(kvec1, kvec2, los) * LV2(kvec3, kvec5, los)
		      + LV2(kvec1, kvec3, los) * LV2(kvec2, kvec5, los)
		      + LV2(kvec1, kvec5, los) * LV2(kvec2, kvec3, los)

		      + LV2(kvec1, kvec2, los) * LV2(kvec5, kvec4, los)
		      + LV2(kvec1, kvec5, los) * LV2(kvec2, kvec4, los)
		      + LV2(kvec1, kvec4, los) * LV2(kvec2, kvec5, los)

		      + LV2(kvec1, kvec5, los) * LV2(kvec3, kvec4, los)
		      + LV2(kvec1, kvec3, los) * LV2(kvec5, kvec4, los)
		      + LV2(kvec1, kvec4, los) * LV2(kvec5, kvec3, los)

		      + LV2(kvec5, kvec2, los) * LV2(kvec3, kvec4, los)
		      + LV2(kvec5, kvec3, los) * LV2(kvec2, kvec4, los)
		      + LV2(kvec5, kvec4, los) * LV2(kvec2, kvec3, los);

	return pow(f, 2) * pow(kn, 2) * result / 15.0;

}


double D1V4(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {
	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV4(kvec1, kvec2, kvec3, kvec4, los)
		      + LV4(kvec1, kvec2, kvec3, kvec5, los)
		      + LV4(kvec1, kvec2, kvec5, kvec4, los)
		      + LV4(kvec1, kvec5, kvec3, kvec4, los)
	              + LV4(kvec5, kvec2, kvec3, kvec4, los);

	return pow(f, 1) * pow(kn, 1) * result / 5.0;

}

double D2V1V1V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f, double b1, double b2, double bK2) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec3, los) * LV1(kvec4, los) * LV1(kvec5, los) * F2_Bias(kvec1, kvec2, b1, b2, bK2)
		      + LV1(kvec2, los) * LV1(kvec4, los) * LV1(kvec5, los) * F2_Bias(kvec1, kvec3, b1, b2, bK2)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec5, los) * F2_Bias(kvec1, kvec4, b1, b2, bK2)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec4, los) * F2_Bias(kvec1, kvec5, b1, b2, bK2)
		      + LV1(kvec1, los) * LV1(kvec4, los) * LV1(kvec5, los) * F2_Bias(kvec2, kvec3, b1, b2, bK2)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV1(kvec5, los) * F2_Bias(kvec2, kvec4, b1, b2, bK2)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV1(kvec4, los) * F2_Bias(kvec2, kvec5, b1, b2, bK2)
		      + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec5, los) * F2_Bias(kvec3, kvec4, b1, b2, bK2)
		      + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec4, los) * F2_Bias(kvec3, kvec5, b1, b2, bK2)
		      + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los) * F2_Bias(kvec4, kvec5, b1, b2, bK2);

	return pow(f, 3) * pow(kn, 3) * result / 10.0;
}

double D2V1V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f, double b1, double b2, double bK2) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV2(kvec2, kvec3, los) * F2_Bias(kvec4, kvec5, b1, b2, bK2)
		      + LV1(kvec1, los) * LV2(kvec2, kvec4, los) * F2_Bias(kvec3, kvec5, b1, b2, bK2)
		      + LV1(kvec1, los) * LV2(kvec2, kvec5, los) * F2_Bias(kvec3, kvec4, b1, b2, bK2)
		      + LV1(kvec1, los) * LV2(kvec3, kvec4, los) * F2_Bias(kvec2, kvec5, b1, b2, bK2)
		      + LV1(kvec1, los) * LV2(kvec3, kvec5, los) * F2_Bias(kvec2, kvec4, b1, b2, bK2)
		      + LV1(kvec1, los) * LV2(kvec4, kvec5, los) * F2_Bias(kvec2, kvec3, b1, b2, bK2)

		      + LV1(kvec2, los) * LV2(kvec1, kvec3, los) * F2_Bias(kvec4, kvec5, b1, b2, bK2)
		      + LV1(kvec2, los) * LV2(kvec1, kvec4, los) * F2_Bias(kvec3, kvec5, b1, b2, bK2)
		      + LV1(kvec2, los) * LV2(kvec1, kvec5, los) * F2_Bias(kvec3, kvec4, b1, b2, bK2)
		      + LV1(kvec2, los) * LV2(kvec3, kvec4, los) * F2_Bias(kvec1, kvec5, b1, b2, bK2)
		      + LV1(kvec2, los) * LV2(kvec3, kvec5, los) * F2_Bias(kvec1, kvec4, b1, b2, bK2)
		      + LV1(kvec2, los) * LV2(kvec4, kvec5, los) * F2_Bias(kvec1, kvec3, b1, b2, bK2)

		      + LV1(kvec3, los) * LV2(kvec2, kvec1, los) * F2_Bias(kvec4, kvec5, b1, b2, bK2)
		      + LV1(kvec3, los) * LV2(kvec2, kvec4, los) * F2_Bias(kvec1, kvec5, b1, b2, bK2)
		      + LV1(kvec3, los) * LV2(kvec2, kvec5, los) * F2_Bias(kvec1, kvec4, b1, b2, bK2)
		      + LV1(kvec3, los) * LV2(kvec1, kvec4, los) * F2_Bias(kvec2, kvec5, b1, b2, bK2)
		      + LV1(kvec3, los) * LV2(kvec1, kvec5, los) * F2_Bias(kvec2, kvec4, b1, b2, bK2)
		      + LV1(kvec3, los) * LV2(kvec4, kvec5, los) * F2_Bias(kvec2, kvec1, b1, b2, bK2)

		      + LV1(kvec4, los) * LV2(kvec2, kvec3, los) * F2_Bias(kvec1, kvec5, b1, b2, bK2)
		      + LV1(kvec4, los) * LV2(kvec2, kvec1, los) * F2_Bias(kvec3, kvec5, b1, b2, bK2)
		      + LV1(kvec4, los) * LV2(kvec2, kvec5, los) * F2_Bias(kvec3, kvec1, b1, b2, bK2)
		      + LV1(kvec4, los) * LV2(kvec3, kvec1, los) * F2_Bias(kvec2, kvec5, b1, b2, bK2)
		      + LV1(kvec4, los) * LV2(kvec3, kvec5, los) * F2_Bias(kvec2, kvec1, b1, b2, bK2)
		      + LV1(kvec4, los) * LV2(kvec1, kvec5, los) * F2_Bias(kvec2, kvec3, b1, b2, bK2)

		      + LV1(kvec5, los) * LV2(kvec2, kvec3, los) * F2_Bias(kvec4, kvec1, b1, b2, bK2)
		      + LV1(kvec5, los) * LV2(kvec2, kvec4, los) * F2_Bias(kvec3, kvec1, b1, b2, bK2)
		      + LV1(kvec5, los) * LV2(kvec2, kvec1, los) * F2_Bias(kvec3, kvec4, b1, b2, bK2)
		      + LV1(kvec5, los) * LV2(kvec3, kvec4, los) * F2_Bias(kvec2, kvec1, b1, b2, bK2)
		      + LV1(kvec5, los) * LV2(kvec3, kvec1, los) * F2_Bias(kvec2, kvec4, b1, b2, bK2)
		      + LV1(kvec5, los) * LV2(kvec4, kvec1, los) * F2_Bias(kvec2, kvec3, b1, b2, bK2);

	return pow(f, 2) * pow(kn, 2) * result / 30.0;
}


double D2V3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f, double b1, double b2, double bK2) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = F2_Bias(kvec1, kvec2, b1, b2, bK2) * LV3(kvec3, kvec4, kvec5, los)
	              + F2_Bias(kvec1, kvec3, b1, b2, bK2) * LV3(kvec2, kvec4, kvec5, los)
	              + F2_Bias(kvec1, kvec4, b1, b2, bK2) * LV3(kvec2, kvec3, kvec5, los)
	              + F2_Bias(kvec1, kvec5, b1, b2, bK2) * LV3(kvec2, kvec3, kvec4, los)
	              + F2_Bias(kvec2, kvec3, b1, b2, bK2) * LV3(kvec1, kvec4, kvec5, los)
	              + F2_Bias(kvec2, kvec4, b1, b2, bK2) * LV3(kvec1, kvec3, kvec5, los)
	              + F2_Bias(kvec2, kvec5, b1, b2, bK2) * LV3(kvec1, kvec3, kvec4, los)
	              + F2_Bias(kvec3, kvec4, b1, b2, bK2) * LV3(kvec1, kvec2, kvec5, los)
	              + F2_Bias(kvec3, kvec5, b1, b2, bK2) * LV3(kvec1, kvec2, kvec4, los)
	              + F2_Bias(kvec4, kvec5, b1, b2, bK2) * LV3(kvec1, kvec2, kvec3, los);

	return pow(f, 1) * pow(kn, 1) * result / 10.0;

}

double D3V1V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV1(kvec2, los) * F3_Bias(kvec3, kvec4, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec1, los) * LV1(kvec3, los) * F3_Bias(kvec2, kvec4, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec1, los) * LV1(kvec4, los) * F3_Bias(kvec2, kvec3, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec1, los) * LV1(kvec5, los) * F3_Bias(kvec2, kvec3, kvec4, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec2, los) * LV1(kvec3, los) * F3_Bias(kvec1, kvec4, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec2, los) * LV1(kvec4, los) * F3_Bias(kvec1, kvec3, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec2, los) * LV1(kvec5, los) * F3_Bias(kvec1, kvec3, kvec4, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec3, los) * LV1(kvec4, los) * F3_Bias(kvec1, kvec2, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec3, los) * LV1(kvec5, los) * F3_Bias(kvec1, kvec2, kvec4, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec4, los) * LV1(kvec5, los) * F3_Bias(kvec1, kvec2, kvec3, b1, b2, b3, bK2, bK3, bDK, bO);

	return pow(f, 2) * pow(kn, 2) * result / 10.0;


}

double D3V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV2(kvec1, kvec2, los) * F3_Bias(kvec3, kvec4, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec1, kvec3, los) * F3_Bias(kvec2, kvec4, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec1, kvec4, los) * F3_Bias(kvec2, kvec3, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec1, kvec5, los) * F3_Bias(kvec2, kvec3, kvec4, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec2, kvec3, los) * F3_Bias(kvec1, kvec4, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec2, kvec4, los) * F3_Bias(kvec1, kvec3, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec2, kvec5, los) * F3_Bias(kvec1, kvec3, kvec4, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec3, kvec4, los) * F3_Bias(kvec1, kvec2, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec3, kvec5, los) * F3_Bias(kvec1, kvec2, kvec4, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec4, kvec5, los) * F3_Bias(kvec1, kvec2, kvec3, b1, b2, b3, bK2, bK3, bDK, bO);

	return pow(f, 1) * pow(kn, 1) * result / 10.0;

}

double D4V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f, double b1) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * F4_Bias(kvec2, kvec3, kvec4, kvec5, b1)
		      + LV1(kvec2, los) * F4_Bias(kvec1, kvec3, kvec4, kvec5, b1)
		      + LV1(kvec3, los) * F4_Bias(kvec1, kvec2, kvec4, kvec5, b1)
		      + LV1(kvec4, los) * F4_Bias(kvec1, kvec2, kvec3, kvec5, b1)
		      + LV1(kvec5, los) * F4_Bias(kvec1, kvec2, kvec3, kvec4, b1);

	return pow(f, 1) * pow(kn, 1) * result / 5.0;

}

double Z5_Bias(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f, 
	       double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO) {

	double Z = V5(kvec1, kvec2, kvec3, kvec4, kvec5, los, f)
		 + V1V4(kvec1, kvec2, kvec3, kvec4, kvec5, los, f)
		 + V2V3(kvec1, kvec2, kvec3, kvec4, kvec5, los, f)
		 + V1V2V2(kvec1, kvec2, kvec3, kvec4, kvec5, los, f) / 2.0
		 + V1V1V3(kvec1, kvec2, kvec3, kvec4, kvec5, los, f) / 2.0
		 + V1V1V1V2(kvec1, kvec2, kvec3, kvec4, kvec5, los, f) / 6.0
		 + V1V1V1V1V1(kvec1, kvec2, kvec3, kvec4, kvec5, los, f) / 120.0

		 + b1 * D1V1V1V1V1(kvec1, kvec2, kvec3, kvec4, kvec5, los, f) / 24.0
		 + b1 * D1V1V1V2(kvec1, kvec2, kvec3, kvec4, kvec5, los, f) / 2.0
		 + b1 * D1V1V3(kvec1, kvec2, kvec3, kvec4, kvec5, los, f)
		 + b1 * D1V2V2(kvec1, kvec2, kvec3, kvec4, kvec5, los, f) / 2.0
		 + b1 * D1V4(kvec1, kvec2, kvec3, kvec4, kvec5, los, f)
		 
		 + D2V1V1V1(kvec1, kvec2, kvec3, kvec4, kvec5, los, f, b1, b2, bK2) / 6.0
		 + D2V1V2(kvec1, kvec2, kvec3, kvec4, kvec5, los, f, b1, b2, bK2) / 6.0
		 + D2V3(kvec1, kvec2, kvec3, kvec4, kvec5, los, f, b1, b2, bK2) 

		 + D3V1V1(kvec1, kvec2, kvec3, kvec4, kvec5, los, f, b1, b2, b3, bK2, bK3, bDK, bO) / 2.0
		 + D3V2(kvec1, kvec2, kvec3, kvec4, kvec5, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 

		 + D4V1(kvec1, kvec2, kvec3, kvec4, kvec5, los, f, b1) 

		 + F5_Bias(kvec1, kvec2, kvec3, kvec4, kvec5, b1);

	return Z;

}

/****************************/
/*     Reconstruction       */
/****************************/

double Z1S1(
        double * kvec1, double * kvec2, double * los, 
        double f, double b1, 
        double one_over_b1_fid, double R) {

    double k1 = NORM(kvec1);
    double k2 = NORM(kvec2);
    double W1 = exp( - pow(k1 * R, 2) / 2.0);
    double W2 = exp( - pow(k2 * R, 2) / 2.0);
    double kvec12[3] = PLUS(kvec1, kvec2);
    double result = (-1.0/2.0) * (one_over_b1_fid) * ( LV1(kvec1, kvec12) * W1 + LV1(kvec2, kvec12) * W2 ) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1);
    return result;
}

double Z1S1_b1_b1(double * kvec1, double * kvec2, double one_over_b1_fid, double R) {
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
    	double W1 = exp( - pow(k1 * R, 2) / 2.0);
	double W2 = exp( - pow(k2 * R, 2) / 2.0);
	double kvec12[3] = PLUS(kvec1, kvec2);
	double result = (-1.0/2.0) * (one_over_b1_fid) * ( LV1(kvec1, kvec12) * W1 + LV1(kvec2, kvec12) * W2 ) * D1() * D1();
	return result;
}

double Z1S1_b1_f(double * kvec1, double * kvec2, double * los, double one_over_b1_fid, double R) {
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
    	double W1 = exp( - pow(k1 * R, 2) / 2.0);
	double W2 = exp( - pow(k2 * R, 2) / 2.0);
	double kvec12[3] = PLUS(kvec1, kvec2);
	double result = (-1.0/2.0) * (one_over_b1_fid) * ( LV1(kvec1, kvec12) * W1 + LV1(kvec2, kvec12) * W2 ) 
	              * ( D1() * V1(kvec1, los) + V1(kvec2, los) * D1() );
	return result;
}

double Z1S1_f_f(double * kvec1, double * kvec2, double * los, double one_over_b1_fid, double R) {
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
    double W1 = exp( - pow(k1 * R, 2) / 2.0);
	double W2 = exp( - pow(k2 * R, 2) / 2.0);
	double kvec12[3] = PLUS(kvec1, kvec2);
	double result = (-1.0/2.0) * (one_over_b1_fid) * ( LV1(kvec1, kvec12) * W1 + LV1(kvec2, kvec12) * W2 ) 
	              * V1(kvec1, los) * V1(kvec2, los);
	return result;
}


double Z2_Bias_Reconstructed(
        double * kvec1, double * kvec2, double * los, 
        double f, double b1, double b2, double bK2, 
        double one_over_b1_fid, double R) {
	return Z2_Bias(kvec1, kvec2, los, f, b1, b2, bK2) 
	     + Z1S1(kvec1, kvec2, los, f, b1, one_over_b1_fid, R);
} 

double Z3_Bias_Reconstructed(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, double one_over_b1_fid, double R) {

	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec23[3] = PLUS(kvec2, kvec3);
	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double k12 = NORM(kvec12);
	double k13 = NORM(kvec13);
	double k23 = NORM(kvec23);

    double W1 = exp( - pow(k1 * R, 2) / 2.0);
	double W2 = exp( - pow(k2 * R, 2) / 2.0);
	double W3 = exp( - pow(k3 * R, 2) / 2.0);

    double W12 = exp( - pow(k12 * R, 2) / 2.0);
	double W23 = exp( - pow(k23 * R, 2) / 2.0);
	double W13 = exp( - pow(k13 * R, 2) / 2.0);

	double Z3 = Z3_Bias(kvec1, kvec2, kvec3, los, f, b1, b2, b3, bK2, bK3, bDK, bO);

	double Z1_1 = Z1_Bias(kvec1, los, f, b1);
	double Z1_2 = Z1_Bias(kvec2, los, f, b1);
	double Z1_3 = Z1_Bias(kvec3, los, f, b1);

	double Z2_12 = Z2_Bias(kvec1, kvec2, los, f, b1, b2, bK2);
	double Z2_13 = Z2_Bias(kvec1, kvec3, los, f, b1, b2, bK2);
	double Z2_23 = Z2_Bias(kvec2, kvec3, los, f, b1, b2, bK2);

	double R_A = ( LV1(kvec23, kvec123) * (W23 * one_over_b1_fid) 
                 + LV1(kvec1, kvec123) *  (W1 * one_over_b1_fid) ) * Z1_1 * Z2_23
		       + ( LV1(kvec13, kvec123) * (W13 * one_over_b1_fid)
                 + LV1(kvec2, kvec123) *  (W2 * one_over_b1_fid) ) * Z1_2 * Z2_13
		       + ( LV1(kvec12, kvec123) * (W12 * one_over_b1_fid) 
                 + LV1(kvec3, kvec123) *  (W3 * one_over_b1_fid) ) * Z1_3 * Z2_12;

	double R_B = ( LV1(kvec1, kvec123) * (W1 * one_over_b1_fid) * LV1(kvec2, kvec123) * (W2 * one_over_b1_fid)
		       +   LV1(kvec1, kvec123) * (W1 * one_over_b1_fid) * LV1(kvec3, kvec123) * (W3 * one_over_b1_fid)
		       +   LV1(kvec2, kvec123) * (W2 * one_over_b1_fid) * LV1(kvec3, kvec123) * (W3 * one_over_b1_fid) )
               * Z1_1 * Z1_2 * Z1_3;

	return Z3 - (1.0/3.0) * R_A + (1.0/6.0) * R_B;
}

double ExpDamping(double * kvec, double * los, double sigma2_perp, double sigma2_para ) {

	double k = NORM(kvec);
	double mu = MU(kvec, los);
	double mu2 = mu * mu;
	double lnD = - k * k * ( (1.0 - mu2) * sigma2_perp + mu2 * sigma2_para ) / 2.0;
	double D = exp(lnD);
	return D;

}

double ExpDamping_RealSpace(double * kvec, double sigma2_perp) {

	double k = NORM(kvec);
	double lnD = - k * k * ( sigma2_perp ) / 2.0;
	double D = exp(lnD);
	return D;

}

#endif

