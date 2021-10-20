#ifndef __common__
#define __common__

/* standard libraries */
#include <iostream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <complex>
#include <cmath>
#include <ctime>
#include <vector>

/* tooles */
#include <fftw3.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_coupling.h>

#define wigner_3j(j1,j2,j3,m1,m2,m3) (gsl_sf_coupling_3j(2*j1,2*j2,2*j3,2*m1,2*m2,2*m3)) 

#ifdef __MPI__
#include <mpi.h>
#endif

/* memory */
double byte = 0;

/* MPI */
int ThisTask = 0;
int NTask = 1;

/* time */
double start = 0.0;
double timesec = 0.0;

void EXIT() {

    #ifdef __MPI__
	MPI_Finalize();
    #endif
    exit(1);
}

#endif

