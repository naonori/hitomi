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
#include <set>
#include <map>
#include <vector>
#include <random>
#include <algorithm>

/* standard tools */
#include <omp.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_bessel.h>
#include <fftw3.h>

/* inner product */
double DOT(double * x, double * y) {
	return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}

/* norm */
double NORM(double * x) {
	return std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
}

/* plus */
#define PLUS(x, y) {x[0]+y[0], x[1]+y[1], x[2]+y[2]}
#define PLUS3(x, y, z) {x[0]+y[0]+z[0], x[1]+y[1]+z[1], x[2]+y[2]+z[2]}
#define PLUS4(x, y, z, l) {x[0]+y[0]+z[0]+l[0], x[1]+y[1]+z[1]+l[1], x[2]+y[2]+z[2]+l[2]}
#define PLUS5(x, y, z, l, m) {x[0]+y[0]+z[0]+l[0]+m[0], x[1]+y[1]+z[1]+l[1]+m[1], x[2]+y[2]+z[2]+l[2]+m[2]}
#define MINUS(x, y) {x[0]-y[0], x[1]-y[1], x[2]-y[2]}
#define MINUS_XVEC(x) { - x[0], - x[1], - x[2]}

#endif
