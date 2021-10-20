#ifndef __bessel__
#define __bessel__

class SphericalBesselClass {
private:
	gsl_interp_accel *acc;
	gsl_spline       *spline;
public:

	SphericalBesselClass(int ell) {

		double xmin = 0;
		double xmax = 10000;
		int NN = 1000000;
		double dx = (xmax - xmin) / double(NN-1);
		double * x = new double[NN];
		double * bessel = new double[NN];
		for(int i = 0; i < NN; i++) {
			x[i] = xmin + dx * double(i);
			bessel[i] = gsl_sf_bessel_jl(ell, x[i]);
		}

		this->acc = gsl_interp_accel_alloc();
		this->spline = gsl_spline_alloc(gsl_interp_cspline, NN);
		gsl_spline_init(this->spline, x, bessel, NN);
		delete [] x; x = NULL;
		delete [] bessel; bessel = NULL;
	}

	~SphericalBesselClass() {
		if(this->acc != NULL) { gsl_interp_accel_free(this->acc); this->acc = NULL; }
		if(this->spline != NULL) { gsl_spline_free(this->spline); this->spline = NULL; }
	}

	double getSphericalBessel(double x) {
		return gsl_spline_eval(this->spline, x, this->acc);
	}

};


#endif

