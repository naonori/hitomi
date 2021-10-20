#ifndef __powerspec__
#define __powerspec__

int calcPowerSpectrum(
        ParticleBOSSClass & P_D, ParticleBOSSClass & P_R, 
        ParameterClass & param, double alpha, double * kbin, double Vsurvey) {

    /*************************************************************************************************/
    /* This function calculates the power spectrum multipoles, expanded by the Legendre function, using FFTs.*/
    /*************************************************************************************************/

	if( !( (param.ELL == param.ell1) && (param.ell2 == 0) ) ) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
		return -1;
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleBOSSClass> dn_00(param);
	dn_00.calcYlmWeightedDensityFluctuation(P_D, P_R, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calcFourierTransform();

	/* define pk_save */
	std::complex<double> * pk_save = new std::complex<double>[param.n_kbin];
	for(int i =0; i < param.n_kbin; i++) { 
		pk_save[i] = 0.0; 
	}

	/* calc power spectrum */
	for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/****************************************/
		/* calc the ylm-weighted density fluctuation */
		DensityFieldClass<ParticleBOSSClass> dn_LM(param);
		dn_LM.calcYlmWeightedDensityFluctuation(P_D, P_R, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM.calcFourierTransform();

		/* calc shotnoise term */
		TwoPointStatisticsClass<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calcShotNoiseForPowerspectrum(P_D, P_R, alpha, param.ELL, _M_);
		
	    for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {

	    	/*****/
	    	/* (-1)**m1 delta_{m1, -M} */
	    	double w = wigner_3j(param.ell1, 0, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, 0, param.ELL, _m1_, 0, _M_);
	    	w *= double(2*param.ELL+1) * double(2*param.ell1+1);
	    	if(fabs(w) < 1.0e-10) {
	    		continue;
	    	}

    		stat.calcPowerSpectrum(dn_LM, dn_00, kbin, shotnoise, param.ell1, _m1_);

		    for(int i =0; i < param.n_kbin; i++) { 
			    pk_save[i] += w * stat.pk[i];
	       	}
	    }
		
        timesec = double(clock() - start);
		if(ThisTask == 0) {
			printf("M = %d | %.3f sec\n", _M_, timesec / CLOCKS_PER_SEC);
		}
	}

    /* calc the normalization factor */
	double norm = ParticleBOSSClass::calcNormalizationForPowerSpectrum(P_D, Vsurvey);
        
    /* save data */
	FILE * fp;
	char buf[1024];
    if(0) {
    } else if (param.flag_recon == "True") {
	    sprintf(buf, "%s/pk%d_recon.dat", param.output_dir.c_str(), param.ELL);
    } else if (param.flag_recon == "False") {
	    sprintf(buf, "%s/pk%d.dat", param.output_dir.c_str(), param.ELL);
    }
	fp = fopen(buf, "w");
	fprintf(fp, "## k [h/Mpc]\t pk_real [(Mpc/h)^3]\t pk_imag [(Mpc/h)^3]\n");
	for(int i =0; i < param.n_kbin; i++) {
		fprintf(fp, "%.5f \t %.7e \t %.7e\n", kbin[i], norm * pk_save[i].real(), norm * pk_save[i].imag());
	}
	fclose(fp);

	delete [] pk_save;
	return 0;

}

int calcTwoPointCorrelationFunction(
        ParticleBOSSClass & P_D, ParticleBOSSClass & P_R, 
        ParameterClass & param, double alpha, double * rbin, double Vsurvey) {

    /*********************************************************************************************/
    /* This function calculates the two-point correlation function multipole, expanded by the Legendre function, using FFTs.*/
    /********************************************************************************************/

	if( !( (param.ELL == param.ell1) && (param.ell2 == 0) ) ) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
        return -1;
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleBOSSClass> dn_00(param);
	dn_00.calcYlmWeightedDensityFluctuation(P_D, P_R, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calcFourierTransform();

	/* define xi_save */
	std::complex<double> * xi_save = new std::complex<double>[param.n_rbin];
	for(int i =0; i < param.n_rbin; i++) { 
		xi_save[i] = 0.0; 
	}

	/* calc power spectrum */
	for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/****************************************/
		/* calc the ylm-weighted density fluctuation */
		DensityFieldClass<ParticleBOSSClass> dn_LM(param);
		dn_LM.calcYlmWeightedDensityFluctuation(P_D, P_R, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM.calcFourierTransform();

		/* calc shotnoise term */
		TwoPointStatisticsClass<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calcShotNoiseForPowerspectrum(P_D, P_R, alpha, param.ELL, _M_);
		
	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {

		/*****/
		/* (-1)**m1 delta_{m1, -M} */
		double w = wigner_3j(param.ell1, 0, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, 0, param.ELL, _m1_, 0, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}

		stat.calcCorrelationFunction(dn_LM, dn_00, rbin, shotnoise, param.ell1, _m1_);

		for(int i =0; i < param.n_rbin; i++) { 
			xi_save[i] += w * stat.xi[i];
	       	}

	}
		timesec = double(clock() - start);
		if(ThisTask == 0) {
			printf("M = %d | %.3f sec\n", _M_, timesec / CLOCKS_PER_SEC);
		}
	}

    /* calc the normalization factor */
	double norm = ParticleBOSSClass::calcNormalizationForPowerSpectrum(P_D, Vsurvey);
 
    /* save data */
	FILE * fp;
	char buf[1024];
    if(0) {
    } else if (param.flag_recon == "True") {
	    sprintf(buf, "%s/xi%d_recon.dat", param.output_dir.c_str(), param.ELL);
    } else if (param.flag_recon == "False") {
	    sprintf(buf, "%s/xi%d.dat", param.output_dir.c_str(), param.ELL);
    }
	fp = fopen(buf, "w");
	fprintf(fp, "## r [Mpc/h] \t xi\n");
	for(int i =0; i < param.n_rbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", rbin[i], norm * xi_save[i].real());
	}
	fclose(fp);

	delete [] xi_save;
	return 0;
}


int calcTwoPointWindowCorrelationFunction(
        ParticleBOSSClass & P_R, 
        ParameterClass & param, double * rbin, double Vsurvey) {

    /**************************************************************************************/
    /* This function calculates the two-point window correlation function multipole, expanded by the Legendre function, using FFTs.*/
    /**************************************************************************************/

	if( !( (param.ELL == param.ell1) && (param.ell2 == 0) ) ) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
        return -1;
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleBOSSClass> dn_00(param);
	dn_00.calcYlmWeightedDensity(P_R, 0, 0);
	/* Fourier transform*/
	dn_00.calcFourierTransform();

	/* define xi_save */
	std::complex<double> * xi_save = new std::complex<double>[param.n_rbin];

	for(int i =0; i < param.n_rbin; i++) { 
		xi_save[i] = 0.0; 
	}

	/* calc power spectrum */
	for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/****************************************/
		/* calc the ylm-weighted density fluctuation */
		DensityFieldClass<ParticleBOSSClass> dn_LM(param);
		dn_LM.calcYlmWeightedDensity(P_R, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM.calcFourierTransform();

		/* calc shotnoise term */
		TwoPointStatisticsClass<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calcShotNoiseForTwoPointWindowCorrelationFunction(P_R, param.ELL, _M_);
		
	    for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {

	    	/*****/
	    	/* (-1)**m1 delta_{m1, -M} */
	    	double w = wigner_3j(param.ell1, 0, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, 0, param.ELL, _m1_, 0, _M_);
	    	w *= double(2*param.ELL+1) * double(2*param.ell1+1);
	    	if(fabs(w) < 1.0e-10) {
	    		continue;
	    	}

	    	stat.calcCorrelationFunction(dn_LM, dn_00, rbin, shotnoise, param.ell1, _m1_);

	    	for(int i =0; i < param.n_rbin; i++) { 
	    		xi_save[i] += w * stat.xi[i];
	        }
	    }

		timesec = double(clock() - start);
		if(ThisTask == 0) {
			printf("M = %d | %.3f sec\n", _M_, timesec / CLOCKS_PER_SEC);
		}

	}

    /* calc the normalization factor */
	double norm = ParticleBOSSClass::calcNormalizationForPowerSpectrum(P_R, Vsurvey);

    /* save data */
	FILE * fp;
	char buf[1024];
    if(0) {
    } else if (param.flag_recon == "True") {
	    sprintf(buf, "%s/xi%d_recon_window.dat", param.output_dir.c_str(), param.ELL);
    } else if (param.flag_recon == "False") {
	    sprintf(buf, "%s/xi%d_window.dat", param.output_dir.c_str(), param.ELL);
    }
	fp = fopen(buf, "w");
	fprintf(fp, "## r [Mpc/h] \t xi\n");
	for(int i =0; i < param.n_rbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", rbin[i], norm * xi_save[i].real());
	}
	fclose(fp);

	delete [] xi_save;
	return 0;
}


int calcPowerSpectrumForBOX(
        ParticleGadgetClass & P_D, ParameterClass & param, double * kbin) {

    /***********************************************************************************************/
    /* This function calculates the power spectrum multipole, expanded by the Legendre function, using FFTs.*/
    /***********************************************************************************************/

	if( !( (param.ELL == param.ell1) && (param.ell2 == 0) ) ) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
        return -1;
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleGadgetClass> dn(param);
	dn.calcNormalDensityFluctuationForBOX(P_D, param);
	/* Fourier transform*/
	dn.calcFourierTransform();

	/* define pk_save */
	std::complex<double> * pk_save = new std::complex<double>[param.n_kbin];
	for(int i =0; i < param.n_kbin; i++) { 
		pk_save[i] = 0.0; 
	}

	/* calc shotnoise term */
	TwoPointStatisticsClass<ParticleGadgetClass> stat(param);
	std::complex<double> shotnoise = double(P_D.n_tot);
	
	stat.calcPowerSpectrum(dn, dn, kbin, shotnoise, param.ELL, 0);

	for(int i =0; i < param.n_kbin; i++) { 
		pk_save[i] += double(2*param.ELL+1) * stat.pk[i];
    }

	timesec = double(clock() - start);
	if(ThisTask == 0) {
		printf("%.3f sec\n", timesec / CLOCKS_PER_SEC);
	}

    /* calc the normalization factor */
	double norm = param.volume / double(P_D.n_tot) / double(P_D.n_tot);

    /* save */
	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/pk%d.dat", param.output_dir.c_str(), param.ELL);
	fp = fopen(buf, "w");
	fprintf(fp, "## k [h/Mpc]\t pk_real [(Mpc/h)^3]\t pk_imag [(Mpc/h)^3]\n");
	for(int i =0; i < param.n_kbin; i++) {
		fprintf(fp, "%.5f \t %.7e \t %.7e\n", kbin[i], norm * pk_save[i].real(), norm * pk_save[i].imag());
	}
	fclose(fp);

	delete [] pk_save;

	return 0;
}

int calcTwoPointCorrelationFunctionForBOX(
        ParticleGadgetClass & P_D, ParameterClass & param, double * rbin) {

    /************************************************************************************/
    /* This function calculates the two-point correlation function multipole, expanded by the Legendre function, using FFTs.*/
    /************************************************************************************/

	if( !( (param.ELL == param.ell1) && (param.ell2 == 0) ) ) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
        return -1;
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleGadgetClass> dn(param);
	dn.calcNormalDensityFluctuationForBOX(P_D, param);
	/* Fourier transform*/
	dn.calcFourierTransform();

	/* define xi_save */
	std::complex<double> * xi_save = new std::complex<double>[param.n_rbin];
	for(int i =0; i < param.n_rbin; i++) { 
		xi_save[i] = 0.0; 
	}

	/* calc power spectrum */
	TwoPointStatisticsClass<ParticleGadgetClass> stat(param);
	std::complex<double> shotnoise = double(P_D.n_tot);
	
	stat.calcCorrelationFunction(dn, dn, rbin, shotnoise, param.ELL, 0);
	
	for(int i =0; i < param.n_rbin; i++) { 
		xi_save[i] += double(2*param.ELL+1) * stat.xi[i];
	}
	
	timesec = double(clock() - start);
	if(ThisTask == 0) {
		printf("%.3f sec\n", timesec / CLOCKS_PER_SEC);
	}

    /* calc. the normalization factor */
	double norm = param.volume/double(P_D.n_tot)/double(P_D.n_tot);

    /* save */
	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/xi%d.dat", param.output_dir.c_str(), param.ELL);
	fp = fopen(buf, "w");
	fprintf(fp, "## r [Mpc/h] \t xi\n");
	for(int i =0; i < param.n_rbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", rbin[i], norm * xi_save[i].real());
	}
	fclose(fp);

	delete [] xi_save;

	return 0;
}


int calcPowerSpectrumForBOXForReconstruction(
        ParticleGadgetClass & P_D, ParticleGadgetClass & P_R, ParameterClass & param, double alpha, double * kbin) {

    /***************************************************************************************/
    /* This function calculates the reconstructed power spectrum multipole, expanded by the Legendre function, using FFTs.*/
    /***************************************************************************************/

	if( !( (param.ELL == param.ell1) && (param.ell2 == 0) ) ) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
        return -1;
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleGadgetClass> dn(param);
	dn.calcNormalDensityFluctuationForBOXForReconstruction(P_D, P_R, alpha);
	/* Fourier transform*/
	dn.calcFourierTransform();

	/* define pk_save */
	std::complex<double> * pk_save = new std::complex<double>[param.n_kbin];
	for(int i =0; i < param.n_kbin; i++) { 
		pk_save[i] = 0.0; 
	}

	/* calc shotnoise term */
	TwoPointStatisticsClass<ParticleGadgetClass> stat(param);
	std::complex<double> shotnoise = stat.calcShotNoiseForPowerspectrumForBOXForReconstruction(P_D, P_R, alpha);
	
	stat.calcPowerSpectrum(dn, dn, kbin, shotnoise, param.ELL, 0);

	for(int i =0; i < param.n_kbin; i++) { 
		pk_save[i] += double(2*param.ELL+1) * stat.pk[i];
    }

	timesec = double(clock() - start);
	if(ThisTask == 0) {
		printf("%.3f sec\n", timesec / CLOCKS_PER_SEC);
	}

    /* calc. the normalization factor */
	double norm = param.volume / double(P_D.n_tot) / double(P_D.n_tot);

    /* save */
	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/pk%d_recon.dat", param.output_dir.c_str(), param.ELL);
	fp = fopen(buf, "w");
	fprintf(fp, "## k [h/Mpc]\t pk_real [(Mpc/h)^3]\t pk_imag [(Mpc/h)^3]\n");
	for(int i =0; i < param.n_kbin; i++) {
		fprintf(fp, "%.5f \t %.7e \t %.7e\n", kbin[i], norm * pk_save[i].real(), norm * pk_save[i].imag());
	}
	fclose(fp);

	delete [] pk_save;

	return 0;

}


int calcTwoPointFunctionForBOXForReconstruction(
        ParticleGadgetClass & P_D, ParticleGadgetClass & P_R,
        ParameterClass & param, double alpha, double * rbin) {

    /*******************************************************************************************/
    /* This function calculates the two-point correlation function multipole, expanded by the Legendre function, using FFTs.*/
    /*******************************************************************************************/

	if( !( (param.ELL == param.ell1) && (param.ell2 == 0) ) ) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
        return -1;
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleGadgetClass> dn(param);
	dn.calcNormalDensityFluctuationForBOXForReconstruction(P_D, P_R, alpha);
	/* Fourier transform*/
	dn.calcFourierTransform();

	/* define xi_save */
	std::complex<double> * xi_save = new std::complex<double>[param.n_rbin];
	for(int i =0; i < param.n_rbin; i++) { 
		xi_save[i] = 0.0; 
	}

	/* calc power spectrum */
	TwoPointStatisticsClass<ParticleGadgetClass> stat(param);
	std::complex<double> shotnoise = stat.calcShotNoiseForPowerspectrumForBOXForReconstruction(P_D, P_R, alpha);
	
	stat.calcCorrelationFunction(dn, dn, rbin, shotnoise, param.ELL, 0);
	
	for(int i =0; i < param.n_rbin; i++) { 
		xi_save[i] += double(2*param.ELL+1) * stat.xi[i];
	}
	
	timesec = double(clock() - start);
	if(ThisTask == 0) {
		printf("%.3f sec\n", timesec / CLOCKS_PER_SEC);
	}

    /* calc. the normalization factor */
	double norm = param.volume / double(P_D.n_tot) / double(P_D.n_tot);

    /* save */
	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/xi%d_recon.dat", param.output_dir.c_str(), param.ELL);
	fp = fopen(buf, "w");
	fprintf(fp, "## r [Mpc/h] \t xi\n");
	for(int i =0; i < param.n_rbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", rbin[i], norm * xi_save[i].real());
	}
	fclose(fp);

	delete [] xi_save;

	return 0;
}


#endif

