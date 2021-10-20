#ifndef __bispec__
#define __bispec__

int calcBiSpectrum(
        ParticleBOSSClass & P_D, ParticleBOSSClass & P_R, 
        ParameterClass & param, double alpha, double * kbin, double Vsurvey) {

    /**************************************************************************************************/
    /* This function calculates the bispectrum multipole, expanded by the TripoSH function, using FFTs.*/
    /**************************************************************************************************/

	if(fabs(wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0)) < 1.0e-10) {
		if(ThisTask == 0) { 
			printf("This multipole combination is not allowed.\n");
			printf("It should be wigner_3j(ell1, ell2, ELL, 0,0,0) != 0\n");
		}
        return -1;
	}

	/*****************************************************************************/
	/* shotnoise */
	/*****************************************************************************/
	/* time */
	timesec = double(clock() - start);
	if(ThisTask == 0) { printf("Computing shotnoise terms.| %.3f sec\n", timesec / CLOCKS_PER_SEC);}

	/* define sn_save:: shotnoise term */
	std::complex<double> * sn_save = new std::complex<double>[param.n_kbin];
	for(int i = 0; i < param.n_kbin; i++) {
		sn_save[i] = 0.0;
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleBOSSClass> dn_00_sn(param);
	dn_00_sn.calcYlmWeightedDensityFluctuation(P_D, P_R, alpha, 0, 0);
	/* Fourier transform*/
	dn_00_sn.calcFourierTransform();

	/* calc shot noise terms */
	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {
	for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
	
		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}
		
		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityFieldClass<ParticleBOSSClass> N_LM_sn(param);
		N_LM_sn.calcYlmWeightedDensityFluctuationForBispectrumShotnoise(P_D, P_R, alpha, param.ELL, _M_);
		/* Fourier transform*/
		N_LM_sn.calcFourierTransform();

		TwoPointStatisticsClass<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calcShotNoiseForBispectrum(P_D, P_R, alpha, param.ELL, _M_);

		if( (param.ell1 == 0) && (param.ell2 == 0) ) {
			for(int i =0; i < param.n_kbin; i++) { 
				sn_save[i] += w * shotnoise;
	       	}
		}

		if( (param.ell2 == 0) ) {
			stat.calcPowerSpectrum(dn_00_sn, N_LM_sn, kbin, shotnoise, param.ell1, _m1_);

                for(int i =0; i < param.n_kbin; i++) { 
                	sn_save[i] += w * stat.pk[param.ith_kbin];
                }
		}

		if( (param.ell1 == 0) ) {
			stat.calcPowerSpectrum(dn_00_sn, N_LM_sn, kbin, shotnoise, param.ell2, _m2_);
			for(int i =0; i < param.n_kbin; i++) { 
				sn_save[i] += w * stat.pk[i];
		    }
		}

		/* time */
		timesec = double(clock() - start);
		if(ThisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}

	}}}
	
	/****/
	dn_00_sn.finalizeDensityField();
	/****************************************/

	/****************************************/
	/* calc F00 */
	DensityFieldClass<ParticleBOSSClass> N_00_sn(param);
	N_00_sn.calcYlmWeightedDensityFluctuationForBispectrumShotnoise(P_D, P_R, alpha, 0, 0);
	/* Fourier transform*/
	N_00_sn.calcFourierTransform();

	/* store spherical bessel functions */
	SphericalBesselClass sj1(param.ell1);
	SphericalBesselClass sj2(param.ell2);

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.n_mesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.n_mesh_tot];
		byte += double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
        if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }
	
		std::string flag = "FALSE";
		for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if(fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if(flag == "TRUE") {
			ToolClass::storeReducedSphericalHarmonicsInConfigurationSpace( param.ell1, _m1_, param, Ylm1 );
			ToolClass::storeReducedSphericalHarmonicsInConfigurationSpace( param.ell2, _m2_, param, Ylm2 );
		}

	for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
	
		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}
		
		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityFieldClass<ParticleBOSSClass> dn_LM_sn(param);
		dn_LM_sn.calcYlmWeightedDensityFluctuation(P_D, P_R, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM_sn.calcFourierTransform();

		TwoPointStatisticsClass<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calcShotNoiseForBispectrum(P_D, P_R, alpha, param.ELL, _M_);

		fftw_complex * xi = fftw_alloc_complex(param.n_mesh_tot);
		byte += double( sizeof(fftw_complex) * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }
        for(long long i = 0; i < param.n_mesh_tot; i++) {
			xi[i][0] = 0.0;
			xi[i][1] = 0.0;
		}

		stat.calcShotNoiseForBispectrum_ijk(dn_LM_sn, N_00_sn, shotnoise, param.ELL, _M_, xi);

		for(int ik = 0; ik < param.n_kbin; ik++) {
			
			double kmag2 = kbin[ik];
			double kmag1 = 0.0;
			kmag1 = kbin[param.ith_kbin];

			/* calc shotnoise */
			std::complex<double> sn_sum = 0.0;
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
				double rmag = sqrt( rvec[0] * rvec[0] +  rvec[1] * rvec[1] +  rvec[2] * rvec[2]);
				std::complex<double> ff(xi[coord][0], xi[coord][1]);
				double j1 = sj1.getSphericalBessel(kmag1 * rmag);
				double j2 = sj2.getSphericalBessel(kmag2 * rmag);

				sn_sum += ( j1 * j2 * ff * Ylm1[coord] * Ylm2[coord] );

			}}}

			std::complex<double> _I_(0.0,1.0);
			double fac = param.volume / double(param.n_mesh_tot);
			sn_sum *= pow(_I_, param.ell1 + param.ell2) * fac;

			sn_save[ik] += ( w * sn_sum );

			/* time */
			timesec = double(clock() - start);
			if(ThisTask == 0) { printf("k2 = %.3f [h/Mpc], m1 = %d, m2 = %d, M = %d| %.3f sec\n", kmag2, _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}
		}

		fftw_free(xi); xi = NULL;
		byte -= double( sizeof(fftw_complex) * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

	}

		delete [] Ylm1; Ylm1 = NULL;
		delete [] Ylm2; Ylm2 = NULL;
		byte -= double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

	}}
	
	/****/
	N_00_sn.finalizeDensityField();
	/****************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* bispectrum */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	if(ThisTask == 0) { printf("Computing the bispectrum multipole.\n");}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleBOSSClass> dn_00(param);
	dn_00.calcYlmWeightedDensityFluctuation(P_D, P_R, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calcFourierTransform();


	/* define bk_save:: bispectrum */
	std::complex<double> * bk_save = new std::complex<double>[param.n_kbin];
	for(int i = 0; i < param.n_kbin; i++) {
		bk_save[i] = 0.0;
	}

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		/****************************************************/
		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.n_mesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.n_mesh_tot];
		byte += double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

		std::string flag = "FALSE";
		for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if(fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if(flag == "TRUE") {
		
			ToolClass::storeReducedSphericalHarmonicsInFourierSpace(param.ell1, _m1_, param, Ylm1);
			ToolClass::storeReducedSphericalHarmonicsInFourierSpace(param.ell2, _m2_, param, Ylm2);
		}

	for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
	
		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}
	
		/* calc the yLM-weighted density perturbation */
		DensityFieldClass<ParticleBOSSClass> dn_LM(param);
		dn_LM.calcYlmWeightedDensityFluctuation(P_D, P_R, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM.calcFourierTransform();
		/* divided by a assignment function */
		dn_LM.calcAssignmentFunctionCorrection();
		/* Inverse Fourier transform*/
		dn_LM.calcInverseFourierTransform();

		/* calc dn_tilde1 */
		DensityFieldClass<ParticleBOSSClass> dn_tilde1(param);
		double kmag1;
		double dk = kbin[1] - kbin[0];
		kmag1 = kbin[param.ith_kbin];
		dn_tilde1.calcInverseFourierTransformForBispectrum(dn_00, kmag1, dk, Ylm1);

		for(int ik = 0; ik < param.n_kbin; ik++) {

			double kmag2 = kbin[ik];

			/* calc dn_tilde2 */
			DensityFieldClass<ParticleBOSSClass> dn_tilde2(param);
			dn_tilde2.calcInverseFourierTransformForBispectrum(dn_00, kmag2, dk, Ylm2);
			
			/* calc bispectrum */
			std::complex<double> bk_sum = 0.0;
			double fac = param.volume / double(param.n_mesh_tot);
			for(long long coord = 0; coord < param.n_mesh_tot; coord++ ) {
				std::complex<double> f1(dn_tilde1[coord][0], dn_tilde1[coord][1]);
				std::complex<double> f2(dn_tilde2[coord][0], dn_tilde2[coord][1]);
				std::complex<double> f3(dn_LM[coord][0], dn_LM[coord][1]);
				bk_sum += (fac * f1 * f2 * f3);
			}

			bk_save[ik] += (w * bk_sum);

			double timesec = double(clock() - start);
			if(ThisTask == 0) { printf("k2 = %.3f [h/Mpc], m1 = %d, m2 = %d, M = %d| %.3f sec\n", kmag2, _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}

		}
	
	}
	
		delete [] Ylm1; Ylm1 = NULL;
		delete [] Ylm2; Ylm2 = NULL;
		byte -= double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }


	}}

    /* the normalization factor */
	double norm = ParticleBOSSClass::calcNormalizationForBiSpectrum(P_D, Vsurvey);

    /* save data */
	FILE * fp;
	char buf[1024];
    if(0) {
    } else if (param.flag_recon == "True") {
        sprintf(buf, "%s/bk%d%d%d_%02d_recon.dat", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_kbin);
    } else if (param.flag_recon == "False") {
        sprintf(buf, "%s/bk%d%d%d_%02d.dat", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_kbin);
    }
    fp = fopen(buf, "w");
    fprintf(fp, "## k1 [h/Mpc] \t k2 [h/Mpc] \t bk_real \t bk_imag\n");
    for(int i = 0; i < param.n_kbin; i++) {
    	fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n", kbin[param.ith_kbin], kbin[i], 
    			norm * ( bk_save[i].real() - sn_save[i].real() ), norm * ( bk_save[i].imag() - sn_save[i].imag() ) );
    }
    fclose(fp);

	delete [] sn_save;
	delete [] bk_save;

	return 0;
}

int calcThreePointCorrelationFunction(
        ParticleBOSSClass & P_D, ParticleBOSSClass & P_R, 
        ParameterClass & param, double alpha, double * rbin, double Vsurvey) {

    /********************************************************************************************/
    /* This function calculates the 3pcf multipole, expanded by the TripoSH function, using FFTs.*/
    /********************************************************************************************/

	if(fabs(wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0)) < 1.0e-10) {
		if(ThisTask == 0) { 
			printf("This multipole combination is not allowed.\n");
			printf("It should be wigner_3j(ell1, ell2, ELL, 0,0,0) != 0\n");
		}
        return -1;
	}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* shotnoise */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	timesec = double(clock() - start);
	if(ThisTask == 0) { printf("Computing shotnoise terms.| %.3f sec\n", timesec / CLOCKS_PER_SEC);}

	/* define sn_save:: shotnoise term */
	std::complex<double> * sn_save = new std::complex<double>[param.n_rbin];
	for(int i = 0; i < param.n_rbin; i++) {
		sn_save[i] = 0.0;
	}

	/****************************************/
	/* calc F00 */
	DensityFieldClass<ParticleBOSSClass> N_00_sn(param);
	N_00_sn.calcYlmWeightedDensityFluctuationForBispectrumShotnoise(P_D, P_R, alpha, 0, 0);
	/* Fourier transform*/
	N_00_sn.calcFourierTransform();

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.n_mesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.n_mesh_tot];
		byte += double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

		std::string flag = "FALSE";
		for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if(fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if(flag == "TRUE") {

			ToolClass::storeReducedSphericalHarmonicsInConfigurationSpace( param.ell1, _m1_, param, Ylm1 );
			ToolClass::storeReducedSphericalHarmonicsInConfigurationSpace( param.ell2, _m2_, param, Ylm2 );
		
		}

	for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
	
		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}
		
		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityFieldClass<ParticleBOSSClass> dn_LM_sn(param);
		dn_LM_sn.calcYlmWeightedDensityFluctuation(P_D, P_R, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM_sn.calcFourierTransform();

		TwoPointStatisticsClass<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calcShotNoiseForBispectrum(P_D, P_R, alpha, param.ELL, _M_);

		stat.calcCorrelationFunctionForThreePointFunction(dn_LM_sn, N_00_sn, rbin, shotnoise, param.ell1, _m1_, Ylm1, Ylm2);

		for(int i =0; i < param.n_rbin; i++) { 
			if(i == param.ith_rbin) {
				sn_save[i] += w * stat.xi[i];
			} else {
				sn_save[i] += 0.0;
			}
	    }

		/* time */
		timesec = double(clock() - start);
		if(ThisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}

	}

		delete [] Ylm1; Ylm1 = NULL;
		delete [] Ylm2; Ylm2 = NULL;
		byte -= double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }
	}}
	
	/****/
	N_00_sn.finalizeDensityField();
	/****************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* three-point function */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	if(ThisTask == 0) { printf("Computing the three-point function multipole.\n");}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleBOSSClass> dn_00(param);
	dn_00.calcYlmWeightedDensityFluctuation(P_D, P_R, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calcFourierTransform();

	/* store spherical bessel functions */
	SphericalBesselClass sj1(param.ell1);
	SphericalBesselClass sj2(param.ell2);


	/* define zeta_save:: bispectrum */
	std::complex<double> * zeta_save = new std::complex<double>[param.n_rbin];
	for(int i = 0; i < param.n_rbin; i++) {
		zeta_save[i] = 0.0;
	}

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {
	
		/****************************************************/
		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.n_mesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.n_mesh_tot];
		byte += double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

		std::string flag = "FALSE";
		for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if(fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if(flag == "TRUE") {

			ToolClass::storeReducedSphericalHarmonicsInFourierSpace(param.ell1, _m1_, param, Ylm1);
			ToolClass::storeReducedSphericalHarmonicsInFourierSpace(param.ell2, _m2_, param, Ylm2);
	
		}

	for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
	
		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}
	
		/* calc the yLM-weighted density perturbation */
		DensityFieldClass<ParticleBOSSClass> dn_LM(param);
		dn_LM.calcYlmWeightedDensityFluctuation(P_D, P_R, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM.calcFourierTransform();
		/* divided by a assignment function */
		dn_LM.calcAssignmentFunctionCorrection();
		/* Inverse Fourier transform*/
		dn_LM.calcInverseFourierTransform();

		/* calc dn_tilde1 */
		DensityFieldClass<ParticleBOSSClass> dn_tilde1(param);
		double rmag1;
	    rmag1 = rbin[param.ith_rbin];
	    dn_tilde1.calcInverseFourierTransformForThreePointFunction(dn_00, rmag1, Ylm1, sj1);

		for(int ir = 0; ir < param.n_rbin; ir++) {
	
			double rmag2 = rbin[ir];

			/* calc dn_tilde2 */
			DensityFieldClass<ParticleBOSSClass> dn_tilde2(param);
			dn_tilde2.calcInverseFourierTransformForThreePointFunction(dn_00, rmag2, Ylm2, sj2);
			
			/* calc bispectrum */
			std::complex<double> zeta_sum = 0.0;
			double fac = param.volume / double(param.n_mesh_tot);
			std::complex<double> _I_(0.0, 1.0);
			for(long long coord = 0; coord < param.n_mesh_tot; coord++ ) {
				std::complex<double> f1(dn_tilde1[coord][0], dn_tilde1[coord][1]);
				std::complex<double> f2(dn_tilde2[coord][0], dn_tilde2[coord][1]);
				std::complex<double> f3(dn_LM[coord][0], dn_LM[coord][1]);
				zeta_sum += (pow(_I_, param.ell1+param.ell2) * fac * f1 * f2 * f3);
			}

			zeta_save[ir] += (w * zeta_sum);

			double timesec = double(clock() - start);
			if(ThisTask == 0) { printf("r2 = %.3f [Mpc/h], m1 = %d, m2 = %d, M = %d| %.3f sec\n", rmag2, _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}

		}
	
	}
	
		delete [] Ylm1; Ylm1 = NULL;
		delete [] Ylm2; Ylm2 = NULL;
		byte -= double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

	}}

    /* the normalization factor */
	double norm = ParticleBOSSClass::calcNormalizationForBiSpectrum(P_D, Vsurvey);

    /* save data */
	FILE * fp;
	char buf[1024];
    if(0) {
    } else if (param.flag_recon == "True") {
        sprintf(buf, "%s/zeta%d%d%d_%02d_recon.dat", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_rbin);
    } else if (param.flag_recon == "False") {
        sprintf(buf, "%s/zeta%d%d%d_%02d.dat", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_rbin);
    }
    fp = fopen(buf, "w");
    fprintf(fp, "## r1 [Mpc/h] \t r2 [Mpc/h] \t zeta \n");
    for(int i = 0; i < param.n_rbin; i++) {
    	fprintf(fp, "%.5f \t %.5f \t %.7e\n",  rbin[param.ith_rbin], rbin[i], norm * (zeta_save[i].real() - sn_save[i].real()));
    }
    fclose(fp);

	delete [] sn_save;
	delete [] zeta_save;

	return 0;
}


int calcThreePointWindowCorrelationFunction(
        ParticleBOSSClass & P_R, ParameterClass & param, double * rbin, double Vsurvey) {

    /***************************************************************************************************/
    /* This function calculates the window 3pcf multipole, expanded by the TripoSH function, using FFTs.*/
    /***************************************************************************************************/

	if(fabs(wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0)) < 1.0e-10) {
		if(ThisTask == 0) { 
			printf("This multipole combination is not allowed.\n");
			printf("It should be wigner_3j(ell1, ell2, ELL, 0,0,0) != 0\n");
		}
        return -1;
	}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* shotnoise */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	timesec = double(clock() - start);
	if(ThisTask == 0) { printf("Computing shotnoise terms.| %.3f sec\n", timesec / CLOCKS_PER_SEC);}

	/* define sn_save:: shotnoise term */
	std::complex<double> * sn_save = new std::complex<double>[param.n_rbin];
	for(int i = 0; i < param.n_rbin; i++) {
		sn_save[i] = 0.0;
	}

	/****************************************/
	/* calc F00 */
	DensityFieldClass<ParticleBOSSClass> N_00_sn(param);
	N_00_sn.calcYlmWeightedDensityForThreePointWindowFunctionShotnoise(P_R, 0, 0);

	/* Fourier transform*/
	N_00_sn.calcFourierTransform();

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.n_mesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.n_mesh_tot];
		byte += double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }
		std::string flag = "FALSE";
		for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if(fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if(flag == "TRUE") {

			ToolClass::storeReducedSphericalHarmonicsInConfigurationSpace( param.ell1, _m1_, param, Ylm1 );
			ToolClass::storeReducedSphericalHarmonicsInConfigurationSpace( param.ell2, _m2_, param, Ylm2 );
		
		}

	for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
	
		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}
		
		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityFieldClass<ParticleBOSSClass> dn_LM_sn(param);
		dn_LM_sn.calcYlmWeightedDensity(P_R, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM_sn.calcFourierTransform();

		TwoPointStatisticsClass<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calcShotNoiseForTwoPointWindowCorrelationFunction(P_R, param.ELL, _M_);

		stat.calcCorrelationFunctionForThreePointFunction(dn_LM_sn, N_00_sn, rbin, shotnoise, param.ell1, _m1_, Ylm1, Ylm2);

		for(int i =0; i < param.n_rbin; i++) { 
			if(i == param.ith_rbin) {
				sn_save[i] += w * stat.xi[i];
			} else {
				sn_save[i] += 0.0;
			}
	    }

		/* time */
		timesec = double(clock() - start);
		if(ThisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}

	}

		delete [] Ylm1; Ylm1 = NULL;
		delete [] Ylm2; Ylm2 = NULL;
		byte -= double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }
	}}
	
	/****/
	N_00_sn.finalizeDensityField();
	/****************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* three-point function */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	if(ThisTask == 0) { printf("Computing the three-point function multipole.\n");}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleBOSSClass> dn_00(param);
	dn_00.calcYlmWeightedDensity(P_R, 0, 0);
	/* Fourier transform*/
	dn_00.calcFourierTransform();

	/* store spherical bessel functions */
	SphericalBesselClass sj1(param.ell1);
	SphericalBesselClass sj2(param.ell2);


	/* define zeta_save:: bispectrum */
	std::complex<double> * zeta_save = new std::complex<double>[param.n_rbin];
	for(int i = 0; i < param.n_rbin; i++) {
		zeta_save[i] = 0.0;
	}

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {
	
		/****************************************************/
		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.n_mesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.n_mesh_tot];
		byte += double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }
		std::string flag = "FALSE";
		for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if(fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if(flag == "TRUE") {

			ToolClass::storeReducedSphericalHarmonicsInFourierSpace(param.ell1, _m1_, param, Ylm1);
			ToolClass::storeReducedSphericalHarmonicsInFourierSpace(param.ell2, _m2_, param, Ylm2);
	
		}

	for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
	
		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}
	
		/* calc the yLM-weighted density perturbation */
		DensityFieldClass<ParticleBOSSClass> dn_LM(param);
		dn_LM.calcYlmWeightedDensity(P_R, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM.calcFourierTransform();
		/* divided by a assignment function */
		dn_LM.calcAssignmentFunctionCorrection();
		/* Inverse Fourier transform*/
		dn_LM.calcInverseFourierTransform();

		/* calc dn_tilde1 */
		DensityFieldClass<ParticleBOSSClass> dn_tilde1(param);
		double rmag1;
		rmag1 = rbin[param.ith_rbin];
		dn_tilde1.calcInverseFourierTransformForThreePointFunction(dn_00, rmag1, Ylm1, sj1);

		for(int ir = 0; ir < param.n_rbin; ir++) {
	
			double rmag2 = rbin[ir];

			/* calc dn_tilde2 */
			DensityFieldClass<ParticleBOSSClass> dn_tilde2(param);
			dn_tilde2.calcInverseFourierTransformForThreePointFunction(dn_00, rmag2, Ylm2, sj2);
			
			/* calc bispectrum */
			std::complex<double> zeta_sum = 0.0;
			double fac = param.volume / double(param.n_mesh_tot);
			std::complex<double> _I_(0.0, 1.0);
			for(long long coord = 0; coord < param.n_mesh_tot; coord++ ) {
				std::complex<double> f1(dn_tilde1[coord][0], dn_tilde1[coord][1]);
				std::complex<double> f2(dn_tilde2[coord][0], dn_tilde2[coord][1]);
				std::complex<double> f3(dn_LM[coord][0], dn_LM[coord][1]);
				zeta_sum += (pow(_I_, param.ell1+param.ell2) * fac * f1 * f2 * f3);
			}

			zeta_save[ir] += (w * zeta_sum);

			double timesec = double(clock() - start);
			if(ThisTask == 0) { printf("r2 = %.3f [Mpc/h], m1 = %d, m2 = %d, M = %d| %.3f sec\n", rmag2, _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}

		}
	
	}
	
		delete [] Ylm1; Ylm1 = NULL;
		delete [] Ylm2; Ylm2 = NULL;
		byte -= double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

	}}

    /* the normalization factor */
	double norm = ParticleBOSSClass::calcNormalizationForBiSpectrum(P_R, Vsurvey);

    /* save data */
	FILE * fp;
	char buf[1024];
    if(0) {
    } else if (param.flag_recon == "True") {
        sprintf(buf, "%s/zeta%d%d%d_%02d_recon_window.dat", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_rbin);
    } else if (param.flag_recon == "False") {
        sprintf(buf, "%s/zeta%d%d%d_%02d_window.dat", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_rbin);
    }
    fp = fopen(buf, "w");
    fprintf(fp, "## r1 [Mpc/h] \t r2 [Mpc/h] \t zeta \n");
    for(int i = 0; i < param.n_rbin; i++) {
    	fprintf(fp, "%.5f \t %.5f \t %.7e\n",  rbin[param.ith_rbin], rbin[i], norm * (zeta_save[i].real() - sn_save[i].real()));
    }
    fclose(fp);

	delete [] sn_save;
	delete [] zeta_save;

	return 0;
}


int calcBiSpectrumForBOX(
        ParticleGadgetClass & P_D, ParameterClass & param, double * kbin) {

    /**************************************************************************************************/
    /* This function calculates the bispectrum multipole, expanded by the TripoSH function, using FFTs.*/
    /**************************************************************************************************/

	if(fabs(wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0)) < 1.0e-10) {
		if(ThisTask == 0) { 
			printf("This multipole combination is not allowed.\n");
			printf("It should be wigner_3j(ell1, ell2, ELL, 0,0,0) != 0\n");
		}
        return -1;
	}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* shotnoise */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	timesec = double(clock() - start);
	if(ThisTask == 0) { printf("Computing shotnoise terms.| %.3f sec\n", timesec / CLOCKS_PER_SEC);}

	/* define sn_save:: shotnoise term */
	std::complex<double> * sn_save = new std::complex<double>[param.n_kbin];
	for(int i = 0; i < param.n_kbin; i++) {
		sn_save[i] = 0.0;
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleGadgetClass> dn_00_sn(param);
	dn_00_sn.calcNormalDensityFluctuationForBOX(P_D, param);
	/* Fourier transform*/
	dn_00_sn.calcFourierTransform();

    /**********************************************/
    /* calc the yLM-weighted density perturbation */
    DensityFieldClass<ParticleGadgetClass> N_LM_sn(param);
    N_LM_sn.calcNormalDensityForBispectrumShotnoiseForBOX(P_D);
    /* Fourier transform*/
    N_LM_sn.calcFourierTransform();

	/* calc shot noise terms */
	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {
		int _M_ = 0;
	
		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}
		
		TwoPointStatisticsClass<ParticleGadgetClass> stat(param);
		std::complex<double> shotnoise = double(P_D.n_tot);

		if( (param.ell1 == 0) && (param.ell2 == 0) ) {
			for(int i =0; i < param.n_kbin; i++) { 
				sn_save[i] += w * shotnoise;
	       		}
		}

		if( (param.ell2 == 0) ) {
		    stat.calcPowerSpectrum(dn_00_sn, N_LM_sn, kbin, shotnoise, param.ell1, _m1_);
			for(int i =0; i < param.n_kbin; i++) { 
				sn_save[i] += w * stat.pk[param.ith_kbin];
			}
		}

		if( (param.ell1 == 0) ) {
			stat.calcPowerSpectrum(dn_00_sn, N_LM_sn, kbin, shotnoise, param.ell2, _m2_);
			for(int i =0; i < param.n_kbin; i++) { 
				sn_save[i] += w * stat.pk[i];
		       	}
		}

		/* time */
		timesec = double(clock() - start);
		if(ThisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}

	}}
	
	/****/
	dn_00_sn.finalizeDensityField();
	N_LM_sn.finalizeDensityField();
	/****************************************/

	/****************************************/
	/* calc F00 */
	DensityFieldClass<ParticleGadgetClass> N_00_sn(param);
	N_00_sn.calcNormalDensityForBispectrumShotnoiseForBOX(P_D);
	/* Fourier transform*/
	N_00_sn.calcFourierTransform();

    /**********************************************/
    /* calc the yLM-weighted density perturbation */
    DensityFieldClass<ParticleGadgetClass> dn_LM_sn(param);
    dn_LM_sn.calcNormalDensityFluctuationForBOX(P_D, param);
    /* Fourier transform*/
    dn_LM_sn.calcFourierTransform();

	/* store spherical bessel functions */
	SphericalBesselClass sj1(param.ell1);
	SphericalBesselClass sj2(param.ell2);

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		int _M_ = 0;

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}

		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.n_mesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.n_mesh_tot];
		byte += double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

		ToolClass::storeReducedSphericalHarmonicsInConfigurationSpace( param.ell1, _m1_, param, Ylm1 );
		ToolClass::storeReducedSphericalHarmonicsInConfigurationSpace( param.ell2, _m2_, param, Ylm2 );

		TwoPointStatisticsClass<ParticleGadgetClass> stat(param);
		std::complex<double> shotnoise = double(P_D.n_tot);

		fftw_complex * xi = fftw_alloc_complex(param.n_mesh_tot);
		byte += double( sizeof(fftw_complex) * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

		for(long long i = 0; i < (param.n_mesh_tot); i++) {
			xi[i][0] = 0.0;
			xi[i][1] = 0.0;
		}

		stat.calcShotNoiseForBispectrum_ijk(dn_LM_sn, N_00_sn, shotnoise, param.ELL, _M_, xi);

		for(int ik = 0; ik < param.n_kbin; ik++) {
			
			double kmag2 = kbin[ik];
			double kmag1 = 0.0;
			kmag1 = kbin[param.ith_kbin];

			/* calc shotnoise */
			std::complex<double> sn_sum = 0.0;
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
				double rmag = sqrt( rvec[0] * rvec[0] +  rvec[1] * rvec[1] +  rvec[2] * rvec[2]);
				std::complex<double> ff(xi[coord][0], xi[coord][1]);
				double j1 = sj1.getSphericalBessel(kmag1 * rmag);
				double j2 = sj2.getSphericalBessel(kmag2 * rmag);

				sn_sum += ( j1 * j2 * ff * Ylm1[coord] * Ylm2[coord] );

			}}}

			std::complex<double> _I_(0.0,1.0);
			double fac = param.volume / double(param.n_mesh_tot);
			sn_sum *= pow(_I_, param.ell1 + param.ell2) * fac;

			sn_save[ik] += ( w * sn_sum );

			/* time */
			timesec = double(clock() - start);
			if(ThisTask == 0) { printf("k2 = %.3f [h/Mpc], m1 = %d, m2 = %d, M = %d| %.3f sec\n", kmag2, _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}
		}

		fftw_free(xi); xi = NULL;
		byte -= double( sizeof(fftw_complex) * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }


		delete [] Ylm1; Ylm1 = NULL;
		delete [] Ylm2; Ylm2 = NULL;
		byte -= double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

	}}
	
	/****/
	N_00_sn.finalizeDensityField();
	dn_LM_sn.finalizeDensityField();
	/****************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* bispectrum */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	if(ThisTask == 0) { printf("computing bispectrum...\n");}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleGadgetClass> dn_00(param);
	dn_00.calcNormalDensityFluctuationForBOX(P_D, param);
	/* Fourier transform*/
	dn_00.calcFourierTransform();

	/* calc the yLM-weighted density perturbation */
	DensityFieldClass<ParticleGadgetClass> dn_LM(param);
	dn_LM.calcNormalDensityFluctuationForBOX(P_D, param);
	/* Fourier transform*/
	dn_LM.calcFourierTransform();
	/* divided by a assignment function */
	dn_LM.calcAssignmentFunctionCorrection();
	/* Inverse Fourier transform*/
	dn_LM.calcInverseFourierTransform();

	/* define bk_save:: bispectrum */
	std::complex<double> * bk_save = new std::complex<double>[param.n_kbin];
	for(int i = 0; i < param.n_kbin; i++) {
		bk_save[i] = 0.0;
	}

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		int _M_ = 0;


		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}

		/****************************************************/
		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.n_mesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.n_mesh_tot];
		byte += double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }
	
		ToolClass::storeReducedSphericalHarmonicsInFourierSpace(param.ell1, _m1_, param, Ylm1);
		ToolClass::storeReducedSphericalHarmonicsInFourierSpace(param.ell2, _m2_, param, Ylm2);

		/* calc dn_tilde1 */
		DensityFieldClass<ParticleGadgetClass> dn_tilde1(param);
		double kmag1;
		double dk = kbin[1] - kbin[0];
		kmag1 = kbin[param.ith_kbin];
		dn_tilde1.calcInverseFourierTransformForBispectrum(dn_00, kmag1, dk, Ylm1);

		for(int ik = 0; ik < param.n_kbin; ik++) {
	
			double kmag2 = kbin[ik];

			/* calc dn_tilde2 */
			DensityFieldClass<ParticleGadgetClass> dn_tilde2(param);
			dn_tilde2.calcInverseFourierTransformForBispectrum(dn_00, kmag2, dk, Ylm2);
			
			/* calc bispectrum */
			std::complex<double> bk_sum = 0.0;
			double fac = param.volume / double(param.n_mesh_tot);
			for(long long coord = 0; coord < param.n_mesh_tot; coord++ ) {
				std::complex<double> f1(dn_tilde1[coord][0], dn_tilde1[coord][1]);
				std::complex<double> f2(dn_tilde2[coord][0], dn_tilde2[coord][1]);
				std::complex<double> f3(dn_LM[coord][0], dn_LM[coord][1]);
				bk_sum += (fac * f1 * f2 * f3);
			}

			bk_save[ik] += (w * bk_sum);

			double timesec = double(clock() - start);
			if(ThisTask == 0) { printf("k2 = %.3f [h/Mpc], m1 = %d, m2 = %d, M = %d| %.3f sec\n", kmag2, _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}

		}
	
	
		delete [] Ylm1; Ylm1 = NULL;
		delete [] Ylm2; Ylm2 = NULL;
		byte -= double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }


	}}

	double norm = param.volume / double(P_D.n_tot) / double(P_D.n_tot);
	norm *= (param.volume/double(P_D.n_tot));

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/bk%d%d%d_%02d.dat", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_kbin);
	fp = fopen(buf, "w");
    fprintf(fp, "## k1 [h/Mpc] \t k2 [h/Mpc] \t bk_real \t bk_imag\n");
	for(int i = 0; i < param.n_kbin; i++) {
		fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n", kbin[param.ith_kbin], kbin[i],
                norm * ( bk_save[i].real() - sn_save[i].real() ), norm * ( bk_save[i].imag() - sn_save[i].imag() ) );
	}
	fclose(fp);

	delete [] sn_save;
	delete [] bk_save;

	return 0;

}

int calcThreePointCorrelationFunctionForBOX(
        ParticleGadgetClass & P_D, ParameterClass & param, double * rbin) {

    /********************************************************************************************/
    /* This function calculates the 3pcf multipole, expanded by the TripoSH function, using FFTs.*/
    /********************************************************************************************/

	if(fabs(wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0)) < 1.0e-10) {
		if(ThisTask == 0) { 
			printf("This multipole combination is not allowed.\n");
			printf("It should be wigner_3j(ell1, ell2, ELL, 0,0,0) != 0\n");
		}
        return -1;
	}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* shotnoise */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	timesec = double(clock() - start);
	if(ThisTask == 0) { printf("Computing shotnoise terms.| %.3f sec\n", timesec / CLOCKS_PER_SEC);}

	/* define sn_save:: shotnoise term */
	std::complex<double> * sn_save = new std::complex<double>[param.n_rbin];
	for(int i = 0; i < param.n_rbin; i++) {
		sn_save[i] = 0.0;
	}

	/****************************************/
	/* calc F00 */
	DensityFieldClass<ParticleGadgetClass> N_00_sn(param);
	N_00_sn.calcNormalDensityForBispectrumShotnoiseForBOX(P_D);
	/* Fourier transform*/
	N_00_sn.calcFourierTransform();

	/**********************************************/
	/* calc the yLM-weighted density perturbation */
	DensityFieldClass<ParticleGadgetClass> dn_LM_sn(param);
	dn_LM_sn.calcNormalDensityFluctuationForBOX(P_D, param);
	/* Fourier transform*/
	dn_LM_sn.calcFourierTransform();


	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		int _M_ = 0;
	
		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}

		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.n_mesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.n_mesh_tot];
		byte += double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }
	
		ToolClass::storeReducedSphericalHarmonicsInConfigurationSpace( param.ell1, _m1_, param, Ylm1 );
		ToolClass::storeReducedSphericalHarmonicsInConfigurationSpace( param.ell2, _m2_, param, Ylm2 );
	
		
		TwoPointStatisticsClass<ParticleGadgetClass> stat(param);
		std::complex<double> shotnoise = double(P_D.n_tot);

		stat.calcCorrelationFunctionForThreePointFunction(dn_LM_sn, N_00_sn, rbin, shotnoise, param.ell1, _m1_, Ylm1, Ylm2);

		for(int i =0; i < param.n_rbin; i++) { 
			if(i == param.ith_rbin) {
				sn_save[i] += w * stat.xi[i];
			} else {
				sn_save[i] += 0.0;
			}
	    }

		/* time */
		timesec = double(clock() - start);
		if(ThisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}


		delete [] Ylm1; Ylm1 = NULL;
		delete [] Ylm2; Ylm2 = NULL;
		byte -= double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

	}}
	
	/****/
	N_00_sn.finalizeDensityField();
    dn_LM_sn.finalizeDensityField();
	/****************************************/
	
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* three-point function */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	if(ThisTask == 0) { printf("Computing the three-point function multipole.\n");}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleGadgetClass> dn_00(param);
	dn_00.calcNormalDensityFluctuationForBOX(P_D, param);
	/* Fourier transform*/
	dn_00.calcFourierTransform();

    /* calc the yLM-weighted density perturbation */
    DensityFieldClass<ParticleGadgetClass> dn_LM(param);
    dn_LM.calcNormalDensityFluctuationForBOX(P_D, param);
    /* Fourier transform*/
    dn_LM.calcFourierTransform();
    /* divided by a assignment function */
    dn_LM.calcAssignmentFunctionCorrection();
    /* Inverse Fourier transform*/
    dn_LM.calcInverseFourierTransform();

	/* store spherical bessel functions */
	SphericalBesselClass sj1(param.ell1);
	SphericalBesselClass sj2(param.ell2);

	/* define zeta_save:: bispectrum */
	std::complex<double> * zeta_save = new std::complex<double>[param.n_rbin];
	for(int i = 0; i < param.n_rbin; i++) {
		zeta_save[i] = 0.0;
	}

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		int _M_ = 0;

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}
	
		/****************************************************/
		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.n_mesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.n_mesh_tot];
		byte += double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }
	
		ToolClass::storeReducedSphericalHarmonicsInFourierSpace(param.ell1, _m1_, param, Ylm1);
		ToolClass::storeReducedSphericalHarmonicsInFourierSpace(param.ell2, _m2_, param, Ylm2);
	


		/* calc dn_tilde1 */
		DensityFieldClass<ParticleGadgetClass> dn_tilde1(param);
		double rmag1;
		rmag1 = rbin[param.ith_rbin];
		dn_tilde1.calcInverseFourierTransformForThreePointFunction(dn_00, rmag1, Ylm1, sj1);

		for(int ir = 0; ir < param.n_rbin; ir++) {
	
			double rmag2 = rbin[ir];

			/* calc dn_tilde2 */
			DensityFieldClass<ParticleGadgetClass> dn_tilde2(param);
			dn_tilde2.calcInverseFourierTransformForThreePointFunction(dn_00, rmag2, Ylm2, sj2);
			
			/* calc bispectrum */
			std::complex<double> zeta_sum = 0.0;
			double fac = param.volume / double(param.n_mesh_tot);
			std::complex<double> _I_(0.0, 1.0);
			for(long long coord = 0; coord < param.n_mesh_tot; coord++ ) {
				std::complex<double> f1(dn_tilde1[coord][0], dn_tilde1[coord][1]);
				std::complex<double> f2(dn_tilde2[coord][0], dn_tilde2[coord][1]);
				std::complex<double> f3(dn_LM[coord][0], dn_LM[coord][1]);
				zeta_sum += (pow(_I_, param.ell1+param.ell2) * fac * f1 * f2 * f3);
			}

			zeta_save[ir] += (w * zeta_sum);

			double timesec = double(clock() - start);
			if(ThisTask == 0) { printf("r2 = %.3f [Mpc/h], m1 = %d, m2 = %d, M = %d| %.3f sec\n", rmag2, _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}

		}
	
	
		delete [] Ylm1; Ylm1 = NULL;
		delete [] Ylm2; Ylm2 = NULL;
		byte -= double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }


	}}

	double norm = param.volume / double(P_D.n_tot) / double(P_D.n_tot);
	norm *= ( param.volume / double(P_D.n_tot) );

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/zeta%d%d%d_%02d.dat", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_rbin);
	fp = fopen(buf, "w");
    fprintf(fp, "## r1 [Mpc/h] \t r2 [Mpc/h] \t zeta \n");
	for(int i = 0; i < param.n_rbin; i++) {
		fprintf(fp, "%.5f \t %.5f \t %.7e\n",  rbin[param.ith_rbin], rbin[i], norm * (zeta_save[i].real() - sn_save[i].real()));
	}
    fclose(fp);

	delete [] sn_save;
	delete [] zeta_save;

	return 0;
}

int calcBiSpectrumForBOXForReconstruction(
        ParticleGadgetClass & P_D, ParticleGadgetClass & P_R, 
        ParameterClass & param, double alpha, double * kbin) {

    /************/
    /* This function calculates the bispectrum multipole, expanded by the TripoSH function, using FFTs after reconstruction. */
    /************/

	if(fabs(wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0)) < 1.0e-10) {
		if(ThisTask == 0) { 
			printf("This multipole combination is not allowed.\n");
			printf("It should be wigner_3j(ell1, ell2, ELL, 0,0,0) != 0\n");
		}
        return -1;
	}

	/*****************************************************************************/
	/* shotnoise */
	/*****************************************************************************/
	/* time */
	timesec = double(clock() - start);
	if(ThisTask == 0) { printf("Computing shotnoise terms.| %.3f sec\n", timesec / CLOCKS_PER_SEC);}

	/* define sn_save:: shotnoise term */
	std::complex<double> * sn_save = new std::complex<double>[param.n_kbin];
	for(int i = 0; i < param.n_kbin; i++) {
		sn_save[i] = 0.0;
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleGadgetClass> dn_00_sn(param);
	dn_00_sn.calcNormalDensityFluctuationForBOXForReconstruction(P_D, P_R, alpha);
	/* Fourier transform */
	dn_00_sn.calcFourierTransform();

	/**********************************************/
	/* calc the yLM-weighted density perturbation */
	DensityFieldClass<ParticleGadgetClass> N_LM_sn(param);
	N_LM_sn.calcNormalDensityForBispectrumShotnoiseForBOXForReconstruction(P_D, P_R, alpha);
	/* Fourier transform*/
	N_LM_sn.calcFourierTransform();

	/* calc shot noise terms */
	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {
		int _M_ = 0;
	
		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}
		

		TwoPointStatisticsClass<ParticleGadgetClass> stat(param);
		std::complex<double> shotnoise = stat.calcShotNoiseForBispectrumForBOXForReconstruction(P_D, P_R, alpha);

		if( (param.ell1 == 0) && (param.ell2 == 0) ) {
			for(int i =0; i < param.n_kbin; i++) { 
				sn_save[i] += w * shotnoise;
	       	}
		}

		if( (param.ell2 == 0) ) {
		    stat.calcPowerSpectrum(dn_00_sn, N_LM_sn, kbin, shotnoise, param.ell1, _m1_);
			for(int i =0; i < param.n_kbin; i++) { 
				sn_save[i] += w * stat.pk[param.ith_kbin];
			}
		}

		if( (param.ell1 == 0) ) {
			stat.calcPowerSpectrum(dn_00_sn, N_LM_sn, kbin, shotnoise, param.ell2, _m2_);
			for(int i =0; i < param.n_kbin; i++) { 
				sn_save[i] += w * stat.pk[i];
		    }
		}

		/* time */
		timesec = double(clock() - start);
		if(ThisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}

	}}
	
	/****/
	dn_00_sn.finalizeDensityField();
    N_LM_sn.finalizeDensityField();
	/****************************************/

	/****************************************/
	/* calc F00 */
	DensityFieldClass<ParticleGadgetClass> N_00_sn(param);
	N_00_sn.calcNormalDensityForBispectrumShotnoiseForBOXForReconstruction(P_D, P_R, alpha);
	/* Fourier transform*/
	N_00_sn.calcFourierTransform();

	/**********************************************/
	/* calc the yLM-weighted density perturbation */
	DensityFieldClass<ParticleGadgetClass> dn_LM_sn(param);
	dn_LM_sn.calcNormalDensityFluctuationForBOXForReconstruction(P_D, P_R, alpha);
	/* Fourier transform*/
	dn_LM_sn.calcFourierTransform();

	/* store spherical bessel functions */
	SphericalBesselClass sj1(param.ell1);
	SphericalBesselClass sj2(param.ell2);

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		int _M_ = 0;

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}

		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.n_mesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.n_mesh_tot];
		byte += double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

		ToolClass::storeReducedSphericalHarmonicsInConfigurationSpace( param.ell1, _m1_, param, Ylm1 );
		ToolClass::storeReducedSphericalHarmonicsInConfigurationSpace( param.ell2, _m2_, param, Ylm2 );


		TwoPointStatisticsClass<ParticleGadgetClass> stat(param);
		std::complex<double> shotnoise = stat.calcShotNoiseForBispectrumForBOXForReconstruction(P_D, P_R, alpha);

		fftw_complex * xi = fftw_alloc_complex(param.n_mesh_tot);
		byte += double( sizeof(fftw_complex) * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

		for(long long i = 0; i < (param.n_mesh_tot); i++) {
			xi[i][0] = 0.0;
			xi[i][1] = 0.0;
		}

		stat.calcShotNoiseForBispectrum_ijk(dn_LM_sn, N_00_sn, shotnoise, param.ELL, _M_, xi);

		for(int ik = 0; ik < param.n_kbin; ik++) {
			
			double kmag2 = kbin[ik];
			double kmag1 = 0.0;
			kmag1 = kbin[param.ith_kbin];

			/* calc shotnoise */
			std::complex<double> sn_sum = 0.0;
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
				double rmag = sqrt( rvec[0] * rvec[0] +  rvec[1] * rvec[1] +  rvec[2] * rvec[2]);
				std::complex<double> ff(xi[coord][0], xi[coord][1]);
				double j1 = sj1.getSphericalBessel(kmag1 * rmag);
				double j2 = sj2.getSphericalBessel(kmag2 * rmag);

				sn_sum += ( j1 * j2 * ff * Ylm1[coord] * Ylm2[coord] );

			}}}

			std::complex<double> _I_(0.0,1.0);
			double fac = param.volume / double(param.n_mesh_tot);
			sn_sum *= pow(_I_, param.ell1 + param.ell2) * fac;

			sn_save[ik] += ( w * sn_sum );

			/* time */
			timesec = double(clock() - start);
			if(ThisTask == 0) { printf("k2 = %.3f [h/Mpc], m1 = %d, m2 = %d, M = %d| %.3f sec\n", kmag2, _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}
		}

		fftw_free(xi); xi = NULL;
		byte -= double( sizeof(fftw_complex) * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }


		delete [] Ylm1; Ylm1 = NULL;
		delete [] Ylm2; Ylm2 = NULL;
		byte -= double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

	}}
	
	/****/
	N_00_sn.finalizeDensityField();
    dn_LM_sn.finalizeDensityField();
	/****************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* bispectrum */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	if(ThisTask == 0) { printf("computing bispectrum...\n");}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleGadgetClass> dn_00(param);
	dn_00.calcNormalDensityFluctuationForBOXForReconstruction(P_D, P_R, alpha);
	/* Fourier transform*/
	dn_00.calcFourierTransform();

    /*******************************************/
	/* calc the yLM-weighted density perturbation */
	DensityFieldClass<ParticleGadgetClass> dn_LM(param);
	dn_LM.calcNormalDensityFluctuationForBOXForReconstruction(P_D, P_R, alpha);
	/* Fourier transform*/
	dn_LM.calcFourierTransform();
	/* divided by a assignment function */
	dn_LM.calcAssignmentFunctionCorrection();
	/* Inverse Fourier transform*/
	dn_LM.calcInverseFourierTransform();

	/* define bk_save:: bispectrum */
	std::complex<double> * bk_save = new std::complex<double>[param.n_kbin];
	for(int i = 0; i < param.n_kbin; i++) {
		bk_save[i] = 0.0;
	}

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		int _M_ = 0;

		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}

		/****************************************************/
		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.n_mesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.n_mesh_tot];
		byte += double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }
	
		ToolClass::storeReducedSphericalHarmonicsInFourierSpace(param.ell1, _m1_, param, Ylm1);
		ToolClass::storeReducedSphericalHarmonicsInFourierSpace(param.ell2, _m2_, param, Ylm2);

		/* calc dn_tilde1 */
		DensityFieldClass<ParticleGadgetClass> dn_tilde1(param);
		double kmag1;
		double dk = kbin[1] - kbin[0];
		kmag1 = kbin[param.ith_kbin];
		dn_tilde1.calcInverseFourierTransformForBispectrum(dn_00, kmag1, dk, Ylm1);

		for(int ik = 0; ik < param.n_kbin; ik++) {
	
			double kmag2 = kbin[ik];

			/* calc dn_tilde2 */
			DensityFieldClass<ParticleGadgetClass> dn_tilde2(param);
			dn_tilde2.calcInverseFourierTransformForBispectrum(dn_00, kmag2, dk, Ylm2);
			
			/* calc bispectrum */
			std::complex<double> bk_sum = 0.0;
			double fac = param.volume / double(param.n_mesh_tot);
			for(long long coord = 0; coord < param.n_mesh_tot; coord++ ) {
				std::complex<double> f1(dn_tilde1[coord][0], dn_tilde1[coord][1]);
				std::complex<double> f2(dn_tilde2[coord][0], dn_tilde2[coord][1]);
				std::complex<double> f3(dn_LM[coord][0], dn_LM[coord][1]);
				bk_sum += (fac * f1 * f2 * f3);
			}

			bk_save[ik] += (w * bk_sum);

			double timesec = double(clock() - start);
			if(ThisTask == 0) { printf("k2 = %.3f [h/Mpc], m1 = %d, m2 = %d, M = %d| %.3f sec\n", kmag2, _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}

		}
	
	
		delete [] Ylm1; Ylm1 = NULL;
		delete [] Ylm2; Ylm2 = NULL;
		byte -= double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }


	}}

    /* the normalization factor */
	double norm = param.volume / double(P_D.n_tot) / double(P_D.n_tot);
	norm *= (param.volume/double(P_D.n_tot));

    /* save data */
	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/bk%d%d%d_%02d_recon.dat", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_kbin);
	fp = fopen(buf, "w");
    fprintf(fp, "## k1 [h/Mpc] \t k2 [h/Mpc] \t bk_real \t bk_imag\n");
	for(int i = 0; i < param.n_kbin; i++) {
		fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n", kbin[param.ith_kbin], kbin[i],
                norm * ( bk_save[i].real() - sn_save[i].real() ), norm * ( bk_save[i].imag() - sn_save[i].imag() ) );
	}
	fclose(fp);

	delete [] sn_save;
	delete [] bk_save;

	return 0;

}

int calcThreePointCorrelationFunctionForBOXForReconstruction(
        ParticleGadgetClass & P_D, ParticleGadgetClass & P_R, 
        ParameterClass & param, double alpha, double * rbin) {

    /********************************************************************************************/
    /* This function calculates the 3pcf multipole, expanded by the TripoSH function, using FFTs after reconstruction.*/
    /********************************************************************************************/

	if(fabs(wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0)) < 1.0e-10) {
		if(ThisTask == 0) { 
			printf("This multipole combination is not allowed.\n");
			printf("It should be wigner_3j(ell1, ell2, ELL, 0,0,0) != 0\n");
		}
        return -1;
	}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* shotnoise */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	timesec = double(clock() - start);
	if(ThisTask == 0) { printf("Computing shotnoise terms.| %.3f sec\n", timesec / CLOCKS_PER_SEC);}

	/* define sn_save:: shotnoise term */
	std::complex<double> * sn_save = new std::complex<double>[param.n_rbin];
	for(int i = 0; i < param.n_rbin; i++) {
		sn_save[i] = 0.0;
	}

	/****************************************/
	/* calc F00 */
	DensityFieldClass<ParticleGadgetClass> N_00_sn(param);
	N_00_sn.calcNormalDensityForBispectrumShotnoiseForBOXForReconstruction(P_D, P_R, alpha);
	/* Fourier transform*/
	N_00_sn.calcFourierTransform();

	/**********************************************/
	/* calc the yLM-weighted density perturbation */
	DensityFieldClass<ParticleGadgetClass> dn_LM_sn(param);
	dn_LM_sn.calcNormalDensityFluctuationForBOXForReconstruction(P_D, P_R, alpha);
	/* Fourier transform*/
	dn_LM_sn.calcFourierTransform();

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		int _M_ = 0;
	
		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}

		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.n_mesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.n_mesh_tot];
		byte += double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }
	
		ToolClass::storeReducedSphericalHarmonicsInConfigurationSpace( param.ell1, _m1_, param, Ylm1 );
		ToolClass::storeReducedSphericalHarmonicsInConfigurationSpace( param.ell2, _m2_, param, Ylm2 );
	

		TwoPointStatisticsClass<ParticleGadgetClass> stat(param);
		std::complex<double> shotnoise = stat.calcShotNoiseForBispectrumForBOXForReconstruction(P_D, P_R, alpha);

		stat.calcCorrelationFunctionForThreePointFunction(dn_LM_sn, N_00_sn, rbin, shotnoise, param.ell1, _m1_, Ylm1, Ylm2);

		for(int i =0; i < param.n_rbin; i++) { 
			if(i == param.ith_rbin) {
				sn_save[i] += w * stat.xi[i];
			} else {
				sn_save[i] += 0.0;
			}
	    }

		/* time */
		timesec = double(clock() - start);
		if(ThisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}


		delete [] Ylm1; Ylm1 = NULL;
		delete [] Ylm2; Ylm2 = NULL;
		byte -= double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

	}}
	
	/****/
	N_00_sn.finalizeDensityField();
    dn_LM_sn.finalizeDensityField();
	/****************************************/
	
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* three-point function */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	if(ThisTask == 0) { printf("Computing the three-point function multipole.\n");}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleGadgetClass> dn_00(param);
	dn_00.calcNormalDensityFluctuationForBOXForReconstruction(P_D, P_R, alpha);
	/* Fourier transform*/
	dn_00.calcFourierTransform();


	/****************************************/
	/* calc the yLM-weighted density perturbation */
	DensityFieldClass<ParticleGadgetClass> dn_LM(param);
	dn_LM.calcNormalDensityFluctuationForBOXForReconstruction(P_D, P_R, alpha);
	/* Fourier transform*/
	dn_LM.calcFourierTransform();
	/* divided by a assignment function */
	dn_LM.calcAssignmentFunctionCorrection();
	/* Inverse Fourier transform*/
	dn_LM.calcInverseFourierTransform();

	/* store spherical bessel functions */
	SphericalBesselClass sj1(param.ell1);
	SphericalBesselClass sj2(param.ell2);

	/* define zeta_save:: bispectrum */
	std::complex<double> * zeta_save = new std::complex<double>[param.n_rbin];
	for(int i = 0; i < param.n_rbin; i++) {
		zeta_save[i] = 0.0;
	}

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for(int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		int _M_ = 0;

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}
	
		/****************************************************/
		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.n_mesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.n_mesh_tot];
		byte += double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }
	
		ToolClass::storeReducedSphericalHarmonicsInFourierSpace(param.ell1, _m1_, param, Ylm1);
		ToolClass::storeReducedSphericalHarmonicsInFourierSpace(param.ell2, _m2_, param, Ylm2);
	

		/* calc dn_tilde1 */
		DensityFieldClass<ParticleGadgetClass> dn_tilde1(param);
		double rmag1;
		rmag1 = rbin[param.ith_rbin];
		dn_tilde1.calcInverseFourierTransformForThreePointFunction(dn_00, rmag1, Ylm1, sj1);

		for(int ir = 0; ir < param.n_rbin; ir++) {
	
			double rmag2 = rbin[ir];

			/* calc dn_tilde2 */
			DensityFieldClass<ParticleGadgetClass> dn_tilde2(param);
			dn_tilde2.calcInverseFourierTransformForThreePointFunction(dn_00, rmag2, Ylm2, sj2);
			
			/* calc bispectrum */
			std::complex<double> zeta_sum = 0.0;
			double fac = param.volume / double(param.n_mesh_tot);
			std::complex<double> _I_(0.0, 1.0);
			for(long long coord = 0; coord < param.n_mesh_tot; coord++ ) {
				std::complex<double> f1(dn_tilde1[coord][0], dn_tilde1[coord][1]);
				std::complex<double> f2(dn_tilde2[coord][0], dn_tilde2[coord][1]);
				std::complex<double> f3(dn_LM[coord][0], dn_LM[coord][1]);
				zeta_sum += (pow(_I_, param.ell1+param.ell2) * fac * f1 * f2 * f3);
			}

			zeta_save[ir] += (w * zeta_sum);

			double timesec = double(clock() - start);
			if(ThisTask == 0) { printf("r2 = %.3f [Mpc/h], m1 = %d, m2 = %d, M = %d| %.3f sec\n", rmag2, _m1_, _m2_, _M_, timesec / CLOCKS_PER_SEC);}

		}
	
	
		delete [] Ylm1; Ylm1 = NULL;
		delete [] Ylm2; Ylm2 = NULL;
		byte -= double( sizeof(std::complex<double>) * 2 * (param.n_mesh_tot) / 1024.0 / 1024.0 / 1024.0);
		if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }


	}}

	double norm = param.volume / double(P_D.n_tot) / double(P_D.n_tot);
	norm *= ( param.volume / double(P_D.n_tot) );

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/zeta%d%d%d_%02d_recon.dat", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_rbin);
	fp = fopen(buf, "w");
    fprintf(fp, "## r1 [Mpc/h] \t r2 [Mpc/h] \t zeta \n");
	for(int i = 0; i < param.n_rbin; i++) {
		fprintf(fp, "%.5f \t %.5f \t %.7e\n",  rbin[param.ith_rbin], rbin[i], norm * (zeta_save[i].real() - sn_save[i].real()));
	}
    fclose(fp);

	delete [] sn_save;
	delete [] zeta_save;

	return 0;
}







#endif

