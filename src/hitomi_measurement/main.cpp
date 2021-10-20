
#include "common.hpp"
#include "parameter.hpp"
#include "tools.hpp"
#include "bessel.hpp"
#include "particle.hpp"
#include "particleGadget.hpp"
#include "particle_reconstruction.hpp"
#include "field.hpp"
#include "powerspec.hpp"
#include "bispec.hpp"

int main(int argc, char *argv[]) {

    /************************************************************/
    /* The default for this code is the serial code */
    /* It is assumed that MPI will be implemented on your own. */
    /************************************************************/

    #ifdef __MPI__
	/* MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
	MPI_Comm_size(MPI_COMM_WORLD, &NTask);
    #endif

    /* time */
	start = clock();
    timesec = double(clock() - start);
    if(ThisTask == 0) { printf("\n"); }
    if(ThisTask == 0) { printf("Start: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }

    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/

    /**************/
    /* PARAMETERS */
    /**************/
    
    timesec = double(clock() - start);
    if(ThisTask == 0) { printf("\n"); }
    if(ThisTask == 0) { printf("Reading input parameters: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }
	/**********************************************************/
	/* If there is no input parameter file, exit the program. */
	/**********************************************************/
	if(argc != 2){
        if(ThisTask == 0) {
	    	printf("The input parameter file is not found.\n");
	    	printf("For example, run './a.out param.ini'\n");
        }
        EXIT();
	}

	/****************************/
	/* Read the parameter file. */
	/****************************/
	ParameterClass param;
	if( param.readParameters(argv) ) {
        if(ThisTask == 0) { printf("'param.readParameters' fails.\n"); }
        EXIT();
	}

    /*************************************************************/
    /* We recommend using MPI parallelism for the input parameter values, because then there is no data communication between MPIs and the program is simpler.
    *  As an example, we consider reading a number of different realizations of mock catalogs using MPI: 
    *  for Ntask=20, NR=0, this code reads the 1-20th mock samples; for NTask=20, NR=1, the code reads the 21-40th mock samples.
    *  In this way, the 'NR' parameter controls the number of mock samples to load.*/
    /*************************************************************/

//    #ifdef __MPI__
//    param.realization = NTask * param.NR + ThisTask + 1;
//    char buf[1024];
//    sprintf(buf, "%s/%s_%04d.dat", param.data_dir.c_str(), param.data_mpi_file.c_str(), param.realization);
//    param.data_file = buf;
//    #endif

    /*************************************************************/
    /* As another example, you can use MPI for the parameter 'ith_kbin' that specifies k1 out of k1 and k2 on which the bispectrum depends. 
     * In this case, when NTask = 20 and n_kbin=20, it is possible to calculate all the scale dependencies of the bispectrum in a single calculation. */
    /*************************************************************/
//    #ifdef __MPI__
//    param.ith_rbin = ThisTask;
//    #endif

    /*************************************************************/
    /* You can also implement MPI parallelism for any other input parameters that are convenient for you.*/
    /*************************************************************/

    /****************************/
	/* Make output directories. */
	/****************************/
    if( param.makeOutputDirectory() ) {
        if(ThisTask == 0) { printf("'param.makeOutputDirectory' fails.\n"); }
        EXIT();
    }

   	/**********************************************************************************************/
	/* Outputs the input parameters read from "param.ini" to "param.output_dir/input_parameters". */
   	/**********************************************************************************************/
    if ( param.printParameters() ) {
        if(ThisTask == 0) { printf("'param.printParameters' fails.\n"); }
        EXIT();
    }

    /*************************/
    /* Set k-bins and r-bins */
    /*************************/
	double kbin[param.n_kbin];
	double rbin[param.n_rbin];
	ToolClass::setKbin(param, kbin);
	ToolClass::setRbin(param, rbin);

    timesec = double(clock() - start);
    if(ThisTask == 0) { printf("Done: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }

    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/

    /*************/
    /* PARTICLES */
    /*************/

	timesec = double(clock() - start);
    if(ThisTask == 0) { printf("\n"); }
    if(ThisTask == 0) { printf("Reading particle data: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }
    /**********************************************************************************************************/
    /* Read the galaxy and random particle data from files "param.data_file" and "random_file", respectively. */
    /**********************************************************************************************************/
	ParticleBOSSClass P_D, P_R;
	if ( P_D.readParticleBOSS(param.data_file) ) { 
        if(ThisTask == 0) { printf("'P_D.readParticleBOSS' fails.\n"); }
        EXIT();
	}
	if ( P_R.readParticleBOSS(param.random_file) ) { 
        if(ThisTask == 0) { printf("'P_R.readParticleBOSS' fails.\n"); }
        EXIT();
	}

    /********************************************************/
    /* Calculate the total number of particles with weights.*/
    /********************************************************/
	if ( P_D.calcTotalNumberOfParticlesWithWeights() ) { 
        if(ThisTask == 0) { printf("'P_D.calcTotalNumberOfParticlesWithWeight()' fails.\n"); }
        EXIT();
    }
	if ( P_R.calcTotalNumberOfParticlesWithWeights() ) { 
        if(ThisTask == 0) { printf("'P_R.calcTotalNumberOfParticlesWithWeight()' fails.\n"); }
        EXIT();
    }

    /***************************************************************************/
    /* Calculate the unit vector of the line-of-sight direction to the galaxy.*/
    /***************************************************************************/
	if ( P_D.calcLOS() ) { 
        if(ThisTask == 0) { printf("'P_D.calcLOS()' fails.\n"); }
        EXIT();
    }
	if ( P_R.calcLOS() ) { 
        if(ThisTask == 0) { printf("'P_R.calcLOS()' fails.\n"); }
        EXIT();
    }

    /*********************************************************************************************************/
    /* Calculate the ratio of the number of weighted galaxies to the number of weighted random particles, i.e., "alpha" */
    /*********************************************************************************************************/
    double alpha = ParticleBOSSClass::calcAlpha(P_D, P_R);
    if( alpha < 0.0 ) {
        if(ThisTask == 0) { printf("'ParticleBOSSClass::calcAlpha(P_D, P_R)' fails.\n");}
        EXIT();
    }
	
    /*****************************************/
    /* Replace the particles in the FFT box. */
    /*****************************************/
    if(ParticleBOSSClass::resetParticleForFFT(P_D, P_R, param)) {
        if(ThisTask == 0) { printf("'ParticleBOSSClass::resetParticleForFFT(P_D, P_R, param)' fails.\n");}
        EXIT();
    }

    /***************************************************************************/
    /* Calculate the survey volume from the number density of random particles */
    /***************************************************************************/
	DensityFieldClass<ParticleBOSSClass> Vol(param);
	double Vsurvey = Vol.calcSurveyVolume(P_R);
	Vol.finalizeDensityField();

    /**************************************************/
    /* Calculate the mean number density of particles */
    /**************************************************/
    double nbar = ParticleBOSSClass::calcMeanNumberDensity(P_D, Vsurvey);
    if(nbar < 0.0) {
        if(ThisTask == 0) { printf("'ParticleBOSSClass::calcMeanNumberDensity(P_D, Vsurvey)' fails.\n");}
        EXIT();
    }

    /************************************************************************************************/
	/* Outputs the values of the survey-related parameters to "param.output_dir/survey_parameters". */
    /************************************************************************************************/
    if(ParticleBOSSClass::printParameters(param, P_D, P_R, alpha, Vsurvey, nbar)) {
        if(ThisTask == 0) {printf("'ParticleBOSSClass::printParameters(param, P_D, P_R, alpha, Vsurvey, nbar))' fails.\n");} 
        EXIT();
    }
    
    timesec = double(clock() - start);
    if(ThisTask == 0) { printf("Done: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }

    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/

    /******************/
    /* RECONSTRUCTION */
    /******************/
	if(param.flag_recon == "True") {
    	timesec = double(clock() - start);
        if(ThisTask == 0) { printf("\n"); }
        if(ThisTask == 0) { printf("Calculating reconstruction: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }
		calcReconstructionParticles(P_D, P_R, param, alpha, Vsurvey);
		ParticleBOSSClass::resetParticleForFFT(P_D, P_R, param);

        timesec = double(clock() - start);
        if(ThisTask == 0) { printf("Done: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }
	}

    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/

    /**********************/
    /* N-POINT STATISTICS */
    /**********************/

    timesec = double(clock() - start);
    if(ThisTask == 0) { printf("\n"); }
    if(ThisTask == 0) { printf("Calculating '%s': %.3f sec.\n", param.measure.c_str(), timesec / CLOCKS_PER_SEC); }
    /************************************************************/
    /* Calculate the power spectrum, 2pcf, bispectrum, or 3pcf. */
    /************************************************************/
    if(0) {
    } else if(param.measure == "pk") {
        if(calcPowerSpectrum(P_D, P_R, param, alpha, kbin, Vsurvey)) {
            if(ThisTask == 0) {printf("'calcPowerSpectrum(P_D, P_R, param, alpha, kbin, Vsurvey)' fails.\n");}   
            EXIT();
        }
    } else if(param.measure == "2PCF") {
        if(calcTwoPointCorrelationFunction(P_D, P_R, param, alpha, rbin, Vsurvey)) {
            if(ThisTask == 0) {printf("'calcTwoPointFunction(P_D, P_R, param, alpha, rbin, Vsurvey)' fails.\n");}   
            EXIT();
        }
    } else if(param.measure == "window2PCF") {
        if(calcTwoPointWindowCorrelationFunction(P_R, param, rbin, Vsurvey)) {
            if(ThisTask == 0) {printf("'calcTwoPointWindowCorrelationFunction(P_R, param, rbin, Vsurvey)' fails.\n");}   
            EXIT();
        }
    } else if(param.measure == "bk") {
        if(calcBiSpectrum(P_D, P_R, param, alpha, kbin, Vsurvey)) {
            if(ThisTask == 0) {printf("'calcBiSpectrum(P_D, P_R, param, alpha, kbin, Vsurvey)' fails.\n");}   
            EXIT();
        }
    } else if(param.measure == "3PCF") {
        if(calcThreePointCorrelationFunction(P_D, P_R, param, alpha, rbin, Vsurvey)) {
            if(ThisTask == 0) {printf("'calcThreePointCorrelationFunction(P_D, P_R, param, alpha, kbin, Vsurvey)' fails.\n");}   
            EXIT();
        }
    } else if(param.measure == "window3PCF") {
        if(calcThreePointWindowCorrelationFunction(P_R, param, rbin, Vsurvey)) {
            if(ThisTask == 0) {printf("'calcThreePointWindowCorrelationFunction(P_R, param, rbin, Vsurvey)' fails.\n");}   
            EXIT();
        }
    } else {
        EXIT();
    }
    timesec = double(clock() - start);
    if(ThisTask == 0) { printf("Done: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }

/******************************************/

    P_D.finalizeParticle();
    P_R.finalizeParticle();
    if(ThisTask == 0) { printf("\n"); }
    if(ThisTask == 0) { printf("End: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }

/******************************************/

    #ifdef __MPI__
	MPI_Finalize();
    #endif

	return 0;
}

