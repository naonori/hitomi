
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
    if(ThisTask == 0) { printf("Start: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }

    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/

    /**************/
    /* PARAMETERS */
    /**************/
    
    timesec = double(clock() - start);
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

    /****************************/
	/* Make an output directory. */
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
    if(ThisTask == 0) { printf("Reading particle data: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }
    
    ParticleGadgetClass P_D;

    if(0) {
    } else if (param.sim_data == "Gadget") {
        /* Read the dark matter particle data from "param.data_file". */
        if ( P_D.readParticlesFromGagdetSnapshot(param.data_file) ) { 
            if(ThisTask == 0) { printf("'P_D.readParticlesFromGadgetSnapshot(param.data_file)' fails.\n");}
            EXIT();
    	}
    } else if (param.sim_data == "Rockstar") {
        /* Read the Rockstar halo data from "param.data_file". */
        if ( P_D.readParticlesFromRockstar(param.data_file, param.log10_Mmin, param.log10_Mmax) ) { 
            if(ThisTask == 0) { printf("'P_D.readParticlesFromRockstar(param.data_file, param.log10_Mmin, param.log10_Mmax)' fails.\n");}
            EXIT();
    	}
    } else {
        EXIT();
    }

    /* Make the boxsize parameter read from the input file the same as the value read from the particle data. */
    P_D.resetParameters(param);

    /* If flag_RSD=True, the RSD effect is added to the particles. */
    if(param.flag_RSD == "True") {
    	if(P_D.setRSD(param)) {
            if(ThisTask == 0) { printf("'P_D.setRSD(param)' fails.\n");}
            EXIT();
        }
    }

	/* Outputs the values of the parameters related to the simulation data. */
    ParticleGadgetClass::printParameters(param, P_D);

    timesec = double(clock() - start);
    if(ThisTask == 0) { printf("Done: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }

    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/

    /******************/
    /* RECONSTRUCTION */
    /******************/
        
    ParticleGadgetClass P_R;
    double alpha = 0.0;
	if(param.flag_recon == "True") {
    	timesec = double(clock() - start);
        if(ThisTask == 0) { printf("Calculating reconstruction: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }
        if(P_R.readRandomParticles(param.random_file)) {
            if(ThisTask == 0) { printf("'P_R.readRandomParticles(param.random_file)' fails.\n"); }
        }
		calcReconstructionParticlesForBOX(P_D, P_R, param);
        P_D.calcPeriodicBoundary(param);
        P_R.calcPeriodicBoundary(param);
		alpha =  double(P_D.n_tot) / double(P_R.n_tot); 
    
        /* Outputs the values of the parameters related to the simulation data. */
        ParticleGadgetClass::printParametersForReconstruction(param, P_D, P_R, alpha);

        timesec = double(clock() - start);
        if(ThisTask == 0) { printf("Done: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }
	}

    /************************************************************************************************************/
    /************************************************************************************************************/
    /************************************************************************************************************/

    timesec = double(clock() - start);
    if(ThisTask == 0) { printf("Calculating '%s': %.3f sec.\n", param.measure.c_str(), timesec / CLOCKS_PER_SEC); }
    /*****************************************************************************/
    /* Calculate the power spectrum multipoles expanded by the Legendre function */
    /*****************************************************************************/
    if(0) {
    } else if (param.flag_recon == "False") {

        if(0) {
        } else if (param.measure == "pk") {
            if(calcPowerSpectrumForBOX(P_D, param, kbin)) {
                if(ThisTask == 0) {printf("'calcPowerSpectrumForBOX(P_D, param, kbin)' fails.\n");}   
                EXIT();
            }
        } else if (param.measure == "2PCF") {
            if(calcTwoPointCorrelationFunctionForBOX(P_D, param, rbin)) {
                if(ThisTask == 0) {printf("'calcTwoPointCorrelationFunctionForBOX(P_D, param, rbin)' fails.\n");}
                EXIT();
            }
        } else if (param.measure == "bk") {
            if(calcBiSpectrumForBOX(P_D, param, kbin)) {
                if(ThisTask == 0) {printf("'calcBiSpectrumForBOX(P_D, param, kbin)' fails.\n");}   
                EXIT();
            }
        } else if (param.measure == "3PCF") {
            if(calcThreePointCorrelationFunctionForBOX(P_D, param, rbin)) {
                if(ThisTask == 0) {printf("'calcThreePointCorrelationFunctionForBOX(P_D, param, rbin)' fails.\n");}
                EXIT();
            }
        } else {
            EXIT();
        }

    } else if (param.flag_recon == "True") {

        if(0) {
        } else if (param.measure == "pk") {
            if(calcPowerSpectrumForBOXForReconstruction(P_D, P_R, param, alpha, kbin)) {
                if(ThisTask == 0) {printf("'calcPowerSpectrumForBOXForReconstruction(P_D, P_R, param, alpha, kbin)' fails.\n");}   
                EXIT();
            }
        } else if (param.measure == "2PCF") {
            if(calcTwoPointFunctionForBOXForReconstruction(P_D, P_R, param, alpha, rbin)) {
                if(ThisTask == 0) {printf("'calcTwoPointFunctionForBOXForReconstruction(P_D, P_R, param, alpha, rbin)' fails.\n");}   
                EXIT();
            }
        } else if (param.measure == "bk") {
            if(calcBiSpectrumForBOXForReconstruction(P_D, P_R, param, alpha, kbin)) {
                if(ThisTask == 0) {printf("'calcBiSpectrumForBOXForReconstruction(P_D, P_R, param, alpha, kbin)' fails.\n");}   
                EXIT();
            }
        } else if (param.measure == "3PCF") {
            if(calcThreePointCorrelationFunctionForBOXForReconstruction(P_D, P_R, param, alpha, rbin)) {
                if(ThisTask == 0) {printf("'calcThreePointCorrelationFunctionForBOXForReconstruction(P_D, P_R, param, alpha, rbin)' fails.\n");}   
                EXIT();
            }
        } else {
            EXIT();
        }

    } else {
        EXIT();
    } 
    
    timesec = double(clock() - start);
    if(ThisTask == 0) { printf("Done: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }

    P_D.finalizeParticle();
    P_R.finalizeParticle();
    if(ThisTask == 0) { printf("End: %.3f sec.\n", timesec / CLOCKS_PER_SEC); }

    #ifdef __MPI__
	MPI_Finalize();
    #endif

	return 0;
}

