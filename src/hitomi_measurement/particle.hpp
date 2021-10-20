#ifndef __particleBOSS__
#define __particleBOSS__

#ifndef __parameter__
#include "parameter.hpp"
#endif

class ParticleBOSSClass {
private:
public:
    /* Particle information */
	struct ParticleInfo {
        double pos[3];
        double los[3];
        double w;
	} * P;

    /* The total number of particles */
	long long n_tot;
	double n_tot_weight;
    /* The maximum value of particle positions */
	double pos_max[3];
    /* The minimum value of particle positions */
	double pos_min[3];

	ParticleInfo & operator [] (long long id) { return this->P[id]; }

	ParticleBOSSClass() {
		/* Initialize */
		this->P = NULL;
		this->n_tot = 0;
		this->n_tot_weight = 0;
		this->pos_max[0] = 0.0; this->pos_min[0] = 0.0;
		this->pos_max[1] = 0.0; this->pos_min[1] = 0.0;
		this->pos_max[2] = 0.0; this->pos_min[2] = 0.0;
	}
	~ParticleBOSSClass() {
		this->finalizeParticle();
	}
	
	int initializeParticle(const long long num) {

        /*************************************************************/
        /* This function allocates memory for the particle information 
         * and then initializes the particle information.*/
        /*************************************************************/

        /* If the total number of particles is less than or equal to zero, 
         * the function returns an error.*/
		if(num <= 0) { 
            printf("The number of particles should be > 0\n"); 
            return -1; 
        }

		/* Insert the total number of particles in this->n_tot */
		this->n_tot = num;

		/* Allocate memory for particles */
		delete [] this->P; this->P = NULL;
		this->P = new ParticleInfo[this->n_tot];
        byte += double( sizeof(ParticleInfo) * (this->n_tot) / 1024.0 / 1024.0 / 1024.0);
        if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }

		/* Initialize the particle information */
		for(long long i = 0; i < this->n_tot; i++) {
			this->P[i].pos[0] = 0.0;
			this->P[i].pos[1] = 0.0;
			this->P[i].pos[2] = 0.0;
			this->P[i].los[0] = 0.0;
			this->P[i].los[1] = 0.0;
			this->P[i].los[2] = 0.0;		
            this->P[i].w = 0.0;
		}

        return 0;

	}

	void finalizeParticle() {

        /***************************************************/
        /* This function finalizes the particle information.*/
        /***************************************************/

		if(this->P != NULL) {
			delete [] this->P; this->P = NULL; 
            byte -= double( sizeof(ParticleInfo) * (this->n_tot) / 1024.0 / 1024.0 / 1024.0);
            if(ThisTask == 0) { printf("Memory = %f Gb\n", byte); }
		}

	}

	int readParticleBOSS(std::string & fname_in) {

        /*************************************************************************************************/
        /* This function reads a particle data. 
         * It requires the particle data to be written in Cartan coordinates, 
         * and the first, second, and third columns of the file represent the x, y, and z coordinates. 
         * The fourth column is the weight function (w) to be entered when calculating the number density.
         * Therefore, the particle data file must be written in four columns: x,y,z,w. */
        /*************************************************************************************************/

		std::ifstream fin;

		/*****************************************/
		/* Count the number of lines (particles) */
		int num_lines = 0; 
		fin.open(fname_in.c_str(), std::ios::in);
		if( fin.fail() ) {
			if(ThisTask == 0) { printf("Cannot open the file '%s'...\n", fname_in.c_str()); }
			fin.close();
			return -1;
		}
		std::string str;
		double x, y, z, w;
		while(getline(fin, str)) {
			if( sscanf(str.c_str(), "%lf %lf %lf %lf", &x, &y, &z, &w) != 4 ) {
				continue;
			}
			num_lines++;
		}
		fin.close();
		/*****************************************/

		/*****************************************/
		/* Initialize the particle information */
        if(this->initializeParticle(num_lines)) {
            if(ThisTask == 0) { printf("'initializeParticle' fails.\n"); }
            return -1;
        }
		/*****************************************/

		/*****************************************/
		/* Read the input particle data */
		num_lines = 0;
		fin.open(fname_in.c_str(), std::ios::in);
		while(getline(fin, str)) {
			if( sscanf(str.c_str(), "%lf %lf %lf %lf", &x, &y, &z, &w) != 4 ) {
				continue;
			}
			this->P[num_lines].pos[0] = x;
			this->P[num_lines].pos[1] = y;
			this->P[num_lines].pos[2] = z;
			this->P[num_lines].w = w;
			num_lines++;
		}
		fin.close();
		/*****************************************/

        this->calcMinAndMax();

		return 0;

	}

	int calcMinAndMax() {
        
        /*************************************************************************************/
        /* This function calculates the maximum and minimum values of the particle positions. */
        /*************************************************************************************/
		
		if( this->P == NULL ) { 
            if(ThisTask == 0) {printf("'calcMinAndMax()' fails.\n"); }     
            return -1;
        }

		double min[3], max[3];

		min[0] = this->P[0].pos[0]; max[0] = this->P[0].pos[0];
		min[1] = this->P[0].pos[1]; max[1] = this->P[0].pos[1];
		min[2] = this->P[0].pos[2]; max[2] = this->P[0].pos[2];

		for (long long i = 0; i < this->n_tot; i++) {
			if(min[0] > this->P[i].pos[0]) {
				min[0] = this->P[i].pos[0];
			}
			if(min[1] > this->P[i].pos[1]) {
				min[1] = this->P[i].pos[1];
			}
			if(min[2] > this->P[i].pos[2]) {
				min[2] = this->P[i].pos[2];
			}

			if(max[0] < this->P[i].pos[0]) {
				max[0] = this->P[i].pos[0];
			}
			if(max[1] < this->P[i].pos[1]) {
				max[1] = this->P[i].pos[1];
			}
			if(max[2] < this->P[i].pos[2]) {
				max[2] = this->P[i].pos[2];
			}
		}
		this->pos_min[0] = min[0]; this->pos_max[0] = max[0];
		this->pos_min[1] = min[1]; this->pos_max[1] = max[1];
		this->pos_min[2] = min[2]; this->pos_max[2] = max[2];
	
		return 0;

	}

	int resetParticle(const double * dP) {

        /*******************************************************************/
        /* This function replaces the particle data in the box for the FFT. */
        /*******************************************************************/

		if( this->P == NULL ) { 
            if(ThisTask == 0) {printf("'resetParticle' fails.\n"); }     
            return -1;
        }

		for (long long p = 0; p < this->n_tot; p++) {
			this->P[p].pos[0] -= dP[0];	
			this->P[p].pos[1] -= dP[1];	
			this->P[p].pos[2] -= dP[2];	
		}

		return 0;
	}

    static int resetParticleForFFT(ParticleBOSSClass & P_D, ParticleBOSSClass & P_R, ParameterClass & param, double factor = 1.0) {

        /***************************************************************************************/
        /* This function replaces the particle data in the box for the FFT.
         * The factor=3 means that the particles will be placed inside the FFT box for 3 meshes. 
         * This is to guarantee that the particles will always be inside the FFT box, 
         * but in principle factor=0 is also acceptable. 
         * factor=3 is just for peace of mind. */
        /***************************************************************************************/

        if(P_D.calcMinAndMax()) {
            if(ThisTask == 0) {printf("'P_D.calcMinAndMax()' fails.\n"); }     
            return -1;
        }
        if(P_R.calcMinAndMax()) {
            if(ThisTask == 0) {printf("'P_R.calcMinAndMax()' fails.\n"); }     
            return -1;
        }

		double dP[3] = {P_R.pos_min[0], P_R.pos_min[1], P_R.pos_min[2]};
		dP[0] -= factor * param.boxsize[0] / double(param.n_mesh[0]);
		dP[1] -= factor * param.boxsize[1] / double(param.n_mesh[1]);
		dP[2] -= factor * param.boxsize[2] / double(param.n_mesh[2]);
	
        if(P_D.resetParticle(dP)) {
            if(ThisTask == 0) {printf("'P_D.resetParticle(dP)' fails.\n"); }     
            return -1;
        }
        if(P_R.resetParticle(dP)) {
            if(ThisTask == 0) {printf("'P_R.resetParticle(dP)' fails.\n"); }     
            return -1;
        }

        if(P_D.calcMinAndMax()) {
            if(ThisTask == 0) {printf("'P_D.calcMinAndMax()' fails.\n"); }     
            return -1;
        }
        if(P_R.calcMinAndMax()) {
            if(ThisTask == 0) {printf("'P_R.calcMinAndMax()' fails.\n"); }     
            return -1;
        }

		return 0;
	}


	static double calcAlpha(ParticleBOSSClass & P_D, ParticleBOSSClass & P_R) {

        /*****************************************************************************************************/
        /* This function calculates the ratio, alpha, of the number of weighted galaxies to the number of weighted random particles. */
        /*****************************************************************************************************/
        
        if( (P_R.n_tot_weight == 0) || (P_D.n_tot_weight == 0)) {
            return -1;
        }

		double alpha =  P_D.n_tot_weight / P_R.n_tot_weight; 
	
		return alpha;

	}

    static double calcMeanNumberDensity(ParticleBOSSClass & P_D, double Vsurvey) {

        if( Vsurvey == 0.0 ) {
            return -1;
        }

        double num_D_weight = 0.0;
        for(long long p = 0; p < P_D.n_tot; p++) {
            num_D_weight += pow(P_D[p].w, 2);
        }

		double nbar =  P_D.n_tot_weight * P_D.n_tot_weight / Vsurvey / num_D_weight; 
	
		return nbar;

	}


	static double calcNormalizationForPowerSpectrum(ParticleBOSSClass & P_D, double Vsurvey) {
	
        /************************************************************************************/
        /* This function calculates the normalization factor for the power spectrum (2PCF). */
        /************************************************************************************/

        if(P_D.n_tot_weight == 0) {
            return -1;
        }

		double norm = Vsurvey / P_D.n_tot_weight / P_D.n_tot_weight; 
		return norm;

	}



	static double calcNormalizationForBiSpectrum(ParticleBOSSClass & P_D, double Vsurvey) {
	
        /********************************************************************************/
        /* This function calculates the normalization factor for the bispectrum (3PCF). */
        /********************************************************************************/
        if(P_D.n_tot_weight == 0) {
            return -1;
        }

        double norm = Vsurvey / P_D.n_tot_weight / P_D.n_tot_weight; 
		norm *= (Vsurvey / P_D.n_tot_weight);
		return norm;

	}

	int calcTotalNumberOfParticlesWithWeights() {

        /************************************************************************/
        /* This function calculates the total number of particles with weights. */
        /************************************************************************/

		if( this->P == NULL ) { 
            if(ThisTask == 0) {printf("'calcTotalNumberOfParticlesWithWeight' fails.\n"); }     
            return -1;
        }
	
		double num_weight = 0.0;
		for(long long p = 0; p < this->n_tot; p++) {
			num_weight += this->P[p].w;
		}

		this->n_tot_weight = num_weight;
	
        return 0;
	}

    int calcLOS() {

        /****************************************************************************************/
        /* This function calculates the unit vector of the line-of-sight direction to the galaxy.*/
        /****************************************************************************************/
 
        if( this->P == NULL ) { 
            if(ThisTask == 0) {printf("'calcLOS' fails.\n"); }     
            return -1;
        }
   
        for (long long p = 0; p < this->n_tot; p++) { 
            double pos_mag = sqrt( this->P[p].pos[0] * this->P[p].pos[0] + this->P[p].pos[1] * this->P[p].pos[1] + this->P[p].pos[2] * this->P[p].pos[2] );
            this->P[p].los[0] = this->P[p].pos[0] / pos_mag;
            this->P[p].los[1] = this->P[p].pos[1] / pos_mag;
            this->P[p].los[2] = this->P[p].pos[2] / pos_mag;
        }

        return 0;
    
    }

	static int printParameters(ParameterClass & param, ParticleBOSSClass & P_D, ParticleBOSSClass & P_R, double alpha, double Vsurvey, double nbar) {

        /*******************************************************************************/
        /* This function outputs the survey-related parameters to the output directory */
        /*******************************************************************************/

        FILE *fp;
        char buf[1024];
        sprintf(buf, "%s/survey_parameters", param.output_dir.c_str());
        if( !(fp = fopen(buf,"w")) ) {
        	if(ThisTask == 0) { printf("output directry '%s' does not exist\n", param.output_dir.c_str()); }
        	return -1;
        }
        
        fprintf(fp, "\n");
        fprintf(fp, "P_D.n_tot = %lld\n", P_D.n_tot);
        fprintf(fp, "P_D.n_tot_weight = %f\n", P_D.n_tot_weight);
        fprintf(fp, "P_R.n_tot = %lld\n", P_R.n_tot);
        fprintf(fp, "P_R.n_tot_weight = %f\n", P_R.n_tot_weight);
        fprintf(fp, "\n");
        fprintf(fp, "alpha = %.5f\n", alpha);
        fprintf(fp, "Vsurvey [(Gpc/h)^3] = %.5f\n", Vsurvey/1.0e9);
        fprintf(fp, "Vsurvey^(1/3) [Mpc/h] = %.5f\n", pow(Vsurvey,1.0/3.0));
        fprintf(fp, "nbar / 1.0e-4 [(Mpc/h)^{-3}] = %.5f\n", nbar/1.0e-4);
        fprintf(fp, "\n");
        fprintf(fp, "P_D.pos_min[0] [Mpc/h] = %.5f\n", P_D.pos_min[0]);
        fprintf(fp, "P_D.pos_min[1] [Mpc/h] = %.5f\n", P_D.pos_min[1]);
        fprintf(fp, "P_D.pos_min[2] [Mpc/h] = %.5f\n", P_D.pos_min[2]);
        fprintf(fp, "P_R.pos_min[0] [Mpc/h] = %.5f\n", P_R.pos_min[0]);
        fprintf(fp, "P_R.pos_min[1] [Mpc/h] = %.5f\n", P_R.pos_min[1]);
        fprintf(fp, "P_R.pos_min[2] [Mpc/h] = %.5f\n", P_R.pos_min[2]);
        fprintf(fp, "\n");
        fprintf(fp, "P_D.pos_max[0] [Mpc/h] = %.5f\n", P_D.pos_max[0]);
        fprintf(fp, "P_D.pos_max[1] [Mpc/h] = %.5f\n", P_D.pos_max[1]);
        fprintf(fp, "P_D.pos_max[2] [Mpc/h] = %.5f\n", P_D.pos_max[2]);
        fprintf(fp, "P_R.pos_max[0] [Mpc/h] = %.5f\n", P_R.pos_max[0]);
        fprintf(fp, "P_R.pos_max[1] [Mpc/h] = %.5f\n", P_R.pos_max[1]);
        fprintf(fp, "P_R.pos_max[2] [Mpc/h] = %.5f\n", P_R.pos_max[2]);
        fprintf(fp, "\n");
        
        fclose(fp);

		return 0;

	}

};


#endif
