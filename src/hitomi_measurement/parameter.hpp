#ifndef __parameter__
#define __parameter__

#ifndef __common__
#include "common.hpp"
#endif

class ParameterClass {

    private:
    public:

        double boxsize[3];
        double volume;
        int n_mesh[3];
        long long n_mesh_tot;

        int ell1;
        int ell2;
        int ELL;

        std::string data_file;
        std::string random_file;
        std::string output_dir;
        std::string data_dir;
        std::string data_mpi_file;

        std::string measure;

        double kmin;
        double kmax;
        int n_kbin;

        double rmin;
        double rmax;
        int n_rbin;

        std::string assign;

        int ith_kbin;
        int ith_rbin;

        int NR;

        double b1_fid;
        double RG;
        std::string flag_recon; 
        int n_mesh_recon[3];

        std::string flag_RSD; 
        std::string sim_data; 
        double log10_Mmin; 
        double log10_Mmax; 

        int realization;

        ParameterClass() {
        
            this->boxsize[0] = 0.0;
            this->boxsize[1] = 0.0;
            this->boxsize[2] = 0.0;
            this->volume = 0.0;
            this->n_mesh[0] = 0;
            this->n_mesh[1] = 0;
            this->n_mesh[2] = 0;
            this->n_mesh_tot = 0;
    
            this->ell1 = -1;
            this->ell2 = -1;
            this->ELL = -1;
    
            this->data_file = "None";
            this->random_file = "None";
            this->output_dir = "None";
            this->data_dir = "None";
            this->data_mpi_file = "None";

            this->measure = "None";
    
            this->kmin = 0.0;
            this->kmax = 0.0;
            this->n_kbin = 1;
    
            this->rmin = 0.0;
            this->rmax = 0.0;
            this->n_rbin = 1;
    
            this->assign = "None";
    
            this->ith_kbin = -1;
            this->ith_rbin = -1;
    
            this->NR = -1;
    
            this->b1_fid = 0.0;
            this->RG = -1.0;
            this->flag_recon = "None"; 
            this->n_mesh_recon[0] = 0;
            this->n_mesh_recon[1] = 0;
            this->n_mesh_recon[2] = 0;
    
            this->flag_RSD = "None"; 
            this->sim_data = "None"; 
            this->log10_Mmin = 0.0; 
            this->log10_Mmax = 0.0; 
    
            this->realization = -1;

        }

        int readParameters(char * argv[]) {

            /*************************************************/
            /* This function reads the input parameter file  */
            /*************************************************/

            std::string param_file = argv[1];
            std::ifstream fin(param_file.c_str());
            
            if(!fin) {
                if(ThisTask == 0) {printf("cannot open '%s'\n", param_file.c_str()); }
                return -1;
            }

            std::string str = "None";
            double boxsize_x = 0.0;
            double boxsize_y = 0.0;
            double boxsize_z = 0.0;
            int n_mesh_x = 0;
            int n_mesh_y = 0;
            int n_mesh_z = 0;
            int n_mesh_recon_x = 0;
            int n_mesh_recon_y = 0;
            int n_mesh_recon_z = 0;
            char dummy[256] = "None";
            char dummy_check_temp[256] = "None";
            char data_file_temp[1024] = "None";
            char data_mpi_file_temp[1024] = "None";
            char random_file_temp[1024] = "None";
            char data_dir_temp[1024] = "None";
            char output_dir_temp[1024] = "None";
            char measure_temp[1024] = "None";
            char assign_temp[1024] = "None";
            char flag_recon_temp[1024] = "None";
            char flag_RSD_temp[1024] = "None";
            char sim_data_temp[1024] = "None";
            unsigned long c = std::string::npos;

            while(getline(fin, str)){

                if(sscanf(str.data(), "%s %s %s", dummy, dummy, dummy) != 3) {
                    continue;
                }

                sscanf(str.data(), "%s %s %s", dummy, dummy_check_temp, dummy);
                std::string dummy_check = dummy_check_temp;
                if(dummy_check != "=") {
                    continue;
                }

                if(str.find("realization") != c)  { sscanf(str.data(), "%s %s %d", dummy, dummy, &this->realization);}

                if(str.find("boxsize_x") != c)    { sscanf(str.data(), "%s %s %lg", dummy, dummy, &boxsize_x);}
                if(str.find("boxsize_y") != c)    { sscanf(str.data(), "%s %s %lg", dummy, dummy, &boxsize_y);}
                if(str.find("boxsize_z") != c)    { sscanf(str.data(), "%s %s %lg", dummy, dummy, &boxsize_z);}

                if(str.find("n_mesh_x") != c)     { sscanf(str.data(), "%s %s %d",  dummy, dummy, &n_mesh_x);}
                if(str.find("n_mesh_y") != c)     { sscanf(str.data(), "%s %s %d",  dummy, dummy, &n_mesh_y);}
                if(str.find("n_mesh_z") != c)     { sscanf(str.data(), "%s %s %d",  dummy, dummy, &n_mesh_z);}

                if(str.find("ell1") != c)         { sscanf(str.data(), "%s %s %d",  dummy, dummy, &this->ell1);}
                if(str.find("ell2") != c)         { sscanf(str.data(), "%s %s %d",  dummy, dummy, &this->ell2);}
                if(str.find("ELL") != c)          { sscanf(str.data(), "%s %s %d",  dummy, dummy, &this->ELL);}

                if(str.find("kmin") != c)         { sscanf(str.data(), "%s %s %lg", dummy, dummy, &this->kmin);}
                if(str.find("kmax") != c)         { sscanf(str.data(), "%s %s %lg", dummy, dummy, &this->kmax);}
                if(str.find("n_kbin") != c)       { sscanf(str.data(), "%s %s %d",  dummy, dummy, &this->n_kbin);}

                if(str.find("rmin") != c)         { sscanf(str.data(), "%s %s %lg", dummy, dummy, &this->rmin);}
                if(str.find("rmax") != c)         { sscanf(str.data(), "%s %s %lg", dummy, dummy, &this->rmax);}
                if(str.find("n_rbin") != c)       { sscanf(str.data(), "%s %s %d",  dummy, dummy, &this->n_rbin);}

                if(str.find("data_file") != c)    { sscanf(str.data(), "%s %s %s",  dummy, dummy, data_file_temp);}
                if(str.find("random_file") != c)  { sscanf(str.data(), "%s %s %s",  dummy, dummy, random_file_temp);}
                if(str.find("data_dir") != c)     { sscanf(str.data(), "%s %s %s",  dummy, dummy, data_dir_temp);}
                if(str.find("output_dir") != c)   { sscanf(str.data(), "%s %s %s",  dummy, dummy, output_dir_temp);}
                if(str.find("data_mpi_file") != c)    { sscanf(str.data(), "%s %s %s",  dummy, dummy, data_mpi_file_temp);}

                if(str.find("measure") != c)    { sscanf(str.data(), "%s %s %s",  dummy, dummy, measure_temp);}

                if(str.find("assignment") != c)   { sscanf(str.data(), "%s %s %s",  dummy, dummy, assign_temp);}
                if(str.find("ith_kbin") != c)     { sscanf(str.data(), "%s %s %d",  dummy, dummy, &this->ith_kbin);}
                if(str.find("ith_rbin") != c)     { sscanf(str.data(), "%s %s %d",  dummy, dummy, &this->ith_rbin);}

                if(str.find("NR") != c)           { sscanf(str.data(), "%s %s %d",  dummy, dummy, &this->NR);}

                if(str.find("b1_fid") != c)       { sscanf(str.data(), "%s %s %lg", dummy, dummy, &this->b1_fid);}
                if(str.find("RG") != c)           { sscanf(str.data(), "%s %s %lg", dummy, dummy, &this->RG);}
                if(str.find("flag_recon") != c)        { sscanf(str.data(), "%s %s %s",  dummy, dummy, flag_recon_temp);}
                if(str.find("n_mesh_recon_x") != c)     { sscanf(str.data(), "%s %s %d",  dummy, dummy, &n_mesh_recon_x);}
                if(str.find("n_mesh_recon_y") != c)     { sscanf(str.data(), "%s %s %d",  dummy, dummy, &n_mesh_recon_y);}
                if(str.find("n_mesh_recon_z") != c)     { sscanf(str.data(), "%s %s %d",  dummy, dummy, &n_mesh_recon_z);}

                if(str.find("flag_RSD") != c)        { sscanf(str.data(), "%s %s %s",  dummy, dummy, flag_RSD_temp);}
                if(str.find("sim_data") != c)        { sscanf(str.data(), "%s %s %s",  dummy, dummy, sim_data_temp);}
            if(str.find("log10_Mmin") != c)      { sscanf(str.data(), "%s %s %lg", dummy, dummy, &this->log10_Mmin);}
            if(str.find("log10_Mmax") != c)      { sscanf(str.data(), "%s %s %lg", dummy, dummy, &this->log10_Mmax);}

        }
		
        /**************/
		this->boxsize[0] = boxsize_x;
		this->boxsize[1] = boxsize_y;
		this->boxsize[2] = boxsize_z;

		this->n_mesh[0] = n_mesh_x;
		this->n_mesh[1] = n_mesh_y;
		this->n_mesh[2] = n_mesh_z;

		this->n_mesh_recon[0] = n_mesh_recon_x;
		this->n_mesh_recon[1] = n_mesh_recon_y;
		this->n_mesh_recon[2] = n_mesh_recon_z;

		this->data_file = data_file_temp;
		this->data_mpi_file = data_mpi_file_temp;
		this->random_file = random_file_temp;
		this->data_dir = data_dir_temp;
		this->output_dir = output_dir_temp;
		this->measure = measure_temp;
		this->assign = assign_temp;

		this->volume = boxsize_x * boxsize_y * boxsize_z;
		this->n_mesh_tot = n_mesh_x * n_mesh_y * n_mesh_z;

		this->flag_recon = flag_recon_temp;
		this->flag_RSD = flag_RSD_temp;
		this->sim_data = sim_data_temp;
        /**************/

        /***************************************************************/
        /* Checks whether the read parameters have reasonable values. */
        /***************************************************************/
        if( this->checkParameters() ) {
            return -1;
        }

        /*************/
        this->data_file = this->data_dir + "/" + this->data_file;
        this->random_file = this->data_dir + "/" + this->random_file;
        /*************/

		return 0;

    }

    int makeOutputDirectory() {
        
        /*******************************************/
        /* This function makes an output directory */
        /*******************************************/

        /* make an output directory, "output_dir" */
        char buf_dir[2048];
        struct stat st;
        sprintf(buf_dir, "%s", this->output_dir.c_str());
        if(ThisTask == 0) {
            if(stat(buf_dir, &st) != 0) {
            	int ret = mkdir(buf_dir, 0777);
            	if(ret == 0) {
            		printf("Directry \"%s\" is successfully made.\n ", buf_dir); 
            	} else {
            		printf("Directry \"%s\" fails to be made\n", buf_dir); 
                    return -1;
            	}
            }

            /* If there are multiple realizations, such as mock data, 
             * make an additional "output_dir/realization" directory under the "output_dir" directory.*/
            if(this->realization > 0) {
    
    			sprintf(buf_dir, "%s/%04d", this->output_dir.c_str(), this->realization);
    			if(stat(buf_dir, &st) != 0) {
    				int ret = mkdir(buf_dir, 0777);
    				if(ret == 0) {
    					if(ThisTask == 0) { printf("Directry \"%s\" is successfully made.\n ", buf_dir); }
    				} else {
    					if(ThisTask == 0) { printf("Directry \"%s\" fails to be made\n", buf_dir); }
                        return -1;
    				}
    			}
            }
        }

        #ifdef __MPI__
        MPI_Barrier( MPI_COMM_WORLD );
        #endif

        /* If there are multiple realizations, such as mock data, 
         * make an additional "output_dir/realization" directory under the "output_dir" directory.*/
        if(this->realization > 0) {

			sprintf(buf_dir, "%s/%04d", this->output_dir.c_str(), this->realization);
			if(stat(buf_dir, &st) != 0) {
				int ret = mkdir(buf_dir, 0777);
				if(ret == 0) {
					if(ThisTask == 0) { printf("Directry \"%s\" is successfully made.\n ", buf_dir); }
				} else {
					if(ThisTask == 0) { printf("Directry \"%s\" fails to be made\n", buf_dir); }
                    return -1;
				}
			}
            this->output_dir = buf_dir;
        }

		return 0;

	}

	int printParameters() {

        /***************************************************************************/
        /* This function outputs the read input parameters to the output directory */
        /***************************************************************************/

        FILE *fp;
        char buf[1024];
        sprintf(buf, "%s/input_parameters", this->output_dir.c_str());
        if( !(fp = fopen(buf,"w")) ) {
        	if(ThisTask == 0) { printf("output directry '%s' does not exist\n", this->output_dir.c_str()); }
        	return -1;
        }
    
        fprintf(fp, "\n");
        fprintf(fp, "data_dir = %s\n",  this->data_dir.c_str());
        fprintf(fp, "data_file = %s\n", this->data_file.c_str());
        fprintf(fp, "random_file = %s\n", this->random_file.c_str());
        fprintf(fp, "output_dir = %s\n", this->output_dir.c_str());
        fprintf(fp, "\n");
        fprintf(fp, "measure = %s\n", this->measure.c_str());
        fprintf(fp, "\n");
        fprintf(fp, "realization = %d\n", this->realization);
        fprintf(fp, "\n");
        fprintf(fp, "ell1 = %d\n", this->ell1);
        fprintf(fp, "ell2 = %d\n", this->ell2);
        fprintf(fp, "ELL = %d\n", this->ELL);
        fprintf(fp, "\n");
        fprintf(fp, "boxsize_x = %.5f\n", this->boxsize[0]);
        fprintf(fp, "boxsize_y = %.5f\n", this->boxsize[1]);
        fprintf(fp, "boxsize_z = %.5f\n", this->boxsize[2]);
        fprintf(fp, "\n");
        fprintf(fp, "n_mesh_x = %d\n", this->n_mesh[0]);
        fprintf(fp, "n_mesh_y = %d\n", this->n_mesh[1]);
        fprintf(fp, "n_mesh_z = %d\n", this->n_mesh[2]);
        fprintf(fp, "\n");
        fprintf(fp, "kmin = %.5f\n", this->kmin);
        fprintf(fp, "kmax = %.5f\n", this->kmax);
        fprintf(fp, "n_kbin = %d\n", this->n_kbin);
        fprintf(fp, "\n");
        fprintf(fp, "rmin = %.5f\n", this->rmin);
        fprintf(fp, "rmax = %.5f\n", this->rmax);
        fprintf(fp, "n_rbin = %d\n", this->n_rbin);
        fprintf(fp, "\n");
        fprintf(fp, "flag_recon = %s\n", this->flag_recon.c_str());
        fprintf(fp, "b1_fid = %.5f\n", this->b1_fid);
        fprintf(fp, "RG = %.5f\n", this->RG);
        fprintf(fp, "n_mesh_recon_x = %d\n", this->n_mesh_recon[0]);
        fprintf(fp, "n_mesh_recon_y = %d\n", this->n_mesh_recon[1]);
        fprintf(fp, "n_mesh_recon_z = %d\n", this->n_mesh_recon[2]);
        fprintf(fp, "\n");
        fprintf(fp, "ith_kbin = %d\n", this->ith_kbin);
        fprintf(fp, "ith_rbin = %d\n", this->ith_rbin);
        fprintf(fp, "\n");
        fprintf(fp, "sim_data = %s\n", this->sim_data.c_str());
        fprintf(fp, "flag_RSD = %s\n", this->flag_RSD.c_str());
        fprintf(fp, "log10_Mmin = %.5f\n", this->log10_Mmin);
        fprintf(fp, "log10_Mmax = %.5f\n", this->log10_Mmax);
        fprintf(fp, "\n");
        fprintf(fp, "NR = %d\n", this->NR);
        fprintf(fp, "data_mpi_file = %s\n", this->data_mpi_file.c_str());
        fprintf(fp, "\n");
        fprintf(fp, "assignment = %s\n", this->assign.c_str());
        fprintf(fp, "\n");
        fclose(fp);

		return 0;

	}

	int checkParameters() {

        if( this->data_file == "None" ) {
    		if(ThisTask == 0) {
				printf("The 'data_file' is not loaded correctly.\n");
                printf("It should NOT be 'None'.\n");
			}
			return -1;
        }

        if( this->random_file == "None" ) {
    		if(ThisTask == 0) {
				printf("The 'random_file' is not loaded correctly.\n");
                printf("It should NOT be 'None'.\n");
			}
			return -1;
        }

        if( this->output_dir == "None" ) {
    		if(ThisTask == 0) {
				printf("The 'output_dir' is not loaded correctly.\n");
                printf("It should NOT be 'None'.\n");
			}
			return -1;
        }

        if( this->data_dir == "None" ) {
    		if(ThisTask == 0) {
				printf("The 'data_dir' is not loaded correctly.\n");
                printf("It should NOT be 'None'.\n");
			}
			return -1;
        }

        if( !( (this->measure == "pk") || (this->measure == "bk") || (this->measure == "2PCF") || (this->measure == "3PCF") || (this->measure == "window2PCF") || (this->measure == "window3PCF") ) ) {
    		if(ThisTask == 0) {
				printf("The 'measure' is not loaded correctly.\n");
                printf("It should be 'pk', 'bk', '2PCF', '3PCF', 'window2PCF', or 'window3PCF'.\n");
			}
			return -1;
        }

		if( !( (this->boxsize[0] > 0.0) && (this->boxsize[1] > 0.0) && (this->boxsize[2] > 0.0) )  ) {
			if(ThisTask == 0) {
				printf("The 'boxsize' is not loaded correctly.\n");
                printf("It should be boxsize_x > 0, boxsize_y > 0, and boxsize_z > 0.\n");
			}
			return -1;
		}

		if( !( (this->n_mesh[0] > 0) && (this->n_mesh[1] > 0) && (this->n_mesh[2] > 0) )  ) {
			if(ThisTask == 0) {
				printf("The 'n_mesh' is not loaded correctly.\n");
                printf("It should be n_mesh_x > 0, n_mesh_y > 0, and n_mesh_z > 0.\n");
			}
			return -1;
		}

		if( !( (this->ell1 >= 0) && (this->ell2 >= 0) && (this->ELL >= 0) )  ) {
			if(ThisTask == 0) {
				printf("The 'ell1', 'ell2', or 'ELL' is not loaded correctly.\n");
                printf("It should be ell1 >= 0, ell2 >= 0, and ELL >= 0.\n");
			}
			return -1;
		}

        if( !(this->kmin > 0.0) ) {
			if(ThisTask == 0) {
				printf("The 'kmin' is not loaded correctly.\n");
                printf("It should be kmin > 0.0.\n");
            }
			return -1;
        }

        if( !(this->kmax > this->kmin) ) {
			if(ThisTask == 0) {
				printf("The 'kmax' is not loaded correctly.\n");
                printf("It should be kmax > kmin.\n");
            }
			return -1;
        }

        if( !(this->n_kbin >= 2) ) {
			if(ThisTask == 0) {
				printf("The 'n_kbin' is not loaded correctly.\n");
                printf("It should be n_kbin >= 2.\n");
            }
			return -1;
        }

         if( !(this->rmin > 0.0) ) {
			if(ThisTask == 0) {
				printf("The 'rmin' is not loaded correctly.\n");
                printf("It should be rmin > 0.0.\n");
            }
			return -1;
        }

        if( !(this->rmax > this->rmin) ) {
			if(ThisTask == 0) {
				printf("The 'rmax' is not loaded correctly.\n");
                printf("It should be rmax > rmin.\n");

            }
			return -1;
        }

        if( !(this->n_rbin >= 2) ) {
			if(ThisTask == 0) {
				printf("The 'n_rbin' is not loaded correctly.\n");
                printf("It should be n_rbin >= 2.\n");
            }
			return -1;
        }

		if( !( (this->assign == "NGP") || (this->assign == "CIC") || (this->assign == "TSC") )  ) {
			if(ThisTask == 0) {
				printf("The 'assign' is not loaded correctly.\n");
				printf("It should be 'NGP', 'CIC' or 'TSC'\n");
			}
			return -1;
		}

        if( !( (this->ith_kbin >= 0) && (this->ith_kbin < this->n_kbin) ) ) {
			if(ThisTask == 0) {
				printf("The 'ith_kbin' is not loaded correctly.\n");
                printf("It should be ith_kbin >= 0 and ith_kbin < n_kbin.\n");
            }
			return -1;
        }

        if( !( (this->ith_rbin >= 0) && (this->ith_rbin < this->n_rbin) ) ) {
			if(ThisTask == 0) {
				printf("The 'ith_rbin' is not loaded correctly.\n");
                printf("It should be ith_rbin >= 0 and ith_rbin < n_rbin.\n");
            }
			return -1;
        }

        if( !(this->NR >= 0) ) {
			if(ThisTask == 0) {
				printf("The 'NR' is not loaded correctly.\n");
                printf("It should be NR >= 0.\n");
            }
			return -1;
        }

        if( this->data_mpi_file == "None" ) {
    		if(ThisTask == 0) {
				printf("The 'data_mpi_file' is not loaded correctly.\n");
                printf("It should NOT be 'None'.\n");
			}
			return -1;
        }

        if( !(this->b1_fid > 0.0) ) {
			if(ThisTask == 0) {
				printf("The 'b1_fid' is not loaded correctly.\n");
                printf("It should be b1_fid > 0.0.\n");
            }
			return -1;
        }   


         if( !(this->RG >= 0.0) ) {
			if(ThisTask == 0) {
				printf("The 'RG' is not loaded correctly.\n");
                printf("It should be RG >= 0.0.\n");
            }
			return -1;
        }   

        if( !( (this->flag_recon == "True") || (this->flag_recon == "False") ) ) {
			if(ThisTask == 0) {
				printf("The 'flag_recon' is not loaded correctly.\n");
                printf("It should be flag_recon = True or False.\n");
            }
			return -1;
        } 

		if( !( (this->n_mesh_recon[0] > 0) && (this->n_mesh_recon[1] > 0) && (this->n_mesh_recon[2] > 0) )  ) {
			if(ThisTask == 0) {
				printf("The 'n_mesh_recon' is not loaded correctly.\n");
                printf("It should be n_mesh_recon_x > 0, n_mesh_recon_y > 0, and n_mesh_recon_z > 0.\n");
			}
			return -1;
		}

        if( !( this->realization >= 0 ) ) {
			if(ThisTask == 0) {
				printf("The 'realization' is not loaded correctly.\n");
                printf("It should be realization >= 0.\n");
            }
			return -1;
        }  

        if( !( (this->flag_RSD == "True") || (this->flag_RSD == "False") ) ) {
			if(ThisTask == 0) {
				printf("The 'flag_RSD' is not loaded correctly.\n");
                printf("It should be flag_RSD = True or False.\n");
            }
			return -1;
        }   

        if( !( (this->sim_data == "Gadget") || (this->sim_data == "Rockstar") ) ) {
			if(ThisTask == 0) {
				printf("The 'sim_data' is not loaded correctly.\n");
                printf("It should be sim_data = Gadget or Rockstar.\n");
            }
			return -1;
        }   



        if( !(this->log10_Mmin > 0.0) ) {
			if(ThisTask == 0) {
				printf("The 'log10_Mmin' is not loaded correctly.\n");
                printf("It should be log10_Mmin > 0.0.\n");
            }
			return -1;
        }   
	
        if( !(this->log10_Mmin < this->log10_Mmax) ) {
			if(ThisTask == 0) {
				printf("The 'log10_Mmax' is not loaded correctly.\n");
                printf("It should be log10_Mmax > log10_Mmin.\n");
            }
			return -1;
        }   
		
        return 0;
	
	}

};




#endif

