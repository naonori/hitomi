#################
## This code generates various input parameter files in PARAM/ and uses them to submit jobs to a PC cluster. If uses want to generate only the parameter files, comment out the following line 287.
#################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess

data = "galaxy" # "galaxy" or "mock"
weight = 0 # 0,1,2,3
measure = "pk" # pk, 2PCF, window2PCF, bk, 3PCF, or window3PCF

flag_recon = "False" # "True" or "False"
RG = 15.0
b1_fid = 2.0

rmin = 30.0
rmax = 150.0
n_rbin = 25

kmin = 0.01 
kmax = 0.2
n_kbin = 20 

try:
    os.mkdir("LOG")
except:
    print("")

try:
    os.mkdir("PARAMS")
except:
    print("")

MULTIPOLE = []
MONO = [(0,0,0),(1,1,0),(2,2,0),(3,3,0),(4,4,0)]
QUAD = [(2,0,2),(1,1,2),(0,2,2),(3,1,2),(2,2,2),(1,3,2),(4,2,2),(3,3,2),(2,4,2)]
MULTIPOLE.extend(MONO)
MULTIPOLE.extend(QUAD)

for NS in ["North", "South"]:
    for zbin in [1,3]:
        for (ell1, ell2, ELL) in [(0,0,0)]:
#        for (ell1, ell2, ELL) in [(0,0,0),(1,1,0),(2,0,2),(1,1,2)]:
#        for (ell1, ell2, ELL) in MULTIPOLE:
            for ith_bin in [0]: # for ith_bin in range(0,13):
                for realization in [0]:

                    if (data == "mock") and (measure == "window2PCF" or measure == "window3PCF"):

                        realization = 0

                    elif (data == "mock") and not (measure == "window2PCF" or measure == "window3PCF"):

                        if realization == 0:

                            print("ERROR")
                            exit()

                        else:
                            
                            pass

                    elif data == "galaxy":
    
                        if realization != 0:

                            print("ERROR")
                            exit()

                        else:
                            
                            pass

                    for NR in [0]:
    
                        fr = open("default_param.ini", "r")
                        AA = fr.readlines()
                        fr.close()
            
                        if NS == "North" and zbin == 1:
            
                            boxsize_x = 1350.0
                            boxsize_y = 2450.0
                            boxsize_z = 1400.0
                            
                            n_mesh_x = 250
                            n_mesh_y = 460
                            n_mesh_z = 260
            	
                        elif NS == "North" and zbin == 3:
                            
                            boxsize_x = 1800.0
                            boxsize_y = 3400.0
                            boxsize_z = 1900.0
                            
                            n_mesh_x = 340
                            n_mesh_y = 650
                            n_mesh_z = 360
                                
                        elif NS == "South" and zbin == 1:
                            
                            boxsize_x = 1000.0
                            boxsize_y = 1900.0
                            boxsize_z = 1100.0
                            
                            n_mesh_x = 190
                            n_mesh_y = 360
                            n_mesh_z = 210
                                
                        elif NS == "South" and zbin == 3:
                            
                            boxsize_x = 1000.0
                            boxsize_y = 2600.0
                            boxsize_z = 1500.0
                            
                            n_mesh_x = 190
                            n_mesh_y = 500
                            n_mesh_z = 280
            	
                        else:
                            print("ERROR")
                            sys.exit()
            
                        if NS == "North":
                            NGC_SGC = "NGC"
                        elif NS == "South":
                            NGC_SGC = "SGC"
                        else:
                            print("ERROR")
                            sys.exit()
            
                        if data == "galaxy":
                        
                            AA = [AA[i].replace("data_dir = /mwork0/sugiymnn/WORK/data/boss/galaxy_DR12v5_CMASSLOWZTOT",\
                                                "data_dir = /mwork0/sugiymnn/WORK/data/boss/galaxy_DR12v5_CMASSLOWZTOT") for i in range(len(AA))]
        
                            if weight == 0:
                                AA = [AA[i].replace("galaxy_DR12v5_CMASSLOWZTOT_North_ZBIN1.dat",\
                                                    "galaxy_DR12v5_CMASSLOWZTOT_%s_ZBIN%d.dat" % (NS, zbin)) for i in range(len(AA))]
                            elif weight == 1:
                                AA = [AA[i].replace("galaxy_DR12v5_CMASSLOWZTOT_North_ZBIN1.dat",\
                                                    "galaxy_DR12v5_CMASSLOWZTOT_%s_ZBIN%d_Weight1_OnlySys.dat" % (NS, zbin)) for i in range(len(AA))]
                            elif weight == 2:
                                AA = [AA[i].replace("galaxy_DR12v5_CMASSLOWZTOT_North_ZBIN1.dat",\
                                                    "galaxy_DR12v5_CMASSLOWZTOT_%s_ZBIN%d_Weight2_NoSys.dat" % (NS, zbin)) for i in range(len(AA))]
                            elif weight == 3:
                                AA = [AA[i].replace("galaxy_DR12v5_CMASSLOWZTOT_North_ZBIN1.dat",\
                                                    "galaxy_DR12v5_CMASSLOWZTOT_%s_ZBIN%d_Weight3_NoWeight.dat" % (NS, zbin)) for i in range(len(AA))]
        
                            AA = [AA[i].replace("random_DR12v5_CMASSLOWZTOT_North_ZBIN1.dat",\
                                                "random_DR12v5_CMASSLOWZTOT_%s_ZBIN%d.dat" % (NS, zbin)) for i in range(len(AA))]
                
                            if flag_recon == "False":
                                
                                if weight == 0 and not (measure == "window2PCF" or measure == "window3PCF"):
                                    AA = [AA[i].replace("output_dir = results",\
                                                        "output_dir = galaxy_%s_zbin%d" %(NS, zbin)) for i in range(len(AA))]
    
                                elif weight == 0 and (measure == "window2PCF" or measure == "window3PCF"):
       
                                    AA = [AA[i].replace("output_dir = results",\
                                                        "output_dir = galaxy_%s_zbin%d_Window" %(NS, zbin)) for i in range(len(AA))]
                                elif weight == 1:
                                    AA = [AA[i].replace("output_dir = results",\
                                                        "output_dir = galaxy_%s_zbin%d_Weight1_OnlySys" %(NS, zbin)) for i in range(len(AA))]
                                elif weight == 2:
                                    AA = [AA[i].replace("output_dir = results",\
                                                        "output_dir = galaxy_%s_zbin%d_Weight2_NoSys" %(NS, zbin)) for i in range(len(AA))]
                                elif weight == 3:
                                    AA = [AA[i].replace("output_dir = results",\
                                                        "output_dir = galaxy_%s_zbin%d_Weight3_NoWeight" %(NS, zbin)) for i in range(len(AA))]
        
                            elif flag_recon == "True":
    
                                if weight == 0 and not (measure == "window2PCF" or measure == "window3PCF"):
                                    AA = [AA[i].replace("output_dir = results",\
                                                        "output_dir = galaxy_%s_zbin%d_recon_R%02d" %(NS, zbin, RG)) for i in range(len(AA))]
    
                                elif weight == 0 and (measure == "window2PCF" or measure == "window3PCF"):
 
                                    AA = [AA[i].replace("output_dir = results",\
                                                        "output_dir = galaxy_%s_zbin%d_recon_R%02d_Window" %(NS, zbin, RG)) for i in range(len(AA))]
    
                                elif weight == 1:
                                    AA = [AA[i].replace("output_dir = results",\
                                                        "output_dir = galaxy_%s_zbin%d_recon_R%02d_Weight1_OnlySys" %(NS, zbin, RG)) for i in range(len(AA))]
                                elif weight == 2:
                                    AA = [AA[i].replace("output_dir = results",\
                                                        "output_dir = galaxy_%s_zbin%d_recon_R%02d_Weight2_NoSys" %(NS, zbin, RG)) for i in range(len(AA))]
                                elif weight == 3:
                                    AA = [AA[i].replace("output_dir = results",\
                                                        "output_dir = galaxy_%s_zbin%d_recon_R%02d_Weight3_NoWeight" %(NS, zbin, RG)) for i in range(len(AA))]
            
                        elif data == "mock":
            
                            AA = [AA[i].replace("data_dir = /mwork0/sugiymnn/WORK/data/boss/galaxy_DR12v5_CMASSLOWZTOT",\
                                                "data_dir = /mwork0/sugiymnn/WORK/data/boss/Patchy-Mocks-DR12%s-COMPSAM_V6C_ZBIN%d" % (NGC_SGC, zbin)) for i in range(len(AA))]
                
                            AA = [AA[i].replace("galaxy_DR12v5_CMASSLOWZTOT_North_ZBIN1.dat",\
                                                "Patchy-Mocks-DR12%s-COMPSAM_V6C_ZBIN%d_%04d.dat" % (NGC_SGC, zbin, realization)) for i in range(len(AA))]
                            AA = [AA[i].replace("random_DR12v5_CMASSLOWZTOT_North_ZBIN1.dat",\
                                                "Patchy-Mocks-Randoms-DR12%s-COMPSAM_V6C_x100_ZBIN%d.dat" % (NGC_SGC, zbin)) for i in range(len(AA))]
        
                            if flag_recon == "False" and not (measure == "window2PCF" or measure == "window3PCF"):

                                AA = [AA[i].replace("output_dir = results",\
                                                    "output_dir = mock_%s_zbin%d" %(NS, zbin)) for i in range(len(AA))]
    
    
                            elif flag_recon == "False" and (measure == "window2PCF" or measure == "window3PCF"):

                                AA = [AA[i].replace("output_dir = results",\
                                                    "output_dir = mock_%s_zbin%d_Window" %(NS, zbin)) for i in range(len(AA))]
    
                            elif flag_recon == "True" and not (measure == "window2PCF" or measure == "window3PCF"):
                                AA = [AA[i].replace("output_dir = results",\
                                                    "output_dir = mock_%s_zbin%d_recon_R%02d" %(NS, zbin, RG)) for i in range(len(AA))]
                            
                            elif flag_recon == "True" and (measure == "window2PCF" or measure == "window3PCF"):
    
                                AA = [AA[i].replace("output_dir = results",\
                                                    "output_dir = mock_%s_zbin%d_recon_R%02d_Window" %(NS, zbin, RG)) for i in range(len(AA))]
    
                        AA = [AA[i].replace("measure = pk", "measure = %s" % measure) for i in range(len(AA))]
    
                        if data == "galaxy":

                            AA = [AA[i].replace("realization = 0", "realization = 0") for i in range(len(AA))]

                        elif data == "mock":
                            
                            if measure == "window2PCF" or measure == "window3PCF":

                                AA = [AA[i].replace("realization = 0", "realization = 0") for i in range(len(AA))]

                            else:
                                
                                AA = [AA[i].replace("realization = 0", "realization = %d" % realization) for i in range(len(AA))]
                        
                        AA = [AA[i].replace("ell1 = 0", "ell1 = %d" % ell1) for i in range(len(AA))]
                        AA = [AA[i].replace("ell2 = 0", "ell2 = %d" % ell2) for i in range(len(AA))]
                        AA = [AA[i].replace("ELL  = 0", "ELL  = %d" %  ELL) for i in range(len(AA))]
                        
                        AA = [AA[i].replace("boxsize_x = 1000.0", "boxsize_x = %3.1f" % boxsize_x) for i in range(len(AA))]
                        AA = [AA[i].replace("boxsize_y = 1000.0", "boxsize_y = %3.1f" % boxsize_y) for i in range(len(AA))]
                        AA = [AA[i].replace("boxsize_z = 1000.0", "boxsize_z = %3.1f" % boxsize_z) for i in range(len(AA))]
                        
                        AA = [AA[i].replace("n_mesh_x = 512", "n_mesh_x = %3.0d" % n_mesh_x) for i in range(len(AA))]
                        AA = [AA[i].replace("n_mesh_y = 512", "n_mesh_y = %3.0d" % n_mesh_y) for i in range(len(AA))]
                        AA = [AA[i].replace("n_mesh_z = 512", "n_mesh_z = %3.0d" % n_mesh_z) for i in range(len(AA))]

                        AA = [AA[i].replace("kmin = 0.01", "kmin = %f" % kmin) for i in range(len(AA))]
                        AA = [AA[i].replace("kmax = 0.2", "kmax = %f" % kmax) for i in range(len(AA))]
                        AA = [AA[i].replace("n_kbin = 20", "n_kbin = %d" % n_kbin) for i in range(len(AA))]

                        AA = [AA[i].replace("rmin = 30", "rmin = %f" % rmin) for i in range(len(AA))]
                        AA = [AA[i].replace("rmax = 150.0", "rmax = %f" % rmax) for i in range(len(AA))]
                        AA = [AA[i].replace("n_rbin = 25", "n_rbin = %d" % n_rbin) for i in range(len(AA))]



                        AA = [AA[i].replace("flag_recon = False", "flag_recon = %s" % flag_recon) for i in range(len(AA))]
                        AA = [AA[i].replace("b1_fid = 1.0", "b1_fid = %1.1f" % b1_fid) for i in range(len(AA))]
                        AA = [AA[i].replace("RG = 15.0", "RG = %2.1f" % RG) for i in range(len(AA))]
                        
                        AA = [AA[i].replace("n_mesh_recon_x = 512", "n_mesh_recon_x = %3.0d" % n_mesh_x) for i in range(len(AA))]
                        AA = [AA[i].replace("n_mesh_recon_y = 512", "n_mesh_recon_y = %3.0d" % n_mesh_y) for i in range(len(AA))]
                        AA = [AA[i].replace("n_mesh_recon_z = 512", "n_mesh_recon_z = %3.0d" % n_mesh_z) for i in range(len(AA))]
     
                        if measure == "pk" or measure == "bk":
                            AA = [AA[i].replace("ith_kbin = 0", "ith_kbin = %d" % ith_bin) for i in range(len(AA))]
                        elif measure == "2PCF" or measure == "3PCF" or measure == "window2PCF" or measure == "window3PCF": 
                            AA = [AA[i].replace("ith_rbin = 0", "ith_rbin = %d" % ith_bin) for i in range(len(AA))]
    
                        AA = [AA[i].replace("NR = 0", "NR = %d" % NR) for i in range(len(AA))]
                        AA = [AA[i].replace("data_mpi_file = Patchy-Mocks-DR12NGC-COMPSAM_V6C_ZBIN1", "data_mpi_file = Patchy-Mocks-DR12%s-COMPSAM_V6C_ZBIN%d" % (NGC_SGC, zbin)) for i in range(len(AA))]
       
                        fname = "%s_%s_zbin%d_l1l2L%d%d%d_weight%d_bin%02d_recon%s_R%04d_NR%02d_%s" % (data, NS, zbin, ell1, ell2, ELL, weight, ith_bin, flag_recon, realization, NR, measure)
                        param_file = "PARAMS/param_%s.ini" % (fname)
                        fw = open(param_file, "w")
                        fw.writelines(AA)
                        fw.close()
            
                        fr = open("run_base.sh", "r")
                        AA = fr.readlines()
                        fr.close()
            	
                        log_file = "LOG/%s.log" % (fname)
                        AA = [AA[i].replace("./a.out default_param.ini > log", "./a.out %s > %s" % (param_file, log_file)) for i in range(len(AA))]
            
                        fw = open("run_new.sh", "w")
                        fw.writelines(AA)
                        fw.close()
            
                        subprocess.call(["chmod", "u+x", "run_new.sh"])
                        subprocess.call(["qsub", "run_new.sh"])
            
           
