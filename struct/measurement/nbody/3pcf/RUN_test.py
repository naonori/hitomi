#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess

try:
    os.mkdir("LOG")
except:
    print("")

try:
    os.mkdir("PARAMS")
except:
    print("")

realization = 1

measure = "3PCF"

boxsize_x = 1000.0
boxsize_y = 1000.0
boxsize_z = 1000.0

n_mesh_x = 256 
n_mesh_y = 256 
n_mesh_z = 256 

n_kbin = 10
kmin = 0.02
kmax = 0.20

n_rbin = 13
rmin = 30.0
rmax = 150.0

flag_recon = "False" # "True" or "False"
b1_fid = 2.0
RG = 15.0

ith_kbin = 0
ith_rbin = 0

# b1_fid_dict = {}
# b1_fid_dict.update({"zbin0_Mmin12.5_Mmax13.0": })
# b1_fid_dict.update({"zbin0_Mmin13.0_Mmax13.5": })
# b1_fid_dict.update({"zbin0_Mmin13.5_Mmax14.0": })
# b1_fid_dict.update({"zbin1_Mmin12.5_Mmax13.0": })
# b1_fid_dict.update({"zbin1_Mmin13.0_Mmax13.5": })
# b1_fid_dict.update({"zbin1_Mmin13.5_Mmax14.0": })
# b1_fid_dict.update({"zbin2_Mmin12.5_Mmax13.0": })
# b1_fid_dict.update({"zbin2_Mmin13.0_Mmax13.5": })
# b1_fid_dict.update({"zbin2_Mmin13.5_Mmax14.0": })
# b1_fid_dict.update({"zbin3_Mmin12.5_Mmax13.0": })
# b1_fid_dict.update({"zbin3_Mmin13.0_Mmax13.5": })
# b1_fid_dict.update({"zbin3_Mmin13.5_Mmax14.0": })
# b1_fid_dict.update({"zbin4_Mmin12.5_Mmax13.0": })
# b1_fid_dict.update({"zbin4_Mmin13.0_Mmax13.5": })
# b1_fid_dict.update({"zbin4_Mmin13.5_Mmax14.0": })

for ith_bin in range(n_rbin):

    if measure == "bk":
        ith_kbin = ith_bin
    elif measure == "3PCF":
        ith_rbin = ith_bin
    else:
        print("ERROR")
        exit()

    for zbin in [0, 1, 2, 3, 4]:
        for flag_RSD in ["False", "True"]:
            
            if flag_RSD == "False":
                MULTIPOLE = [(0,0,0), (1,1,0), (2,2,0)]
            else:
                MULTIPOLE = [(0,0,0), (1,1,0), (2,2,0), (2,0,2), (1,1,2), (0,2,2)]
    
            for (ell1, ell2, ELL) in MULTIPOLE:
    
                for sim_data in ["Gadget", "Rockstar"]:
                    
                    if sim_data == "Gadget":
                        Mmin_Mmax = [(12.9, 13.1)]
                    else:
                        Mmin_Mmax = [(12.5, 13.0), (13.0, 13.5), (13.5, 14.0)]
                    
                    for (log10_Mmin, log10_Mmax) in Mmin_Mmax:
    
                        if flag_recon == "True":
                            b1_fid = b1_fid_dict["zbin%d_Mmin%2.1f_Mmax%2.1f" % (zbin, log10_Mmin, log10_Mmax)]
                        
                        data_dir = "/mwork0/sugiymnn/WORK/data/nbody/%04d" % realization
                        
                        if sim_data == "Gadget":
                            data_file = "snapdir_%03d/snapshot_%03d" % (zbin, zbin)
                            random_file = "AAA"
                        elif sim_data == "Rockstar":
                            data_file = "rockstar/out_%d.list" % (zbin)
                            if flag_recon == "False":
                                random_file = "AAA"
                            elif flag_recon == "True":
                                random_file = "rockstar/random/random_Mmin%s_Mmax%s_zbin%d.dat" % (log10_Mmin, log10_Mmax, zbin)
                        else:
                            print("ERROR")
                            exit()
                        
                        if sim_data == "Gadget":
                            if flag_recon == "False":
                                output_dir = "results_test_%s_zbin%d_RSD%s" % (sim_data, zbin, flag_RSD)
                            elif flag_recon == "True":
                                output_dir = "results_test_%s_zbin%d_RSD%s_recon_R%02d" % (sim_data, zbin, flag_RSD, RG)
                            else:
                                print("ERROR")
                                exit()
                        
                        elif sim_data == "Rockstar":
                            if flag_recon == "False":
                                output_dir = "results_test_%s_Mmin%2.1f_Mmax%2.1f_zbin%d_RSD%s" % (sim_data, log10_Mmin, log10_Mmax, zbin, flag_RSD)
                            elif flag_recon == "True":
                                output_dir = "results_test_%s_Mmin%2.1f_Mmax%2.1f_zbin%d_RSD%s_recon_R%02d" % (sim_data, log10_Mmin, log10_Mmax, zbin, flag_RSD, RG)
                            else:
                                print("ERROR")
                                exit()
                        else:
                            print("ERROR")
                            exit()
                        
                        
                        fr = open("default_param.ini", "r")
                        AA = fr.readlines()
                        fr.close()
                        
                        AA = [AA[i].replace("data_dir = /mwork0/sugiymnn/WORK/data/boss/galaxy_DR12v5_CMASSLOWZTOT",\
                                            "data_dir = %s" % data_dir) for i in range(len(AA))]
                        
                        AA = [AA[i].replace("data_file   = galaxy_DR12v5_CMASSLOWZTOT_North_ZBIN1.dat",\
                                            "data_file   = %s" % data_file) for i in range(len(AA))]
    
                        AA = [AA[i].replace("random_file = random_DR12v5_CMASSLOWZTOT_North_ZBIN1.dat",\
                                            "random_file = %s" % random_file) for i in range(len(AA))]
    
                        AA = [AA[i].replace("output_dir = results",\
                                            "output_dir = %s" % output_dir) for i in range(len(AA))]
                        
                        AA = [AA[i].replace("measure = pk", "measure = %s" % measure) for i in range(len(AA))]
                        
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
                        
                        AA = [AA[i].replace("kmin = 0.01", "kmin = %1.3f" % kmin) for i in range(len(AA))]
                        AA = [AA[i].replace("kmax = 0.2",  "kmax = %1.3f" % kmax) for i in range(len(AA))]
                        AA = [AA[i].replace("n_kbin = 20", "n_kbin = %02d" % n_kbin) for i in range(len(AA))]
                        
                        AA = [AA[i].replace("rmin = 30", "rmin = %2.1f" % rmin) for i in range(len(AA))]
                        AA = [AA[i].replace("rmax = 150.0", "rmax = %3.1f" % rmax) for i in range(len(AA))]
                        AA = [AA[i].replace("n_rbin = 25", "n_rbin = %02d" % n_rbin) for i in range(len(AA))]
                        
                        AA = [AA[i].replace("flag_recon = False", "flag_recon = %s" % flag_recon) for i in range(len(AA))]
                        AA = [AA[i].replace("b1_fid = 1.0", "b1_fid = %1.3f" % b1_fid) for i in range(len(AA))]
                        AA = [AA[i].replace("RG = 15.0", "RG = %2.1f" % RG) for i in range(len(AA))]
                        AA = [AA[i].replace("n_mesh_recon_x = 512", "n_mesh_recon_x = %3.0d" % n_mesh_x) for i in range(len(AA))]
                        AA = [AA[i].replace("n_mesh_recon_y = 512", "n_mesh_recon_y = %3.0d" % n_mesh_y) for i in range(len(AA))]
                        AA = [AA[i].replace("n_mesh_recon_z = 512", "n_mesh_recon_z = %3.0d" % n_mesh_z) for i in range(len(AA))]
                        
                        AA = [AA[i].replace("ith_kbin = 0", "ith_kbin = %d" % ith_kbin) for i in range(len(AA))]
                        AA = [AA[i].replace("ith_rbin = 0", "ith_rbin = %d" % ith_rbin) for i in range(len(AA))]
                        
                        AA = [AA[i].replace("sim_data = Gadget", "sim_data = %s" % sim_data) for i in range(len(AA))]
                        AA = [AA[i].replace("flag_RSD = True", "flag_RSD = %s" % flag_RSD) for i in range(len(AA))]
                        AA = [AA[i].replace("log10_Mmin = 12.9", "log10_Mmin = %2.1f" % log10_Mmin) for i in range(len(AA))]
                        AA = [AA[i].replace("log10_Mmax = 13.5", "log10_Mmax = %2.1f" % log10_Mmax) for i in range(len(AA))]
                        
                        if sim_data == "Gadget":
                            if flag_recon == "False":
                                fname = "%s%d%d%d_%s_zbin%d_RSD%s_%02d" % (measure, ell1, ell2, ELL, sim_data, zbin, flag_RSD, ith_bin)
                            elif flag_recon == "True":
                                fname = "%s%d%d%d_%s_zbin%d_RSD%s_recon_R%02d_%02d" % (measure, ell1, ell2, ELL, sim_data, zbin, flag_RSD, RG, ith_bin)
                        
                        elif sim_data == "Rockstar":
                            if flag_recon == "False":
                                fname = "%s%d%d%d_%s_Mmin%2.1f_Mmax%2.1f_zbin%d_RSD%s_%02d" % (measure, ell1, ell2, ELL, sim_data, log10_Mmin, log10_Mmax, zbin, flag_RSD, ith_bin)
                            elif flag_recon == "True":
                                fname = "%s%d%d%d_%s_Mmin%2.1f_Mmax%2.1f_zbin%d_RSD%s_recon_R%02d_%02d" % (measure, ell1, ell2, ELL, sim_data, log10_Mmin, log10_Mmax, zbin, flag_RSD, RG, ith_bin)
                        
                        else:
                            print("ERROR")
                            exit()
                        
                        
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
                                
                           
