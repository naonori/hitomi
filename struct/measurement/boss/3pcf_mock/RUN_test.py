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

output_dir = "results_test"
measure = "3PCF"

ell1 = 0
ell2 = 0
ELL  = 0

boxsize_x = 1350.0
boxsize_y = 2450.0
boxsize_z = 1400.0
n_mesh_x = 250
n_mesh_y = 460
n_mesh_z = 260

n_rbin = 13
rmin = 30.0
rmax = 150.0

flag_recon = "False" # "True" or "False"
b1_fid = 2.0
RG = 15.0

realization = 1

for ith_bin in range(0,13):

    fr = open("default_param.ini", "r")
    AA = fr.readlines()
    fr.close()

    AA = [AA[i].replace("data_file   = Patchy-Mocks-DR12NGC-COMPSAM_V6C_ZBIN1_0001.dat", "data_file   = Patchy-Mocks-DR12NGC-COMPSAM_V6C_ZBIN1_%04d.dat" % realization) for i in range(len(AA))]
    AA = [AA[i].replace("realization = 0", "realization = %d" % realization) for i in range(len(AA))]

    AA = [AA[i].replace("output_dir = results", "output_dir = %s" % output_dir) for i in range(len(AA))]
    AA = [AA[i].replace("measure = pk", "measure = %s" % measure) for i in range(len(AA))]
    
    AA = [AA[i].replace("ell1 = 0", "ell1 = %d" % ell1) for i in range(len(AA))]
    AA = [AA[i].replace("ell2 = 0", "ell2 = %d" % ell2) for i in range(len(AA))]
    AA = [AA[i].replace("ELL  = 0", "ELL  = %d" %  ELL) for i in range(len(AA))]
     
    AA = [AA[i].replace("boxsize_x = 1000.0", "boxsize_x = %3.1f" % boxsize_x) for i in range(len(AA))]
    AA = [AA[i].replace("boxsize_y = 1000.0", "boxsize_y = %3.1f" % boxsize_y) for i in range(len(AA))]
    AA = [AA[i].replace("boxsize_z = 1000.0", "boxsize_z = %3.1f" % boxsize_z) for i in range(len(AA))]
    AA = [AA[i].replace("n_mesh_x = 512", "n_mesh_x = %3.0d" % n_mesh_x) for i in range(len(AA))]
    AA = [AA[i].replace("n_mesh_y = 512", "n_mesh_y = %3.0d" % n_mesh_y) for i in range(len(AA))]
    AA = [AA[i].replace("n_mesh_z = 512", "n_mesh_z = %3.0d" % n_mesh_z) for i in range(len(AA))]
    
    AA = [AA[i].replace("rmin = 30", "rmin = %f" % rmin) for i in range(len(AA))]
    AA = [AA[i].replace("rmax = 150.0", "rmax = %f" % rmax) for i in range(len(AA))]
    AA = [AA[i].replace("n_rbin = 25", "n_rbin = %d" % n_rbin) for i in range(len(AA))]

    AA = [AA[i].replace("ith_rbin = 0", "ith_rbin = %d" % ith_bin) for i in range(len(AA))]

    AA = [AA[i].replace("flag_recon = False", "flag_recon = %s" % flag_recon) for i in range(len(AA))]
    AA = [AA[i].replace("b1_fid = 1.0", "b1_fid = %f" % b1_fid) for i in range(len(AA))]
    AA = [AA[i].replace("RG = 15.0", "RG = %f" % RG) for i in range(len(AA))]
    AA = [AA[i].replace("n_mesh_recon_x = 512", "n_mesh_recon_x = %3.0d" % n_mesh_x) for i in range(len(AA))]
    AA = [AA[i].replace("n_mesh_recon_y = 512", "n_mesh_recon_y = %3.0d" % n_mesh_y) for i in range(len(AA))]
    AA = [AA[i].replace("n_mesh_recon_z = 512", "n_mesh_recon_z = %3.0d" % n_mesh_z) for i in range(len(AA))]

    fname = "zeta%d%d%d_%02d_%04d" % (ell1, ell2, ELL, ith_bin, realization)
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
            
       
