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
measure = "bk"

ell1 = 0
ell2 = 0
ELL  = 0

boxsize_x = 1350.0
boxsize_y = 2450.0
boxsize_z = 1400.0
n_mesh_x = 250
n_mesh_y = 460
n_mesh_z = 260

n_kbin = 10
kmin = 0.02
kmax = 0.20

flag_recon = "False" # "True" or "False"
b1_fid = 2.0
RG = 15.0

for ith_bin in range(0,10):

    fr = open("default_param.ini", "r")
    AA = fr.readlines()
    fr.close()

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
    
    AA = [AA[i].replace("kmin = 0.01", "kmin = %f" % kmin) for i in range(len(AA))]
    AA = [AA[i].replace("kmax = 0.2",  "kmax = %f" % kmax) for i in range(len(AA))]
    AA = [AA[i].replace("n_kbin = 20", "n_kbin = %d" % n_kbin) for i in range(len(AA))]

    AA = [AA[i].replace("ith_kbin = 0", "ith_kbin = %d" % ith_bin) for i in range(len(AA))]

    AA = [AA[i].replace("flag_recon = False", "flag_recon = %s" % flag_recon) for i in range(len(AA))]
    AA = [AA[i].replace("b1_fid = 1.0", "b1_fid = %f" % b1_fid) for i in range(len(AA))]
    AA = [AA[i].replace("RG = 15.0", "RG = %f" % RG) for i in range(len(AA))]
    AA = [AA[i].replace("n_mesh_recon_x = 512", "n_mesh_recon_x = %3.0d" % n_mesh_x) for i in range(len(AA))]
    AA = [AA[i].replace("n_mesh_recon_y = 512", "n_mesh_recon_y = %3.0d" % n_mesh_y) for i in range(len(AA))]
    AA = [AA[i].replace("n_mesh_recon_z = 512", "n_mesh_recon_z = %3.0d" % n_mesh_z) for i in range(len(AA))]
    
    fname = "bk%d%d%d_%02d_recon%s_R%02d" % (ell1, ell2, ELL, ith_bin, flag_recon, RG)
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
            
       
