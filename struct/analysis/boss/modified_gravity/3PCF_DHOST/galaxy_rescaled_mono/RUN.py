#######################
# Purpose of this code
#######################
# This code executes MontePython.py. This code runs the case -N =200000 in "chains_base" and computes the best fit and covariance matrix of the parameters. Then it uses them to perform the parameter estimation again.
#######################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
import scipy as sp
import subprocess

output_dir = "chains_base" #
#output_dir = "chains_main" # chains_base or chains_main
try:
    os.mkdir(output_dir)
except:
    pass

MOCK = []
for i in range(100):
    MOCK.append("mock_%04d" % (i+1))

#for param in ["DHOST_North_zbin1", "DHOST_South_zbin1", "DHOST_North_zbin3", "DHOST_South_zbin3",\
#              "DHOST_zbin1", "DHOST_zbin3", "DHOST_North", "DHOST_South", "DHOST"]:
for param in ["DHOST"]:

    try:
        os.mkdir("%s/%s" % (output_dir, param))
    except:
        pass
    
#    for data_type in ["galaxy"] + MOCK:
    for data_type in ["galaxy"]:
    
        print("%s_%s" % (data_type, param))

        fr = open("run_base.sh", "r")
        AA = fr.readlines()
        fr.close()
        
        if output_dir == "chains_base":
            AA[8] = "python MontePython.py run -o %s/%s/%s_%s -p input/%s_%s.param --conf=my.conf --superupdate 20 -N 200000" % (output_dir, param, data_type, param, data_type, param)
    
        elif output_dir == "chains_main":
            AA[8] = "python MontePython.py run -o %s/%s/%s_%s -p input/%s_%s.param -b chains_base/%s/%s_%s/%s_%s.bestfit -c chains_base/%s/%s_%s/%s_%s.covmat --conf=my.conf --superupdate 20 -N 200000" % (output_dir, param, data_type, param, data_type, param, param, data_type, param, data_type, param, param, data_type, param, data_type, param)
        
        fw = open("run_new.sh", "w")
        fw.writelines(AA)
        fw.close()
        
        subprocess.call(["chmod", "u+x", "run_new.sh"])

        if output_dir == "chains_main" and data_type == "galaxy":
#            for i in range(8):
            for i in range(1):
#            for i in range(7):
                subprocess.call(["qsub", "run_new.sh"])

        else:
            subprocess.call(["qsub", "run_new.sh"])
        
