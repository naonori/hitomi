#################
## This code submit jobs to a PC cluster with various input parameters.
#################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess

data = "galaxy" # galaxy or mock

zbin = 1 # 1 or 3

MULTIPOLE = []
MONO = [(0,0,0),(1,1,0),(2,2,0),(3,3,0),(4,4,0)]
QUAD = [(2,0,2),(1,1,2),(0,2,2),(3,1,2),(2,2,2),(1,3,2),(4,2,2),(3,3,2),(2,4,2)]
MULTIPOLE.extend(MONO)
MULTIPOLE.extend(QUAD)

PARAM_NAME = []
PARAM_NAME.append("FG_b3_f0")
PARAM_NAME.append("FS_b3_f0")
PARAM_NAME.append("FT_b3_f0")
PARAM_NAME.append("FG_b2_f1")
PARAM_NAME.append("FS_b2_f1")
PARAM_NAME.append("FT_b2_f1")
PARAM_NAME.append("FG_b1_f2")
PARAM_NAME.append("FS_b1_f2")
PARAM_NAME.append("FT_b1_f2")

PARAM_NAME.append("GG_b2_f1")
PARAM_NAME.append("GS_b2_f1")
PARAM_NAME.append("GT_b2_f1")
PARAM_NAME.append("GG_b1_f2")
PARAM_NAME.append("GS_b1_f2")
PARAM_NAME.append("GT_b1_f2")
PARAM_NAME.append("GG_b0_f3")
PARAM_NAME.append("GS_b0_f3")
PARAM_NAME.append("GT_b0_f3")

PARAM_NAME.append("b3_f1")
PARAM_NAME.append("b2_f2")
PARAM_NAME.append("b1_f3")
PARAM_NAME.append("b0_f4")

for (ell1, ell2, ELL) in MULTIPOLE:
    for param_name in PARAM_NAME:

        fr = open("run_base.sh", "r")
        AA = fr.readlines()
        fr.close()
        
        AA = [AA[i].replace("python calc_3pcf.py", "python calc_3pcf.py %s %d %d %d %d %s" % (data, zbin, ell1, ell2, ELL, param_name)) for i in range(len(AA))]
        
        fw = open("run_new.sh", "w")
        fw.writelines(AA)
        fw.close()
        
        subprocess.call(["chmod", "u+x", "run_new.sh"])
        subprocess.call(["qsub", "run_new.sh"])
                    
         
