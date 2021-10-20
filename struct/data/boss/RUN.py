#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess

## For galaxies
#for NS in ["North", "South"]:
#    for zbin in [1,2,3]:
#        for Weight in [0,1,2,3]:
##        for Weight in [0]:
#
#            fr = open("run_base.sh", "r")
#            AA = fr.readlines()
#            fr.close()
#            
#            AA = [AA[i].replace("NS zbin Weight", "%s %d %d" % (NS, zbin, Weight)) for i in range(len(AA))]
#            
#            fw = open("run_new.sh", "w")
#            fw.writelines(AA)
#            fw.close()
#
#            subprocess.call(["chmod", "u+x", "run_new.sh"])
#            subprocess.call(["qsub", "run_new.sh"])
#
## For mocks
#for NS in ["North", "South"]:
#    for zbin in [1,2,3]:
#       for NR in range(0, 21):
##        for NR in [0]:
#
#            fr = open("run_base.sh", "r")
#            AA = fr.readlines()
#            fr.close()
#            
#            AA = [AA[i].replace("NS zbin NR", "%s %d %d" % (NS, zbin, NR)) for i in range(len(AA))]
#            
#            fw = open("run_new.sh", "w")
#            fw.writelines(AA)
#            fw.close()
#
#            subprocess.call(["chmod", "u+x", "run_new.sh"])
#            subprocess.call(["qsub", "run_new.sh"])
#
