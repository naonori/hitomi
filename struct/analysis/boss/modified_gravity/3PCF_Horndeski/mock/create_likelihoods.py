#######################
# Purpose of this code
#######################
# The code makes the likelihood directory for mock catalogues.
#######################

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os
import numpy as np
import shutil

MOCK = []
for i in range(100):
    MOCK.append("mock_%04d" % (i+1))

for data_type in MOCK:

    model_type = "mock"

    for lkl in ["Horndeski_North_zbin1_case7", "Horndeski_North_zbin1_case7_NS", "Horndeski_North_zbin1_case7_ZBIN", "Horndeski_North_zbin1_case7_TOTAL",\
                "Horndeski_South_zbin1_case7", "Horndeski_South_zbin1_case7_NS", "Horndeski_South_zbin1_case7_ZBIN", "Horndeski_South_zbin1_case7_TOTAL",\
                "Horndeski_North_zbin3_case7", "Horndeski_North_zbin3_case7_NS", "Horndeski_North_zbin3_case7_ZBIN", "Horndeski_North_zbin3_case7_TOTAL",\
                "Horndeski_South_zbin3_case7", "Horndeski_South_zbin3_case7_NS", "Horndeski_South_zbin3_case7_ZBIN", "Horndeski_South_zbin3_case7_TOTAL"]:

        ##############################
        # Make a likelihoods directory
        
        likelihoods = "%s_%s"  % (data_type, lkl)
        
        print(likelihoods)
        
        try:
            shutil.rmtree("likelihoods/%s" % likelihoods)
        except:
            pass
        
        shutil.copytree("likelihoods/galaxy_%s" % lkl, "likelihoods/%s" % likelihoods)
        shutil.move("likelihoods/%s/galaxy_%s.data" % (likelihoods, lkl), "likelihoods/%s/%s.data" % (likelihoods, likelihoods))
        
        ##############################
        # Edit the .data file in the likelihood directory.
        
        fr = open("likelihoods/%s/%s.data" % (likelihoods, likelihoods), "r")
        AA = fr.readlines()
        fr.close()
 
        for i in range(len(AA)):
            AA[i] = AA[i].replace("galaxy_%s" % lkl, likelihoods)
        
        for i in range(len(AA)):
            if AA[i].find("%s.data_type" % likelihoods) >= 0:
                AA[i] = '%s.data_type  = "%s"\n' % (likelihoods, data_type)
                
            if AA[i].find("%s.model_type" % likelihoods) >= 0:
                AA[i] = '%s.model_type  = "%s"\n' % (likelihoods, model_type)
        
        fw = open("likelihoods/%s/%s.data" % (likelihoods, likelihoods), "w")
        fw.writelines(AA)
        fw.close()
         
        ##############################
        # Edit the "__init__.py" in the likelihood directory.
        
        fr = open("likelihoods/%s/__init__.py" % (likelihoods), "r")
        AA = fr.readlines()
        fr.close()

        for i in range(len(AA)):
            AA[i] = AA[i].replace("galaxy_%s" % lkl, likelihoods)

        fw = open("likelihoods/%s/__init__.py" % (likelihoods), "w")
        fw.writelines(AA)
        fw.close()
 

##############################
# Copy and edit "input/xxx.param".

for data_type in MOCK:

    for param in ["Horndeski.param", "Horndeski_North.param", "Horndeski_South.param", "Horndeski_zbin1.param", "Horndeski_zbin3.param",\
                  "Horndeski_North_zbin1.param", "Horndeski_South_zbin1.param", "Horndeski_North_zbin3.param", "Horndeski_South_zbin3.param"]:

        shutil.copy("input/galaxy_%s" % param, "input/%s_%s" % (data_type, param))
        
        fr = open("input/%s_%s" % (data_type, param), "r")
        AA = fr.readlines()
        fr.close()
         
        for i in range(len(AA)):
            AA[i] = AA[i].replace("galaxy", data_type)
        
        fw = open("input/%s_%s" % (data_type, param), "w")
        fw.writelines(AA)
        fw.close()
                
                
