#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess

output_dir = "chains_base"
#output_dir = "chains_main"

GALAXY = ["galaxy"]

MOCK = []
for i in range(100):
    MOCK.append("mock_%04d" % (i+1))

for param in ["GENERAL_North_zbin1", "GENERAL_South_zbin1", "GENERAL_North_zbin3", "GENERAL_South_zbin3"]:

#    for data_type in ["galaxy"] + MOCK:
    for data_type in ["galaxy"]:
 
        print("%s_%s" % (data_type, param))

        subprocess.call(["python", "MontePython.py", "info", "%s/%s/%s_%s" % (output_dir, param, data_type, param), "--noplot", "--noplot-2d", "--want-covmat"])
    
