#!/usr/bin/env python
# -*- coding: utf-8 -*-

from classy import Class

class FiducialHubbleParameterAndAngularDiameterDistance():
    def __init__(self, redshift, params_fid):
        self.H_fid = 0.0
        self.Da_fid = 0.0
        if redshift < 1.0e-10:
            self.redshift = 1.0e-10
        else:
            self.redshift = redshift
        self._c_ = 2.99792458e5  # [ km/s ] ##
        self.params_fid = params_fid

    def calcHubbleParameterAndAngularDiameterDistance(self):

        cosmo_temp = Class()

        ##########################################
        ########## fiducial parameter ############
        ##########################################

        h = self.params_fid["h"]
        ### fiducial Dz and H ###
        cosmo_temp.set(self.params_fid)
        cosmo_temp.compute()
        self.Da_fid = cosmo_temp.angular_distance(self.redshift) * h
        self.H_fid = cosmo_temp.Hubble(self.redshift) * self._c_
        cosmo_temp.struct_cleanup()
        cosmo_temp.empty()

    def getFiducialAngularDiameterDistance(self):
        return self.Da_fid

    def getFiducialHubbleParameter(self):
        return self.H_fid

