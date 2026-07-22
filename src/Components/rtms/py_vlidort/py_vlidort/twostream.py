#!/usr/bin/env python3

"""
    Parent class with utilities needed to 
    setup input data
    and run python wrappers for twostream
"""

from   py_vlidort import TWOSTREAM_ 
import numpy as np 

#---
def LAMBERTIAN_run(args):

    # Call TWO_STREAM wrapper function
    I, reflectance, rc = TWOSTREAM_.twostream_lambert_driver(*args)

    return I,reflectance

#---
def MODIS_BRDF_run(args):

    # Call TWO_STREAM wrapper function
    I, reflectance, rc = TWOSTREAM_.twostream_brdf_rtls(*args)

    return I,reflectance

"""
Utilities to help with packing and upacking variables
"""
def ts_pack_args_LAMB(self,channel,args,npts):
    # Atmospheric Optical Property Arguments
    rot,depol_ratio,alpha,tau,ssa,g,pmatrix,tauI,ssaI,gI,pmatrixI,tauL,ssaL,gL,pmatrixL,pe,te,ze,vza,sza,raa,flux_factor,albedo = args


    # create list of input arguments
    packed_args = [(channel, self.plane_parallel, self.angles,
            rot[:,i:i+1,:], depol_ratio, alpha[:,:,i:i+1],
            tau[:,:,i:i+1], ssa[:,:,i:i+1], g[:,:,i:i+1], pmatrix[:,:,i:i+1,:,:],
            tauI[:,:,i:i+1], ssaI[:,:,i:i+1], gI[:,:,i:i+1], pmatrixI[:,:,i:i+1,:,:],
            tauL[:,:,i:i+1], ssaL[:,:,i:i+1], gL[:,:,i:i+1], pmatrixL[:,:,i:i+1,:,:],
            pe[:,i:i+1], ze[:,i:i+1], te[:,i:i+1],
            albedo,
            [sza[i:i+1]], [raa[i:i+1]], [vza[i:i+1]],
            flux_factor,
            self.MISSING,
            self.verbose, self.debug) for i in range(npts)]



    return packed_args

def ts_pack_args_MODIS_BRDF(self,channel,args,npts):
    # Atmospheric Optical Property Arguments
    rot,depol_ratio,alpha,tau,ssa,g,pmatrix,tauI,ssaI,gI,pmatrixI,tauL,ssaL,gL,pmatrixL,pe,te,ze,param,kernel_wt,vza,sza,raa,flux_factor = args

    # create list of input arguments
    packed_args = [(channel, self.plane_parallel, self.angles,
            rot[:,i:i+1,:], depol_ratio, alpha[:,:,i:i+1],
            tau[:,:,i:i+1], ssa[:,:,i:i+1], g[:,:,i:i+1], pmatrix[:,:,i:i+1,:,:],
            tauI[:,:,i:i+1], ssaI[:,:,i:i+1], gI[:,:,i:i+1], pmatrixI[:,:,i:i+1,:,:],
            tauL[:,:,i:i+1], ssaL[:,:,i:i+1], gL[:,:,i:i+1], pmatrixL[:,:,i:i+1,:,:],
            pe[:,i:i+1], ze[:,i:i+1], te[:,i:i+1],
            kernel_wt[:,:,i:i+1], param[:,:,i:i+1],
            [sza[i:i+1]], [raa[i:i+1]], [vza[i:i+1]],
            flux_factor,
            self.MISSING,
            self.verbose, self.debug) for i in range(npts)]



    return packed_args

def ts_unpack_result(result):
    # initialize temporary outputs
    I = []
    reflectance = []

    # reshape outputs
    for r in result:
        I_r,reflectance_r = r
        I.append(I_r)
        reflectance.append(reflectance_r)
    I = np.concatenate(I)
    reflectance = np.concatenate(reflectance)

    return I, reflectance

#---
def ts_initOutputs(self):
    # Initiate output arrays
    nch    = self.nch
    nobs   = self.nobs

    self.ts_I = np.ones([nobs,nch])*self.MISSING
    self.ts_reflectance = np.ones([nobs,nch])*self.MISSING
