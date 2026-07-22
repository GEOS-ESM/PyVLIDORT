#!/usr/bin/env python3

"""
    Parent class with utilities needed to 
    setup input data
    and run python wrappers for vlidort
"""

import os
from   py_vlidort import VLIDORT_POLAR_ 
import numpy   as np


WrapperFuncs = {'MODIS_BRDF_PMATRIX'  : VLIDORT_POLAR_.brdf_modis_pmatrix,
                'ROT_CALC'            : VLIDORT_POLAR_.rot_calc,
                'MODIS_BRDF_PMOM'     : VLIDORT_POLAR_.brdf_modis_pmom,
#                'MODIS_BRDF_BPDF': VLIDORT_POLAR_.vector_brdf_modis_bpdf,
                'LAMBERTIAN_PMOM'     : VLIDORT_POLAR_.lambert_pmom,
                'LAMBERTIAN_PMATRIX'  : VLIDORT_POLAR_.lambert_pmatrix}
#                'GissCX'         : VLIDORT_POLAR_.vector_gisscx,
#                'CX'             : VLIDORT_POLAR_.vector_cx,
#                'CX_CLOUD'             : VLIDORT_POLAR_.vector_cx_cloud,

#---
def CX_run(args):
       
    # Call VLIDORT wrapper function
    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = VLIDORT_POLAR_.vector_cx(*args)

    return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U
#---
def CX_CLOUD_run(args):

    # Call VLIDORT wrapper function
    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = VLIDORT_POLAR_.vector_cx_cloud(*args)

    return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U

#---
def LAMBERTIAN_PMOM_run(args):

    # Call VLIDORT wrapper function
    I, reflectance, Q, U, rc = VLIDORT_POLAR_.lambert_pmom(*args)

    surf_reflectance = None
    BR_Q = None
    BR_U = None
    return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U
#---
def LAMBERTIAN_PMATRIX_run(args):

    # Call VLIDORT wrapper function
    I, reflectance, Q, U, rc = VLIDORT_POLAR_.lambert_pmatrix(*args)

    surf_reflectance = None
    BR_Q = None
    BR_U = None
    return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U

#---
def MODIS_BRDF_PMOM_run(args):

    # Call VLIDORT wrapper function
    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = VLIDORT_POLAR_.brdf_modis_pmom(*args)

    return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U
#---
def MODIS_BRDF_PMATRIX_run(args):

    # Call VLIDORT wrapper function
    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = VLIDORT_POLAR_.brdf_modis_pmatrix(*args)

    return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U

""" 
Utilities to help with packing and upacking variables
"""
def vl_pack_args_LAMB(self,channel,args,NSTOKES,npts):
    # Atmospheric Optical Property Arguments
    rot,depol_ratio,alpha,tau,ssa,g,pmatrix,tauI,ssaI,gI,pmatrixI,tauL,ssaL,gL,pmatrixL,pe,te,ze,vza,sza,raa,flux_factor,albedo = args

    # create list of input arguments
    packed_args = [(channel, self.nstreams, self.plane_parallel,NSTOKES,self.angles,
            rot[:,i:i+1,:], depol_ratio, alpha[:,:,i:i+1],
            tau[:,:,i:i+1], ssa[:,:,i:i+1], pmatrix[:,:,i:i+1,:,:],
            tauI[:,:,i:i+1], ssaI[:,:,i:i+1], pmatrixI[:,:,i:i+1,:,:],
            tauL[:,:,i:i+1], ssaL[:,:,i:i+1], pmatrixL[:,:,i:i+1,:,:],
            pe[:,i:i+1], ze[:,i:i+1], te[:,i:i+1],
            albedo,
            [sza[i:i+1]], [raa[i:i+1]], [vza[i:i+1]],
            flux_factor,
            self.MISSING,
            self.verbose, self.debug) for i in range(npts)]


    return packed_args

def vl_pack_args_MODIS_BRDF(self,channel,args,NSTOKES,npts):
    # Atmospheric Optical Property Arguments
    rot,depol_ratio,alpha,tau,ssa,g,pmatrix,tauI,ssaI,gI,pmatrixI,tauL,ssaL,gL,pmatrixL,pe,te,ze,param,kernel_wt,vza,sza,raa,flux_factor = args


    # create list of input arguments
    packed_args = [(channel, self.nstreams, self.plane_parallel,NSTOKES,self.angles,
            rot[:,i:i+1,:], depol_ratio, alpha[:,:,i:i+1],
            tau[:,:,i:i+1], ssa[:,:,i:i+1], pmatrix[:,:,i:i+1,:,:],
            tauI[:,:,i:i+1], ssaI[:,:,i:i+1], pmatrixI[:,:,i:i+1,:,:],
            tauL[:,:,i:i+1], ssaL[:,:,i:i+1], pmatrixL[:,:,i:i+1,:,:],
            pe[:,i:i+1], ze[:,i:i+1], te[:,i:i+1],
            kernel_wt[:,:,i:i+1], param[:,:,i:i+1],
            [sza[i:i+1]], [raa[i:i+1]], [vza[i:i+1]],
            flux_factor,
            self.MISSING,
            self.verbose, self.debug) for i in range(npts)]


    return packed_args


def vl_unpack_result(result):
    # Unpack pooled results
    # initialize temporary outputs
    I = []
    Q = []
    U = []
    reflectance = []
    surf_reflectance = []
    BR_Q = []
    BR_U = []

    # reshape outputs
    for r in result:
        I_r,Q_r,U_r,reflectance_r,surf_reflectance_r,BR_Q_r,BR_U_r = r
        I.append(I_r)
        Q.append(Q_r)
        U.append(U_r)
        reflectance.append(reflectance_r)
        surf_reflectance.append([surf_reflectance_r])
        BR_Q.append([BR_Q_r])
        BR_U.append([BR_U_r])
    I = np.concatenate(I)
    Q = np.concatenate(Q)
    U = np.concatenate(U)
    reflectance = np.concatenate(reflectance)
    surf_reflectance = np.concatenate(surf_reflectance)
    BR_Q = np.concatenate(BR_Q)
    BR_U = np.concatenate(BR_U)

    return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U

#---
def vl_initOutputs(self):
    # Initiate output arrays
    nch    = self.nch
    nobs   = self.nobs

    self.I = np.ones([nobs,nch])*self.MISSING
    self.Q = np.ones([nobs,nch])*self.MISSING
    self.U = np.ones([nobs,nch])*self.MISSING
    self.reflectance = np.ones([nobs,nch])*self.MISSING
    self.surf_reflectance = np.ones([nobs,nch])*self.MISSING
    self.BR_Q = np.ones([nobs,nch])*self.MISSING
    self.BR_U = np.ones([nobs,nch])*self.MISSING


def get_chd(channel):
    chd = '%.2f'%channel
    chd = chd.replace('.','d')

    return chd
