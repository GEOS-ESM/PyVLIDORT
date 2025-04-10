#!/usr/bin/env python3

"""
    Calculates optical property inputs for VLIDORT

    Adapted from sbg_vlidort.py
    Patricia Castellanos, Feb 2025

    Adapted from accp_polar_vlidort.py
    Patricia Castellanos, Jul 2024

    Adapted from polar_vlidort.py and lidar_vlidort.py
    Patricia Castellanos, Jan 2020

"""
import numpy   as np
from pyobs import mietable as mt
from pyobs.aop import G2GAOP
import yaml
from pyobs.constants import MAPL_GRAV as GRAV
from pyobs.constants import MAPL_RDRY as RGAS
from pyobs.constants import MAPL_KAPPA as KAPPA
from py_leo_vlidort.vlidort import WrapperFuncs


class INPUTS_VLIDORT(G2GAOP):
    """
    Routines for calculating optical properties for running vlidort
    aer: GEOS-5 has already been sampled on satellite track and read into aer as a dataset
    wavelength: [nm]
    mtFile = yaml file with aerosol optics tables
    """

    #---
    def getROT(self,wavelength):
        """
        Use f2py ROT_CALC utility to calculate ROT
        wavelength: float or iterable of floats [nm]
        ROT: [nlev,npts,nch]
        depol_ratio: [nch]
        """
        if isinstance(wavelength,float):
            wavelength = [wavelength]
        args = [wavelength, self.pe.astype('float64'), self.ze.astype('float64'), self.te.astype('float64'), self.MISSING, self.verbose]
        vlidortWrapper = WrapperFuncs['ROT_CALC']
        self.ROT, self.depol_ratio, rc = vlidortWrapper(*args)        


    #---
    def getpyobsAOP(self,wavelength):
        """
        use pyobs utilities to get AOP
        wavelength: [nm]
        """
        aop = self.getAOPrt(wavelength=wavelength,vector=True)

        # need to reshape these to [nlev,nch,nobs]
        if isinstance(wavelength,float):
            aop = aop.expand_dims(dim={"ch": 1},axis=1)

        self.tau = aop.AOT.astype('float64').transpose().to_numpy()
        self.ssa = aop.SSA.astype('float64').transpose().to_numpy()
        self.pmom = aop.PMOM.astype('float64').transpose('lev','ch','nobs','m','p').to_numpy()


    #---
    def getEdgeVars(self):
        """
        Calculate atmospheric profile properties needed for Rayleigh calc
        Get altitude, pressure, and temperature and edge of layers
        Assumes dimensions of aer dataset is npts,nlev
        """

        # Get layer thicnkness from DELP & AIRDENS
        # DELP: Pa = kg m-1 s-2
        # AIRDENS: kg m-3
        # GRAV: m s-2
        # -----------------------------------------
        rhodz = self.aer['DELP'] / GRAV
        dz = rhodz / self.aer['AIRDENS']       # column thickness in m

        # add up the thicknesses to get edge level altitudes and pressures
        npts, nlev = dz.shape
        ze = np.array([dz[:,i:].sum(axis=1) for i in range(nlev)])
        # append surface level, altitude = 0
        # [nlev+1,nacross]
        self.ze = np.append(ze,np.zeros([1,npts]),axis=0)

        ptop = 1. # Pa
        pe = np.zeros([nlev+1,npts])
        pe[0,:] = ptop
        for ilev in range(nlev):
            pe[ilev+1,:] = pe[ilev,:] + self.aer['DELP'][:,ilev]

        self.pe = pe

        # get the mid-level pressures and temperatures
        self.pm = (self.pe[:-1,:] + self.pe[1:,:])*0.5
        self.tm = self.pm/(self.aer['AIRDENS'].T*RGAS)

        # get the edge level temperature
        self.te = np.zeros([nlev+1,npts])
        self.te[0,:] = self.tm[0,:]  #isothermal a highest level
        for ilev in range(1,nlev):
            alpha = np.log(self.pe[ilev,:]/self.pm[ilev-1,:])/np.log(self.pm[ilev,:]/self.pm[ilev-1,:])
            self.te[ilev,:] = self.tm[ilev-1] + alpha*(self.tm[ilev,:] - self.tm[ilev-1,:])

        # dry adiabatic
        self.te[nlev,:] = self.tm[nlev-1,:]*(self.pe[nlev,:]/self.pm[nlev-1,:])**KAPPA

    #---
    def getMie(self):
        """
        Use pyobs utilities to load aerosol optics tables
        mtFile = yaml file with aerosol optics tables
        """
        self.mieTable = yaml.safe_load(open(self.mtFile))
        self.vector = True

        # load mietables
        # ---------------------------
        self.p, self.m = 0,0
        for s in self.mieTable:
            m = self.mieTable[s]
            m['mie'] = mt.MIETABLE(m['monoFile'])

        
            dims = dict(self.mieTable[s]['mie'].ds.sizes)
            self.p = max(self.p,dims['p'])
            self.m = max(self.m,dims['m'])


#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":


    pass
