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
import xarray as xr
from pyobs import mietable as mt
from pyobs.aop import G2GAOP
import yaml
from pyobs.constants import MAPL_GRAV as GRAV
from pyobs.constants import MAPL_RDRY as RGAS
from pyobs.constants import MAPL_KAPPA as KAPPA
from .vlidort import WrapperFuncs


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
        args = [wavelength, self.EDGES.pe.astype('float64'), self.EDGES.ze.astype('float64'), self.EDGES.te.astype('float64'), self.MISSING, self.verbose]
        vlidortWrapper = WrapperFuncs['ROT_CALC']
        ROT, depol_ratio, rc = vlidortWrapper(*args)        

        # Format wavelength as a 1D array for the coordinate
        # (e.g., 550.0 becomes [550.0])
        wav_coord = np.atleast_1d(wavelength)


        # Pack ROT and depol_ratio into an xarray Dataset called RAYLEIGH
        self.RAYLEIGH = xr.Dataset(
            {
                'ROT': (('lev', 'nobs', 'nch'), ROT),
                'depol_ratio': (('nch',), depol_ratio)
            },
            coords={
                # Inherit exactly the same coordinates for lev and nobs
                'lev': self.AER.coords['lev'],
                'nobs': self.AER.coords['nobs'],
                'nch': wav_coord 
            },
            attrs={
                'description': 'Rayleigh optical depth and depolarization ratio'
            }
        )


    #---
    def getpyobsAOP(self,wavelength):
        """
        use pyobs utilities to get AOP
        wavelength: [nm]
        """
        aop = self.getAOPrt(wavelength=wavelength,vector=True,do_g=True)

        # need to reshape these to [nlev,nch,nobs]
        aop = aop.expand_dims(dim={"ch": 1},axis=1)

        self.tau = aop.AOT.astype('float64').transpose('lev','ch','nobs').to_numpy()
        self.ssa = aop.SSA.astype('float64').transpose('lev','ch','nobs').to_numpy()
        self.g   = aop.G.astype('float64').transpose('lev','ch','nobs').to_numpy()
        self.pmatrix = aop.PMATRIX.astype('float64').transpose('lev','ch','nobs','ang','p').to_numpy()
        self.angles = aop.ang


    #---
    def getEdgeVars(self):
        """
        Calculate atmospheric profile properties needed for Rayleigh calc
        Get altitude, pressure, and temperature and edge of layers
        Assumes dimensions of aer dataset is npts,nlev
        """

        # Get layer thicnkness from DELP & AIRDENS
        # DELP: Pa = kg m-1 s-2 [npts,nlev]
        # AIRDENS: kg m-3       [npts,nlev]
        # GRAV: m s-2

        # output
        # te [nlev+,npts]
        # pe [nlev+,npts]
        # ze [nlev+1,npts]
        # -----------------------------------------
        rhodz = self.AER['DELP'] / GRAV
        dz = rhodz / self.AER['AIRDENS']       # column thickness in m

        # add up the thicknesses to get edge level altitudes and pressures
        npts, nlev = dz.shape
        ze = np.array([dz[:,i:].sum(axis=1) for i in range(nlev)])
        # append surface level, altitude = 0
        # [nlev+1,nacross]
        ze = np.append(ze,np.zeros([1,npts]),axis=0)

        ptop = 1. # Pa
        pe = np.zeros([nlev+1,npts])
        pe[0,:] = ptop
        for ilev in range(nlev):
            pe[ilev+1,:] = pe[ilev,:] + self.AER['DELP'][:,ilev]

        # get the mid-level pressures and temperatures
        pm = (pe[:-1,:] + pe[1:,:])*0.5
        tm = pm/(self.AER['AIRDENS'].T*RGAS)

        # get the edge level temperature
        te = np.zeros([nlev+1,npts])
        te[0,:] = tm[0,:]  #isothermal a highest level
        for ilev in range(1,nlev):
            alpha = np.log(pe[ilev,:]/pm[ilev-1,:])/np.log(pm[ilev,:]/pm[ilev-1,:])
            te[ilev,:] = tm[ilev-1] + alpha*(tm[ilev,:] - tm[ilev-1,:])

        # dry adiabatic
        te[nlev,:] = tm[nlev-1,:]*(pe[nlev,:]/pm[nlev-1,:])**KAPPA

        # Create a new coordinate for the edge levels (1 to nlev+1)
        leve_coord = np.arange(1, nlev + 2, dtype=np.float32)
        
        self.EDGES = xr.Dataset(
            {
                'pe': (('leve', 'nobs'), pe),
                'ze': (('leve', 'nobs'), ze),
                'te': (('leve', 'nobs'), te)
            },
            coords={
                'leve': leve_coord,
                # Inherit the nobs MultiIndex (time, ncross) directly from self.AER
                'nobs': self.AER.coords['nobs']
            },
            attrs={
                'description': 'Edge level profiles for pressure (pe), altitude (ze), and temperature (te)'
            }
        )


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
        self.p, self.m, self.ang = 0,0,0
        for s in self.mieTable:
            m = self.mieTable[s]
            m['mie'] = mt.MIETABLE(m['monoFile'])

        
            dims = dict(self.mieTable[s]['mie'].ds.sizes)
            self.p = max(self.p,dims['p'])
            self.ang = max(self.ang,dims['ang'])


#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":


    pass
