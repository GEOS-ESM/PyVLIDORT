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
from py_vlidort import AOP_
from pyobs.mietable import _rh, _drh  # Import the RH array and step size

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
    def getAOP(self,wavelength,vector=True,Species=None,do_g=True):

        # Tables must have be consistent across species
        # ---------------------------------------------
        if vector:
            if not self.vector:
                print('Warning: will not calculate PMOM because of inconsistent Mie Tables.')
                vector = False

        # All species on file or a subset
        # -------------------------------
        if Species is None:
            Species = list(self.mieTable.keys())
        if isinstance(Species,str):
            Species = [Species,]

        # Identify only the variables we strictly need
        req_vars = ['RH', 'AIRDENS']
        req_vars.append('DELP' if 'DELP' in self.aer else 'delp')

        for s in Species:
            req_vars.extend(self.mieTable[s]['tracers'])

        # Subset the dataset and load ONLY those variables into RAM
        a = self.aer[req_vars].load()

        # pre-load RH, AIRDENS, and DELP and Species so you don't hit dask
        # repeatedly looping through  AOP calculations
        # Extract underlying numpy arrays to skip Xarray overhead in the loop
        # -------------------------------------------------------
        rh = a['RH'].values.copy()
        airdens = a['AIRDENS'].values
        dp = a['DELP'].values if 'DELP' in a else a['delp'].values

        # Save original dims and coords to build the final Xarray Dataset at the end
        rh_dims = a['RH'].dims
        rh_coords = a['RH'].coords

        # Handy arrays for extensive properties
        # -------------------------------------
        rhodz = dp / GRAV
        dz = rhodz / a['AIRDENS']       # column thickness

        # Relevant dimensions
        # -------------------
        space = rh.shape
        aot, sca, g = np.zeros(space, dtype=np.float32), np.zeros(space, dtype=np.float32), np.zeros(space, dtype=np.float32)
        if vector:
            ns = np.prod(space)
            p = self.p
            ang = self.ang
            pmatrix = np.zeros((ns,p,ang), dtype=np.float32) #flatten space dimensions for convenience

        ns = np.prod(space) # flattened size
        rh_flat = rh.ravel().astype(np.float64)
        
        # Flatten accumulators for the Fortran call
        aot_flat = np.zeros(ns, dtype=np.float64,order='F')
        sca_flat = np.zeros(ns, dtype=np.float64,order='F')
        g_flat = np.zeros(ns, dtype=np.float64,order='F')

        if vector:
            do_vector = 1
            pmatrix_flat = np.zeros((self.p, self.ang, ns), dtype=np.float64,order='F')
        else:
            do_vector = 0
            # Dummy arrays to satisfy Fortran arguments without wasting memory
            pmatrix_flat = np.zeros((1, 1, 1), dtype=np.float64,order='F')


        for s in Species:
            Tracers = self.mieTable[s]['tracers']
            mie = self.mieTable[s]['mie']
            bins = self.mieTable[s].get('bin', list(range(1, len(Tracers)+1)))
            
            for q, bin in zip(Tracers, bins):
                # Flatten tracer mass
                q_mass_flat = np.asfortranarray((rhodz * a[q].values).ravel().astype(np.float64))
                
                # Extract 1D arrays from optics file for this bin & wavelength
                # (You extract these once per bin from mie.ds)
                bin_idx = bin - 1
                lut_bext = np.asfortranarray(mie.ds['bext'].sel(wavelength=wavelength/1e9, method='nearest').isel(bin=bin_idx).interp(rh=_rh).values.astype(np.float64))
                lut_bsca = np.asfortranarray(mie.ds['bsca'].sel(wavelength=wavelength/1e9, method='nearest').isel(bin=bin_idx).interp(rh=_rh).values.astype(np.float64))
                lut_g = np.asfortranarray(mie.ds['g'].sel(wavelength=wavelength/1e9, method='nearest').isel(bin=bin_idx).interp(rh=_rh).values.astype(np.float64))
                
                if vector:
                    # Concat pmatrix LUTs
                    p11 = mie.ds['p11'].sel(wavelength=wavelength/1e9, method='nearest').isel(bin=bin_idx)
                    p22 = mie.ds['p22'].sel(wavelength=wavelength/1e9, method='nearest').isel(bin=bin_idx)
                    p33 = mie.ds['p33'].sel(wavelength=wavelength/1e9, method='nearest').isel(bin=bin_idx)
                    p44 = mie.ds['p44'].sel(wavelength=wavelength/1e9, method='nearest').isel(bin=bin_idx)
                    p12 = mie.ds['p12'].sel(wavelength=wavelength/1e9, method='nearest').isel(bin=bin_idx)
                    p34 = mie.ds['p34'].sel(wavelength=wavelength/1e9, method='nearest').isel(bin=bin_idx)
                
                    lut_pmatrix = xr.concat((p11, p22, p33, p44, p12, p34), 'p').interp(rh=_rh).values.astype(np.float64)
                    # Ensure shape is (p, ang, nrh) for Fortran
                    lut_pmatrix = np.asfortranarray(np.moveaxis(lut_pmatrix, 1, 2) )
                else:
                    # Dummy LUT to satisfy Fortran arguments
                    lut_pmatrix = np.zeros((1, 1, len(_rh)), dtype=np.float64,order='F')


                # Call Fortran 
                AOP_.aop_accum(
                    rh_flat, _drh, q_mass_flat, 
                    lut_bext, lut_bsca, lut_g, lut_pmatrix, 
                    do_vector, aot_flat, sca_flat, g_flat, pmatrix_flat
                )

        # Reshape back to 3D arrays after loop
        aot = aot_flat.reshape(space)
        sca = sca_flat.reshape(space)
        g = g_flat.reshape(space)
        if vector:
            pmatrix = np.moveaxis(pmatrix_flat, -1, 0).reshape(space + (self.p, self.ang))


        # Final normalization of SSA and g
        # protect against divide by zero
        # this can happen if you ask for the AOP of an individual species
        # and its' concentration in a layer is zero
        # --------------------------------
        ssa_flat = np.full(ns, np.nan, dtype=np.float64)
        
        # Normalize SSA
        I = np.where(aot_flat != 0.0)[0]
        ssa_flat[I] = sca_flat[I] / aot_flat[I]

        # Normalize Pmatrix
        if vector:
            I = np.where(sca_flat != 0.0)[0]
            # pmatrix_flat is (p, ang, ns). Divide along the last axis
            pmatrix_flat[:, :, I] = pmatrix_flat[:, :, I] / sca_flat[I]
            
            I_zero = np.where(sca_flat == 0.0)[0]
            pmatrix_flat[:, :, I_zero] = np.nan

        # Normalize Asymmetry Parameter
        if do_g or not vector:
             I = np.where(sca_flat != 0.0)[0]
             g_flat[I] = g_flat[I] / sca_flat[I]
             
             I_zero = np.where(sca_flat == 0.0)[0]
             g_flat[I_zero] = np.nan

        # Reshape everything back to the original grid dimensions
        aot = aot_flat.reshape(space)
        sca = sca_flat.reshape(space)
        ssa = ssa_flat.reshape(space)
        g = g_flat.reshape(space)
        
        if vector:
            # Move the spatial 'ns' axis back to the front, then reshape
            pmatrix = np.moveaxis(pmatrix_flat, -1, 0).reshape(space + (self.p, self.ang))

        # Pack results into a Dataset
        # ---------------------------
        A = dict (AOT = {'long_name':'Aerosol Optical Thickness', 'units':'1'},
                  SSA = {'long_name':'Aerosol Single Scattering Albedo', 'units':'1'},
                  G = {'long_name':'Aerosol Asymmetry Parameter', 'units':'1'},
                  PMOM = {'long_name':'Aerosol Phase Matrix Expansion Coefficients', 'units':'1'},
                  PMATRIX = {'long_name':'Aerosol Phase Matrix (non-zero elements)', 'units':'1'}
                  )

        DA = dict( AOT = xr.DataArray(aot,dims=rh_dims,coords=rh_coords,attrs=A['AOT']),
                   SSA = xr.DataArray(ssa,dims=rh_dims,coords=rh_coords,attrs=A['SSA']),
                 )

        DA['DELP'] = xr.DataArray(dp,dims=rh_dims,coords=rh_coords)
        DA['AIRDENS'] = xr.DataArray(airdens,dims=rh_dims,coords=rh_coords)

        if vector:
            coords = dict(rh_coords).copy()
            coords['p'] = mie.ds.coords['p']
            coords['ang'] = mie.ds.coords['ang']
            DA['PMATRIX'] = xr.DataArray(pmatrix, dims=rh_dims+('p','ang'),coords=coords,attrs=A['PMATRIX'])

        if do_g or not vector:
            DA['G'] = xr.DataArray(g,dims=rh_dims,coords=rh_coords,attrs=A['G'])

        aop =  xr.Dataset(DA)

        self.tau = aop.AOT.expand_dims('ch').astype('float64').transpose('lev','ch','nobs').to_numpy()
        self.ssa = aop.SSA.expand_dims('ch').astype('float64').transpose('lev','ch','nobs').to_numpy()
        self.g   = aop.G.expand_dims('ch').astype('float64').transpose('lev','ch','nobs').to_numpy()
        self.pmatrix = aop.PMATRIX.expand_dims('ch').astype('float64').transpose('lev','ch','nobs','ang','p').to_numpy()
        self.angles = aop.ang.to_numpy()

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
        self.angles = aop.ang.to_numpy()


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
