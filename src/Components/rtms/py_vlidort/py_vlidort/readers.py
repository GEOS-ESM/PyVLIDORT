#!/usr/bin/env python3

"""
    Collection of reader helper subroutines

    Patricia Castellanos, Jul 2026

"""
import numpy   as np
import xarray as xr

class READERS(object):
    #--
    def getDims(self):
        """
        Get granule dimensions
        """
        filename = self.inFile.replace('%col', 'aer_Nv')
        if self.verbose:
            print(f'opening file {filename}')

        with xr.open_dataset(filename) as ds:
            self.ntyme = ds.sizes.get('time', 1)
            self.nlev = ds.sizes.get('lev', 1)

            # get across/along dimensions if they exist
            self.ncross = ds.sizes.get('ncross', 1)
            self.nalong = ds.sizes.get('nalong', 1)

        self.nobs = self.ntyme * self.ncross * self.nalong
        self.nbatch = self.nproc

    #--
    def readAngles(self):
        """
        Read in viewing and solar Geometry from angFile
        """
        filename = self.inFile.replace('%col', self.instname)
        if self.verbose:
            print(f'opening file {filename}')

        ds = xr.open_dataset(filename, chunks="auto").squeeze()

        # figure out what to stack into 'nobs'
        stack_dims = [d for d in ['time', 'ncross', 'nalong'] if d in ds.dims]
        if len(stack_dims) > 1:
            ds = ds.stack(nobs=stack_dims)
        elif len(stack_dims) == 1:
            ds = ds.rename({stack_dims[0]: 'nobs'})

        # Map variables dynamically (using dict.get for alias fallback)
        for sds in self.SDS_ANG:
            setattr(self, sds, ds[self.ncALIAS.get(sds, sds)])

        # define RAA according to photon travel direction
        saa = self.SAA + 180.0
        I = saa >= 360.
        saa[I.compute()] = saa[I.compute()] - 360.

        RAA = self.VAA - saa
        I = RAA < 0
        RAA[I.compute()] = RAA[I.compute()]+360.0
        self.RAA = RAA

        # Limit SZAs
        iGood = self.SZA < 80
        self.iGood = self.iGood & iGood.values
        self.nobs = np.sum(self.iGood)

    #---
    def readSampledGEOS(self):
        """
        Read in model sampled track
        """

        cols = ['aer_Nv'] + (['met_Nv'] if len(self.SDS_MET) > 0 else [])
        inList = [self.inFile.replace('%col', col) for col in cols]

        if self.verbose:
            for f in inList:
                print(f'opening file {f}')

        self.AER = xr.open_mfdataset(inList, chunks="auto").squeeze()

        # make arrays [nobs,nlev]
        # dynamically stack dimensions into 'nobs'
        stack_dims = [d for d in ['time', 'ncross', 'nalong'] if d in self.AER.dims]
        if len(stack_dims) > 1:
            self.AER = self.AER.stack(nobs=stack_dims)
        elif len(stack_dims) == 1:
            self.AER = self.AER.rename({stack_dims[0]: 'nobs'})

        # the ellipsis (...) safely ignores variables that might not have a 'lev' dimension
        self.AER = self.AER.transpose("nobs", "lev", ...)
        iGood = np.arange(len(self.iGood))[self.iGood]
        self.AER = self.AER.isel(nobs=iGood)

    # ---
    def LandSeaMask(self):
        """
        Read in invariant dataset
        """
        filename = self.inFile.replace('%col', 'asm_Nx')
        if self.verbose:
            print(f'opening file {filename}')

        ds = xr.open_dataset(filename, chunks="auto").squeeze()

        # stack dimensions into 'nobs'
        stack_dims = [d for d in ['time', 'ncross', 'nalong'] if d in ds.dims]
        if len(stack_dims) > 1:
            ds = ds.stack(nobs=stack_dims)
        elif len(stack_dims) == 1:
            ds = ds.rename({stack_dims[0]: 'nobs'})

        # map variables dynamically (using dict.get for alias fallback)
        for sds in self.SDS_INV:
            setattr(self, sds, ds[self.ncALIAS.get(sds, sds)])

        iGood = self.FRLAND >= 0.99
        self.iLand = iGood.values
        iGood = self.FRLAND < 0.99
        self.iSea  = iGood.values

        # self.iGood = self.iGood & iGood
        self.nobsLand  = np.sum(self.iGood & self.iLand)
        self.nobsSea   = np.sum(self.iGood & self.iSea)
