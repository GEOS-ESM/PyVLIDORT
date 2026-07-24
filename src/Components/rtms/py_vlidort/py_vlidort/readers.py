#!/usr/bin/env python3

"""
    Collection of reader helper subroutines

    Patricia Castellanos, Jul 2026

"""
import numpy   as np
import xarray as xr


class READERS(object):

    ALBEDO_READER_MAP = {
        'AMES_BRDF':  'readSampledAMESBRDF',
        #'MODIS_BRDF': 'readSampledMODISBRDF',
        #'LAMBERTIAN': 'readSampledLambertian',
        # Add more as they are created
    }

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
        ang_vars = {}
        for sds in self.SDS_ANG:
            ang_vars[sds] = ds[self.ncALIAS.get(sds, sds)]        

        # define RAA according to photon travel direction
        saa = ang_vars['SAA'] + 180.0
        I = saa >= 360.
        saa[I.compute()] = saa[I.compute()] - 360.

        RAA = ang_vars['VAA'] - saa
        I = RAA < 0
        RAA[I.compute()] = RAA[I.compute()]+360.0

        # Limit SZAs
        iGood = ang_vars['SZA'] < 80
        self.iGood = self.iGood & iGood.values
        self.nobs = np.sum(self.iGood)

        # Pack angles into an xarray Dataset
        self.GEOMS = xr.Dataset({
            'RAA': RAA,
            'SZA': ang_vars['SZA'],
            'SAA': ang_vars['SAA'],
            'VZA': ang_vars['VZA']
        })


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

    #---
    def readSampledAMESBRDF(self):
        """
        Read in AMES BRDF kernel weights
        that have already been sampled on swath
        """
        if self.verbose:
            print('opening BRDF file ',self.brdfFile)
        ds = xr.open_dataset(self.brdfFile,chunks="auto").squeeze()

        # concatenate kernel weights into one data array
        da = xr.concat([ds.Ki,ds.Kg,ds.Kv],dim="nkernel")
        da = da.rename("kernel_wt")

        # identify the spatial dimensions (e.g., ['time', 'ncross'] or ['nalong', 'ncross'])
        spatial_dims = [dim for dim in da.dims if dim not in ["nkernel", "nwav"]]

        # reshape to [nkernel,nch,nobs]
        kernel_wt = da.stack(nobs=spatial_dims).transpose("nkernel", "nwav", "nobs")

        # Create RTLSparam scaled to the full number of observations
        full_nobs = kernel_wt.sizes["nobs"]
        nwav_size = kernel_wt.sizes["nwav"]
 
        # Make array of shape [nparam, nwav, nobs] initialized to [2, 1, ..., full_nobs]
        param = np.zeros((2, nwav_size, full_nobs), dtype=np.float64)
        param[0, :, :] = 2
        param[1, :, :] = 1

        # Convert it to an xarray DataArray so it matches kernel_wt's coordinates
        RTLSparam = xr.DataArray(
            param, 
            dims=["nparam", "nwav", "nobs"],
            coords={
                "nobs": kernel_wt.coords["nobs"],
                "nwav": kernel_wt.coords["nwav"]
            } 
        )

        # filter for nans
        iGood = ~np.isnan(kernel_wt[0,0,:])
        self.iGood = self.iGood & iGood.values
        self.nobs = np.sum(self.iGood)

        # Pack into SURFACE dataset
        self.SURFACE = xr.Dataset({
            'kernel_wt': kernel_wt,
            'RTLSparam': RTLSparam
        })

