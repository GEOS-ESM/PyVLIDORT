"""
Helper writing routines
"""

import os
import numpy as np
from py_vlidort.io_atts import ARGS_ATTS, OUT_ATTS
import xarray as xr
import pandas as pd
import shutil
from netCDF4 import Dataset as ncDataset
from numcodecs import Blosc
from glob import glob

class WRITERS(object):
    def expand_dimensions(self):
            """
            Reads the final written file with the 'nobs' dimension, 
            reshapes it back to (time, nalong, ncross), and overwrites.
            """
            print(f"Expanding dimensions for {self.argsFile}...")
            temp_file = self.argsFile + '.tmp'

            # Load the completed dataset
            with xr.open_dataset(self.argsFile) as ds:

                # build the lists of dimensions to use
                dim_ranges = [range(self.ntyme)]
                dim_names = ['time']
                
                if getattr(self, 'nalong', 1) > 1:
                    dim_ranges.append(range(self.nalong))
                    dim_names.append('nalong')
                    
                if getattr(self, 'ncross', 1) > 1:
                    dim_ranges.append(range(self.ncross))
                    dim_names.append('ncross')


                # Create the appropriate index (MultiIndex if >1 dim, regular Index if just time)
                if len(dim_names) > 1:
                    full_index = pd.MultiIndex.from_product(dim_ranges, names=dim_names)
                else:
                    full_index = pd.Index(dim_ranges[0], name=dim_names[0])
                    
                # Filter the index using iGood so it matches the length of the current 'nobs'
                valid_index = full_index[self.iGood]
                
                # Assign the index to 'nobs'
                ds['nobs'] = valid_index

                # Reconstruct the grid: unstack if 2D/3D, rename if 1D
                if len(dim_names) > 1:
                    ds_expanded = ds.unstack('nobs')
                else:
                    ds_expanded = ds.rename({'nobs': dim_names[0]})


                # Save it back out with the expanded dimensions
                ds_expanded.to_netcdf(temp_file, format='NETCDF4')

            # Replace the original flat file with the new expanded file
            shutil.move(temp_file, self.argsFile)
            

    #---
    def writeArgs(self,ich,sob,eob,iobs,args,surface_type):
        """
        Write out the input arguments
        """


        os.makedirs(os.path.dirname(self.argsFile), exist_ok=True)

        # Atmospheric Optical Property Inputs
        if surface_type == 'RTLS':
            rot,depol_ratio,alpha,tau,ssa,g,pmatrix,tauI,ssaI,gI,pmatrixI,tauL,ssaL,gL,pmatrixL,pe,te,ze,param,kernel_wt,vza,sza,raa,flux_factor = args
        elif surface_type == 'LAMBERTIAN':
            rot,depol_ratio,alpha,tau,ssa,g,pmatrix,tauI,ssaI,gI,pmatrixI,tauL,ssaL,gL,pmatrixL,pe,te,ze,vza,sza,raa,flux_factor,albedo = args

        # Number of obs in this batch
        nbatch = eob - sob

        # Define global dataset coordinates
        coords = {
            'lev':     np.arange(self.nlev),
            'leve':    np.arange(self.nlev + 1),
            'nobs':    np.arange(nbatch),
            'ang':     np.arange(self.ang),
            'ch':      [ich],
            'nkernel': np.arange(3),
            'nparam':  np.arange(2),
            'npol':    np.arange(6),
        }

        # Map variable names to (dimensions, data, attributes)
        data_vars = {

            # 1-D Variables
            'VZA':            (['nobs'], vza, ARGS_ATTS['vza']),
            'SZA':            (['nobs'], sza, ARGS_ATTS['sza']),
            'RAA':            (['nobs'], raa, ARGS_ATTS['raa']),
            'IOBS':           (['nobs'], iobs, {'long_name': 'Original indices of valid observations'}), 

            # 2-D Variables
            'PE':             (['leve', 'nobs'], pe, ARGS_ATTS['pe']),
            'TE':             (['leve', 'nobs'], te, ARGS_ATTS['te']),
            'ZE':             (['leve', 'nobs'], ze, ARGS_ATTS['ze']),

            # 3-D Variables
            'ROT':            (['lev', 'nobs', 'ch'], rot, ARGS_ATTS['rot']),
            'TAU':            (['lev', 'nobs', 'ch'], tau.transpose(0, 2, 1), ARGS_ATTS['tau']),
            'SSA':            (['lev', 'nobs', 'ch'], ssa.transpose(0, 2, 1), ARGS_ATTS['ssa']),
            'G':              (['lev', 'nobs', 'ch'], g.transpose(0, 2, 1), ARGS_ATTS['g']),

            # 5-D Variable
            'PMATRIX':        (['lev', 'ang', 'npol', 'nobs', 'ch'], pmatrix.transpose(0, 3, 4, 2, 1), ARGS_ATTS['pmatrix']),
        }

        if surface_type == 'RTLS':
            data_vars['RTLS_PARAM'] = (['nparam', 'nobs', 'ch'], param.transpose(0,2,1), ARGS_ATTS['param'])
            data_vars['RTLS_KERNEL_WT'] = (['nkernel', 'nobs', 'ch'], kernel_wt.transpose(0,2,1), ARGS_ATTS['kernel_wt'])
        elif surface_type == 'LAMBERTIAN':
            albedo_array = np.full((self.nbatch, 1), albedo)
            data_vars['ALBEDO'] = (['nobs', 'ch'], albedo_array, ARGS_ATTS['albedo'])

        # Variables WITHOUT nobs dimension (only written once)
        static_vars = {
            'wavelength': (['ch'], [self.channels[ich]], ARGS_ATTS['ch']),
            'angle':      (['ang'], self.angles.data, ARGS_ATTS['ang']),
            'DEPOL_RATIO': (['ch'], depol_ratio, ARGS_ATTS['depol']),
        }



        # create xarray data arrays of inputs
        # ------------------------------------
        if sob == 0:        
            all_vars = {**data_vars, **static_vars}
        else:
            all_vars = data_vars

        # Create dataset
        ds = xr.Dataset(
            data_vars=all_vars,
            coords=coords,
            attrs={
                'title': 'VLIDORT-GEOS-SBG Simulator Input Arguments',
                'institution': 'NASA/Goddard Space Flight Center',
                'source': 'Global Model and Assimilation Office',
                'history': 'VLIDORT inputs derived from a GEOS simulation',
                'references': 'n/a',
                'contact': 'Patricia Castellanos <patricia.castellanos@nasa.gov>',
                'Conventions': 'CF',
                'inFile': self.inFile,
            }
        )


        # Generate encoding dictionary for everything except wavelength & angle
        compressor = Blosc(cname='lz4', clevel=1, shuffle=Blosc.BITSHUFFLE)
        vars_to_encode = [k for k in data_vars.keys() if k not in ('wavelength', 'angle')]
        var_encoding = {"compressor": compressor, "dtype": "f4"}

        encoding = {var: var_encoding for var in vars_to_encode}

        # chunk data
        chunk_dict = {
            'nobs': 500, 
            'ang': 371, 
            'lev': 72, 
            'leve': 73, 
            'npol': 6, 
            'nkernel': 3, 
            'nparam': 2,
        }
        ds_chunked = ds.chunk({k: v for k, v in chunk_dict.items() if k in ds.dims})
        # Remove conflicting missing_value from attrs
        for var in ds_chunked.data_vars:
            if 'missing_value' in ds_chunked[var].attrs:
                del ds_chunked[var].attrs['missing_value']

        channel_file = self.argsFile.replace('.zarr', f'_ch{ich:03d}.zarr')

        if sob == 0:
            # First batch: create new file
            ds_chunked.to_zarr(channel_file, mode='w', encoding=encoding)
        else:
            # Subsequent batches: append along nobs
            ds_chunked.to_zarr(channel_file, append_dim='nobs')

    #---
    def writeNC (self, ich_start=None, ich_end=None):
        """
        Write a NetCDF file vlidort output
        """

        os.makedirs(os.path.dirname(self.outFile),exist_ok=True)

        # Determine channel range
        if ich_start is None:
            ich_start = 0
        if ich_end is None:
            ich_end = self.nch

        # Subset to only processed channels
        ch_slice = slice(ich_start, ich_end)
        nch_out = ich_end - ich_start
        channels_out = self.channels[ich_start:ich_end]

        # key -> (Data source, OUT_ATTS key)
        variables_map = {
            'toa_reflectance':    (self.reflectance[:, ch_slice],    'toa'),
            'I':                  (self.I[:, ch_slice],              'I'),
            'Q':                  (self.Q[:, ch_slice],              'Q'),
            'U':                  (self.U[:, ch_slice],              'U'),
            'surf_reflectance':   (self.surf_reflectance[:, ch_slice], 'sref'),
            'surf_reflectance_Q': (self.BR_Q[:, ch_slice],           'srefq'),
            'surf_reflectance_U': (self.BR_U[:, ch_slice],           'srefu'),
            'ts_toa_reflectance': (self.ts_reflectance[:, ch_slice], 'ts_toa'),
            'ts_I':               (self.ts_I[:, ch_slice],           'ts_I'),
        }

        # build data_vars dictionary
        data_vars = {
            name: (["nobs", "ch"], data, OUT_ATTS[attr_key])
            for name, (data, attr_key) in variables_map.items()
        }
        data_vars['wavelength'] = (['ch'], channels_out, OUT_ATTS['ch'])

        #  create output dataset
        ds = xr.Dataset(
            data_vars=data_vars,
            coords={
                'ch': np.arange(nch_out)
            },
            attrs={
                'title': 'VLIDORT-GEOS-SBG Simulator',
                'institution': 'NASA/Goddard Space Flight Center',
                'source': 'Global Model and Assimilation Office',
                'history': 'VLIDORT simulation run on sampled GEOS',
                'references': 'n/a',
                'contact': 'Patricia Castellanos <patricia.castellanos@nasa.gov>',
                'Conventions': 'CF',
                'inFile': self.inFile,
                'ich_start': ich_start,
                'ich_end': ich_end,
            }
        )

        # expand nobs back to original dimensions
        dim_ranges = [range(self.ntyme)]
        dim_names = ['time']
        
        if getattr(self, 'nalong', 1) > 1:
            dim_ranges.append(range(self.nalong))
            dim_names.append('nalong')
            
        if getattr(self, 'ncross', 1) > 1:
            dim_ranges.append(range(self.ncross))
            dim_names.append('ncross')

        # Create the appropriate index (MultiIndex if >1 dim, regular Index if just time)
        if len(dim_names) > 1:
            full_index = pd.MultiIndex.from_product(dim_ranges, names=dim_names)
        else:
            full_index = pd.Index(dim_ranges[0], name=dim_names[0])

        # Filter the index using iGood so it matches the length of the current 'nobs'
        valid_index = full_index[self.iGood]

        # Assign the index to 'nobs'
        ds['nobs'] = valid_index

        # Reconstruct the grid: unstack if 2D/3D, rename if 1D
        if len(dim_names) > 1:
            ds = ds.unstack('nobs')
        else:
            ds = ds.rename({'nobs': dim_names[0]})

        # Generate encoding dict
        var_encoding = {"zlib": True, "_FillValue": self.MISSING, "dtype": "f4"}
        encoding = {name: var_encoding for name in variables_map.keys()}

        # Write netcdf file
        ds.to_netcdf(
            path=self.outFile,
            format='NETCDF4',
            engine='netcdf4',
            encoding=encoding,
            unlimited_dims=['ch']
        )

        if self.verbose:
            print(" <> wrote %s "%(self.outFile))
