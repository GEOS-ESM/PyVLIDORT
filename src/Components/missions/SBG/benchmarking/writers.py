"""
Helper writing routines
"""

import os
import numpy as np
from py_vlidort.io_atts import ARGS_ATTS, OUT_ATTS
import xarray as xr

class WRITERS(object):

    #---
    def writeArgs(self,ich,sob,eob,args,surface_type):
        """
        Write out the input arguments
        """


        os.makedirs(os.path.dirname(self.argsFile), exist_ok=True)

        # Atmospheric Optical Property Inputs
        if surface_type == 'RTLS':
            rot,depol_ratio,alpha,tau,ssa,g,pmatrix,tauI,ssaI,gI,pmatrixI,tauL,ssaL,gL,pmatrixL,pe,te,ze,param,kernel_wt,vza,sza,raa,flux_factor = args
        elif surface_type == 'LAMBERTIAN':
            rot,depol_ratio,alpha,tau,ssa,g,pmatrix,tauI,ssaI,gI,pmatrixI,tauL,ssaL,gL,pmatrixL,pe,te,ze,vza,sza,raa,flux_factor,albedo = args

        if (ich == 0) and (sob == 0):
            # Define global dataset coordinates
            coords = {
                'lev':     np.arange(self.nlev),
                'leve':    np.arange(self.nlev + 1),
                'nobs':    np.arange(self.nbatch),
                'ang':     np.arange(self.ang),
                'ch':      [ich],
                'nkernel': np.arange(3),
                'nparam':  np.arange(2),
                'npol':    np.arange(6),
            }

            # Map variable names to (dimensions, data, attributes)
            data_vars = {
                'wavelength':     (['ch'], [self.channels[ich]], ARGS_ATTS['ch']),
                'angle':          (['ang'], self.angles.data, ARGS_ATTS['ang']),

                # 1-D Variables
                'DEPOL_RATIO':    (['ch'], depol_ratio, ARGS_ATTS['depol']),
                'VZA':            (['nobs'], vza, ARGS_ATTS['vza']),
                'SZA':            (['nobs'], sza, ARGS_ATTS['sza']),
                'RAA':            (['nobs'], raa, ARGS_ATTS['raa']),

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
                data_vars['RTLS_PARAM'] = (['nparam', 'nobs', 'ch'], param.transpose(0, 2, 1), ARGS_ATTS['param'])
                data_vars['RTLS_KERNEL_WT'] = (['nkernel', 'nobs', 'ch'], kernel_wt.transpose(0, 2, 1), ARGS_ATTS['kernel_wt'])
            elif surface_type == 'LAMBERTIAN':
                albedo_array = np.full((self.nbatch, 1), albedo)
                data_vars['ALBEDO'] = (['nobs', 'ch'], albedo_array, ARGS_ATTS['albedo'])

            # create xarray data arrays of inputs
            # ------------------------------------

            # Create dataset
            ds = xr.Dataset(
                data_vars=data_vars,
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
            vars_to_encode = [k for k in data_vars.keys() if k not in ('wavelength', 'angle')]
            var_encoding = {"zlib": True, "_FillValue": self.MISSING, "dtype": "f4"}
            encoding = {var: var_encoding for var in vars_to_encode}

            # Write netcdf file
            ds.to_netcdf(
                path=self.argsFile,
                format='NETCDF4',
                engine='netcdf4',
                encoding=encoding,
                unlimited_dims=['nobs', 'ch']
            )

        else:
            #  Append along channel dimension
            nc = ncDataset(self.argsFile,mode='a')

            # 1-D variables
            nc.variables['DEPOL_RATIO'][ich] = depol_ratio

            for name, data in [('VZA', vza), ('SZA', sza), ('RAA', raa)]:
                nc.variables[name][sob:eob] = data

            # 2-D variables
            for name, data in [('PE', pe), ('TE', te), ('ZE', ze)]:
                nc.variables[name][:, sob:eob] = data

            # 3-D variables
            nc.variables['ROT'][:, sob:eob, ich] = rot

            # 3-D variables requiring transpose(0, 2, 1)
            transposed_3d = [
                ('TAU', tau), ('SSA', ssa), ('G', g),
                ('RTLS_PARAM', param), ('RTLS_KERNEL_WT', kernel_wt)
            ]
            for name, data in transposed_3d:
                nc.variables[name][:, sob:eob, ich] = data.transpose(0, 2, 1)

            # 5-D variables
            nc.variables['PMATRIX'][:, :, :, sob:eob, ich] = pmatrix.transpose(0, 2, 3, 4, 1)

            nc.close()

    #---
    def writeNC (self):
        """
        Write a NetCDF file vlidort output
        """

        os.makedirs(os.path.dirname(self.outFile),exist_ok=True)

        # key -> (Data source, OUT_ATTS key)
        variables_map = {
            'toa_reflectance':    (self.reflectance,    'toa'),
            'I':                  (self.I,              'I'),
            'Q':                  (self.Q,              'Q'),
            'U':                  (self.U,              'U'),
            'surf_reflectance':   (self.surf_reflectance, 'sref'),
            'surf_reflectance_Q': (self.BR_Q,           'srefq'),
            'surf_reflectance_U': (self.BR_U,           'srefu'),
            'ts_toa_reflectance': (self.ts_reflectance, 'ts_toa'),
            'ts_I':               (self.ts_I,           'ts_I'),
        }

        # build data_vars dictionary
        data_vars = {
            name: (["nobs", "ch"], data, OUT_ATTS[attr_key])
            for name, (data, attr_key) in variables_map.items()
        }
        data_vars['wavelength'] = (['ch'], self.channels, OUT_ATTS['ch'])

        #  create output dataset
        ds = xr.Dataset(
            data_vars=data_vars,
            coords={
                'nobs': np.arange(self.nobs),
                'ch': np.arange(self.nch)
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
            }
        )

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
