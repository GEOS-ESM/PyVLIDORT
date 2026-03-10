#!/usr/bin/env python3

"""
    Calculates polarized TOA radiance for sbg imager
    Model fields have already been sampled using trj_sampler
    Uses POLAR_VLIDORT as parent class

    Adapted from accp_polar_vlidort.py
    Patricia Castellanos, Jul 2024

    Adapted from polar_vlidort.py and lidar_vlidort.py
    Patricia Castellanos, Jan 2020

"""
import os
import argparse
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
import numpy   as np
import xarray  as xr
from netCDF4 import Dataset as ncDataset
import yaml
from py_vlidort.inputs_vlidort import INPUTS_VLIDORT


from py_vlidort.vlidort import LAMBERTIAN_PMATRIX_run
from py_vlidort.twostream import MODIS_BRDF_run as ts_MODIS_BRDF_run # for now, needs to be updated
from multiprocessing import Pool

# Generic Lists of Varnames and Units
VNAMES_DU = ['DU001','DU002','DU003','DU004','DU005']
VNAMES_SS = ['SS001','SS002','SS003','SS004','SS005']
VNAMES_BC = ['BCPHOBIC','BCPHILIC']
VNAMES_OC = ['OCPHOBIC','OCPHILIC']
VNAMES_SU = ['SO4']

META    = ['DELP','PS','RH','AIRDENS','LONGITUDE','LATITUDE','isotime']
AERNAMES = VNAMES_SU + VNAMES_SS + VNAMES_OC + VNAMES_BC + VNAMES_DU
SDS_AER = META + AERNAMES
SDS_MET = [] #[CLDTOT]
SDS_INV = ['FRLAND']
SDS_ANG = ['SZA','SAA','VZA','VAA']
ncALIAS = {'LONGITUDE': 'longitude',
           'LATITUDE' : 'latitude',
           'SZA'      : 'sza',
           'SAA'      : 'saa',
           'VZA'      : 'vza',
           'VAA'      : 'vaa'}

MISSING = np.float32(-1.e+20)


class SBG_VLIDORT(INPUTS_VLIDORT):
    """
    Everything needed for calling VLIDORT
    GEOS-5 has already been sampled on satellite track

    Reads in an SBG granule, but uses some set lambertian albedos for testing
    Gets VLIDORT inputs
    Runs VLIDORT

    inFile        : string template for inputs
    outFile       : where to write vlidort outputs
    argsFile      : where to write vlidort inputs
    mtFile        : aerosol optics tables
    albedoType    : what kind of surface albedo model to use
    instname      : instrument name
    dryrun        : set everything up but don't calculate AOPs or run VLIDORT
    do_vlidort    : calculate AOPs but don't run VLIDORT
    nstreams      : number of vlidort streams
    plane_parallel: use plane_parallel assumption in vlidort
    brdfFile      : string template for file with brdf parameters
    verbose       : write debugging outputs
    """
    def __init__(self,inFile,outFile,argsFile,mtFile,albedoType,
                instname,dryrun,
                nstreams=12,
                plane_parallel=True,
                brdfFile=None,
                do_vlidort=True,
                verbose=False,
                nproc=125):
        self.SDS_AER     = SDS_AER
        self.SDS_MET     = SDS_MET
        self.SDS_INV     = SDS_INV
        self.SDS_ANG     = SDS_ANG
        self.AERNAMES    = AERNAMES
        self.inFile      = inFile
        self.outFile     = outFile
        self.argsFile    = argsFile
        self.albedoType  = albedoType
        self.mtFile      = mtFile
        self.verbose     = verbose
        self.brdfFile    = brdfFile
        self.nstreams    = nstreams
        self.plane_parallel = plane_parallel
        self.instname   = instname
        self.MISSING    = MISSING
        self.nproc      = nproc


        # load optics tables
        # get:
        # p = number of phase function elements
        # m = number of phase function expansion moments
        self.getMie()

        # get granule dimensions
        self.getDims()

        # get channels
        self.getChannels()
        
        # Start out with all good obs
        self.iGood = np.ones([self.nobs]).astype(bool)

        # Read in precalculated Scene Geometry
        # limit iGood to sza < 80
        self.readAngles()

        if self.nobs > 0:

            # Read in surface data
            self.readSampledAMESBRDF()

            # Land-Sea Mask
            # limit iGood to land pixels
            self.LandSeaMask()      

            # Do only 1 batch of  pixels
            self.nobs = self.nbatch
            self.iGood[np.where(self.iGood)[0][self.nobs:]] = False

            # Read in model data
            self.readSampledGEOS() 

            # Calculate P,T atmospheric profile properties needed for Rayleigh calc
            self.getEdgeVars()

            if (self.nobs > 0) and not dryrun:
                # Initiate Output Arrays
                self.initOutputs()

                # Loop through channels
                for ich,channel in enumerate(self.channels):
                    print('ich: ',ich,' channel: ',channel)
                    # Get Rayleigh optical depth profile
                    self.getROT(channel)

                    if do_vlidort:

                        # Get the index of good obs
                        iGood  = np.arange(len(self.iGood))[self.iGood]

                        # loop through nobs in batches
                        for sob in range(0,self.nobs,self.nbatch):
                            print('sob, nobs',sob, self.nobs)
                            eob = min([self.nobs,sob + self.nbatch])
                            iobs = iGood[sob:eob]
                            npts = eob - sob

                            # Subset inputs for batch
                            # And get optical property inputs
                            args  = self.getargs(ich,sob,eob,iobs,npts)

                            # Get pool of processors
                            p = Pool(self.nproc)
                            # Run TWOSTREAM
                            print('   - run TWOSTREAM')
                            self.runTWOSTREAM(p,ich,sob,args)
                            # Run VLIDORT using multiprocessing
                            print('   - run VLIDORT')
                            self.runVLIDORT(p,ich,sob,args)
                            # close pool of processors
                            p.close()

                # Write outputs
                self.writeNC()
 
    #--
    def getDims(self):
        """
        Get granule dimensions
        """
        col = 'aer_Nv'
        if self.verbose:
            print('opening file',self.inFile.replace('%col',col))
        ds = xr.open_dataset(self.inFile.replace('%col',col)) 
        self.ntyme,self.nlev,self.nacross = ds.sizes['time'],ds.sizes['lev'],ds.sizes['ncross']
        self.nobs = self.nacross*self.ntyme
        self.nbatch = self.nproc

    #---
    def getChannels(self):
        """
        Read in AMES BRDF to get the channels
        """
        if self.verbose:
            print('opening BRDF file ',self.brdfFile)
        ds = xr.open_mfdataset(os.path.dirname(self.brdfFile)+'/*ames*',chunks="auto")

        self.channels = ds.nwav.values  # microns
        self.nch = len(self.channels)
        self.channels = self.channels*1e3  # nm
        ds.close()        

    #---
    def readSampledGEOS(self):
        """
        Read in model sampled track
        """
        col = 'aer_Nv'
        if self.verbose:
            print('opening file',self.inFile.replace('%col',col))

        inList = [self.inFile.replace('%col',col)]

        if len(self.SDS_MET) > 0:
            col = 'met_Nv'
            inList.append(self.inFile.replace('%col',col))
            if self.verbose:
                print('opening file',self.inFile.replace('%col',col))

        self.AER = xr.open_mfdataset(inList,chunks="auto")
        # make arrays [nobs,nlev]
        self.AER = self.AER.squeeze()
        self.AER = self.AER.stack(nobs=("time","ncross"))
        self.AER = self.AER.transpose("nobs","lev")
        iGood = np.arange(len(self.iGood))[self.iGood]
        self.AER = self.AER.isel(nobs=iGood)

    # ---
    def LandSeaMask(self):
        """
        Read in invariant dataset
        """
        col = 'asm_Nx'
        if self.verbose:
            print('opening file',self.inFile.replace('%col',col))
        ds = xr.open_dataset(self.inFile.replace('%col',col),chunks="auto")

        for sds in self.SDS_INV:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = ds[sds_].squeeze().stack(nobs=("time","ncross"))
            self.__dict__[sds] = var

        iGood = self.FRLAND >= 0.99
        self.iLand = iGood.values
        iGood = self.FRLAND < 0.99
        self.iSea  = iGood.values

        # self.iGood = self.iGood & iGood
        self.nobsLand  = np.sum(self.iGood & self.iLand)
        self.nobsSea   = np.sum(self.iGood & self.iSea)

        self.iGood = self.iGood & self.iLand
        self.nobs = np.sum(self.iGood)



    def readAngles(self):
        """
        Read in viewing and solar Geometry from angFile
        """

        col = self.instname
        if self.verbose: 
            print('opening file',self.inFile.replace('%col',col))
        ds = xr.open_dataset(self.inFile.replace('%col',col),chunks="auto")

        for sds in self.SDS_ANG:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = ds[sds_].squeeze().stack(nobs=("time","ncross"))
            self.__dict__[sds] = var

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
    def initOutputs(self):
        # Initiate output arrays
        nch    = self.nch
        nobs   = self.nobs
        nlev   = self.nlev
        self.I = np.ones([nobs,nch])*MISSING
        self.Q = np.ones([nobs,nch])*MISSING
        self.U = np.ones([nobs,nch])*MISSING
        self.reflectance = np.ones([nobs,nch])*MISSING
        self.surf_reflectance = np.ones([nobs,nch])*MISSING
        self.BR_Q = np.ones([nobs,nch])*MISSING
        self.BR_U = np.ones([nobs,nch])*MISSING

        self.ts_I = np.ones([nobs,nch])*MISSING
        self.ts_reflectance = np.ones([nobs,nch])*MISSING

    #---
    def getargs(self,ich,sob,eob,iobs,npts):
        """
        Generic call to subset args
        """
        # Subset ROT for good obs only. dims are [nlev,nobs]
        rot  = self.ROT[:,sob:eob,:]
        depol_ratio = self.depol_ratio

        # Subset aerosol fields for good obs only.
        # calculate AOPs dims are [nlev,nch,nobs]
        self.aer = self.AER.isel(nobs=slice(sob,eob))
        self.getpyobsAOP(self.channels[ich])
        tau  = self.tau
        ssa  = self.ssa
        pmatrix = self.pmatrix
        g    = self.g

        # Subset vertical levels for good obs only. dims are [nlev+1,nobs]
        pe   = self.pe[:,sob:eob].astype('float64')
        ze   = self.ze[:,sob:eob].astype('float64')
        te   = self.te[:,sob:eob].astype('float64')


        # Surbset surface data for good obs only. dims are [nparam,nch,nobs]
        param     = self.RTLSparam[:,:,0:npts].astype('float64')
        kernel_wt = self.kernel_wt[:,ich:ich+1,iobs].astype('float64').to_numpy()

        # subset angles for good obs only. dims are [nobs]
        vza = self.VZA[iobs].astype('float64').to_numpy()
        sza = self.SZA[iobs].astype('float64').to_numpy()
        raa = self.RAA[iobs].astype('float64').to_numpy()

        # cloud properties
        # empty for now
        dims = tau.shape
        tauI, tauL = np.zeros(dims).astype('float64'), np.zeros(dims).astype('float64')
        ssaI, ssaL = np.zeros(dims).astype('float64'), np.zeros(dims).astype('float64')
        gI, gL     = np.zeros(dims).astype('float64'), np.zeros(dims).astype('float64')
        dims = pmatrix.shape
        pmatrixI, pmatrixL = np.zeros(dims).astype('float64'), np.zeros(dims).astype('float64')

        # trace gas absorption
        # empty for now
        dims = tau.shape
        alpha = np.zeros(dims).astype('float64')

        # solar flux (f0)
        # constant for now
        flux_factor = np.ones([1]).astype('float64')

        args = rot,depol_ratio,alpha,tau,ssa,g,pmatrix,tauI,ssaI,gI,pmatrixI,tauL,ssaL,gL,pmatrixL,pe,te,ze,param,kernel_wt,vza,sza,raa,flux_factor

        self.writeArgs(ich,sob,eob,args)

        return args

    #---
    def runTWOSTREAM(self,p,ich,sob,args):
        """
        Calls TWOSTREAM
        """
        # Get pool of processors
        p = Pool(self.nproc)

        # get wavlength
        channel = [self.channels[ich]]

        # get last index and number of pts
        eob = min([self.nobs,sob + self.nbatch])
        npts = eob - sob

        # Atmospheric Optical Property Arguments
        rot,depol_ratio,alpha,tau,ssa,g,pmatrix,tauI,ssaI,gI,pmatrixI,tauL,ssaL,gL,pmatrixL,pe,te,ze,param,kernel_wt,vza,sza,raa,flux_factor = args

        # create list of input arguments
        args = [(channel, self.plane_parallel, rot[:,i:i+1,:], depol_ratio,
                alpha[:,:,i:i+1],
                tau[:,:,i:i+1], ssa[:,:,i:i+1], g[:,:,i:i+1],
                pe[:,i:i+1], ze[:,i:i+1], te[:,i:i+1],
                kernel_wt[:,:,i:i+1], param[:,:,i:i+1],
                [sza[i:i+1]], [raa[i:i+1]], [vza[i:i+1]],
                flux_factor,
                MISSING,
                self.verbose) for i in range(npts)]

        # run vlidort distributed across processors
        result = p.map(ts_MODIS_BRDF_run,args)

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

        # store outputs
        self.ts_I[sob:eob,ich] = np.squeeze(I)
        self.ts_reflectance[sob:eob,ich] = np.squeeze(reflectance)

        # close pool of processors
        p.close()

    #---
    def runVLIDORT(self,p,ich,sob,args):
        """
        Calls VLIDORT 
        """

        # get wavlength
        channel = [self.channels[ich]]

        # get last ob and number of pts
        eob = min([self.nobs,sob + self.nbatch])
        npts = eob - sob

        # Atmospheric Optical Property Inputs
        rot,depol_ratio,alpha,tau,ssa,g,pmatrix,tauI,ssaI,gI,pmatrixI,tauL,ssaL,gL,pmatrixL,pe,te,ze,param,kernel_wt,vza,sza,raa,flux_factor = args
      
        # vector or scalar
        NSTOKES = 3
 
        # create list of input arguments
        args = [(channel, self.nstreams, self.plane_parallel,NSTOKES,self.angles,
                rot[:,i:i+1,:], depol_ratio, alpha[:,:,i:i+1],
                tau[:,:,i:i+1], ssa[:,:,i:i+1], pmatrix[:,:,i:i+1,:,:],
                tauI[:,:,i:i+1], ssaI[:,:,i:i+1], pmatrixI[:,:,i:i+1,:,:],
                tauL[:,:,i:i+1], ssaL[:,:,i:i+1], pmatrixL[:,:,i:i+1,:,:],
                pe[:,i:i+1], ze[:,i:i+1], te[:,i:i+1],
                kernel_wt[:,:,i:i+1], param[:,:,i:i+1],
                [sza[i:i+1]], [raa[i:i+1]], [vza[i:i+1]],
                flux_factor,
                MISSING,
                self.verbose) for i in range(npts)]

        # run vlidort distributed across processors
        result = p.map(MODIS_BRDF_PMATRIX_run,args)

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
            surf_reflectance.append(surf_reflectance_r)
            BR_Q.append(BR_Q_r)
            BR_U.append(BR_U_r)
        I = np.concatenate(I)
        Q = np.concatenate(Q)
        U = np.concatenate(U)
        reflectance = np.concatenate(reflectance)
        surf_reflectance = np.concatenate(surf_reflectance)
        BR_Q = np.concatenate(BR_Q)
        BR_U = np.concatenate(BR_U)

        # store outputs 
        self.I[sob:eob,ich] = np.squeeze(I)
        self.reflectance[sob:eob,ich] = np.squeeze(reflectance)
        self.surf_reflectance[sob:eob,ich] = np.squeeze(surf_reflectance)
        self.Q[sob:eob,ich] = np.squeeze(Q)
        self.U[sob:eob,ich] = np.squeeze(U) 
        self.BR_Q[sob:eob,ich] = np.squeeze(BR_Q)
        self.BR_U[sob:eob,ich] = np.squeeze(BR_U)

    #---
    def writeArgs(self,ich,sob,eob,args):
        """
        Write out the input arguments
        """


        if not os.path.exists(os.path.dirname(self.argsFile)):
            os.makedirs(os.path.dirname(self.argsFile))

        # Atmospheric Optical Property Inputs
        rot,depol_ratio,alpha,tau,ssa,g,pmatrix,tauI,ssaI,gI,pmatrixI,tauL,ssaL,gL,pmatrixL,pe,te,ze,param,kernel_wt,vza,sza,raa,flux_factor = args


        if (ich == 0) and (sob == 0):
            # create xarray data arrays of inputs
            # ------------------------------------

            # Define variable attributes
            ch_att = dict(
                        standard_name="wavelength",
                        long_name="wavelength",
                        units="nm",
                        missing_value=MISSING,
                        )

            ang_att = dict(
                        standard_name="angle",
                        long_name="scattering angle",
                        units="degrees",
                        missing_value=MISSING,
                        )

            rot_att = dict(
                        standard_name="ROT",
                        long_name="rayleigh optical thickenss",
                        units="None",
                        missing_value=MISSING,
                        )

            depol_att = dict(
                        standard_name="depolarization ratio",
                        long_name="Rayleigh depolarization ratio",
                        units="None",
                        missing_value=MISSING,
                        )

            tau_att = dict(
                        standard_name="TAU",
                        long_name="aerosol optical thickenss",
                        units="None",
                        missing_value=MISSING,
                        )

            ssa_att = dict(
                        standard_name="SSA",
                        long_name="single scattering albedo",
                        units="None",
                        missing_value=MISSING,
                        )

            g_att = dict(
                        standard_name="G",
                        long_name="assymetry parameter",
                        units="None",
                        missing_value=MISSING,
                        )

            pmatrix_att = dict(
                        standard_name="PMATRIX",
                        long_name="aerosol scattering matrix (p11,22,33,44,12,34)",
                        units="None",
                        missing_value=MISSING,
                        )

            pe_att = dict(
                        standard_name="PE",
                        long_name="pressure at layer edges",
                        units="Pa",
                        missing_value=MISSING,
                        )

            ze_att = dict(
                        standard_name="ZE",
                        long_name="height above sea level at layer edges",
                        units="m",
                        missing_value=MISSING,
                        )

            te_att = dict(
                        standard_name="TE",
                        long_name="temperature at layer edges",
                        units="K",
                        missing_value=MISSING,
                        )

            kernel_wt_att = dict(
                        standard_name="KERNEL_WT",
                        long_name="RTLS BRDF kernel weights",
                        units="None",
                        missing_value=MISSING,
                        )

            param_att = dict(
                        standard_name="RTLS_PARAM",
                        long_name="RTLS parameters",
                        units="K",
                        missing_value=MISSING,
                        )

            sza_att = dict(
                        standard_name="SZA",
                        long_name="solar zenith angle",
                        units="degrees",
                        missing_value=MISSING,
                        )

            vza_att = dict(
                        standard_name="VZA",
                        long_name="sensor zenith angle",
                        units="degrees",
                        missing_value=MISSING,
                        )

            raa_att = dict(
                        standard_name="RAA",
                        long_name="relative azimuth angle",
                        units="degrees",
                        missing_value=MISSING,
                        )

            flux_att = dict(
                        standard_name="FLUX",
                        long_name="solar flux (F0)",
                        units="Done",
                        missing_value=MISSING,
                        )

            # Create dataset
            ds = xr.Dataset(
                data_vars=dict(
                    wavelength=(['ch'],[self.channels[ich]],ch_att),
                    angle=(['ang'],self.angles.data,ang_att),
                ),
                coords=dict(
                    lev=np.arange(self.nlev),
                    leve=np.arange(self.nlev+1),
                    nobs=np.arange(self.nbatch),
                    ang=np.arange(self.ang),
                    ch=[ich],
                    nkernel=np.arange(3),
                    nparam=np.arange(2),
                    npol=np.arange(6),
                ),
                attrs=dict(
                    title='VLIDORT-GEOS-SBG Simulator Input Arguments',
                    institution = 'NASA/Goddard Space Flight Center',
                    source = 'Global Model and Assimilation Office',
                    history = 'VLIDORT inputs derived from a GEOS simulation',
                    references = 'n/a',
                    contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>',
                    Conventions = 'CF',
                    inFile = self.inFile,
                    ),
            )

            # Add 1-D Variables
            dims = ["ch"]
            coords = dict(
                    ch=[ich],
                )

            da = xr.DataArray(depol_ratio,dims=dims,coords=coords,attrs=depol_att)
            ds['DEPOL_RATIO'] = da


            dims = ["nobs"]
            coords = dict(
                    nobs=np.arange(self.nbatch),
                )

            da = xr.DataArray(vza,dims=dims,coords=coords,attrs=vza_att)
            ds['VZA'] = da

            da = xr.DataArray(sza,dims=dims,coords=coords,attrs=sza_att)
            ds['SZA'] = da

            da = xr.DataArray(raa,dims=dims,coords=coords,attrs=raa_att)
            ds['RAA'] = da

            # Add 2-D Variables
            dims = ["leve", "nobs"]
            coords = dict(
                    nobs=np.arange(self.nbatch),
                    leve=np.arange(self.nlev+1),
                )

            da = xr.DataArray(pe,dims=dims,coords=coords,attrs=pe_att)
            ds['PE'] = da

            da = xr.DataArray(te,dims=dims,coords=coords,attrs=te_att)
            ds['TE'] = da

            da = xr.DataArray(ze,dims=dims,coords=coords,attrs=ze_att)
            ds['ZE'] = da

            # Add 3-D varaibles
            dims = ["lev","nobs","ch"]
            coords=dict(
                    lev=np.arange(self.nlev),
                    nobs=np.arange(self.nbatch),
                    ch=[ich],
                )

            da = xr.DataArray(rot,dims=dims,coords=coords,attrs=rot_att)
            ds['ROT'] = da

            da = xr.DataArray(tau.transpose(0,2,1),dims=dims,coords=coords,attrs=tau_att)
            ds['TAU'] = da

            da = xr.DataArray(ssa.transpose(0,2,1),dims=dims,coords=coords,attrs=ssa_att)
            ds['SSA'] = da

            da = xr.DataArray(g.transpose(0,2,1),dims=dims,coords=coords,attrs=g_att)
            ds['G'] = da      

            dims = ["nparam","nobs","ch"]
            coords=dict(
                    nparam=np.arange(2),
                    nobs=np.arange(self.nbatch),
                    ch=[ich],
                )  

            da = xr.DataArray(param.transpose(0,2,1),dims=dims,coords=coords,attrs=param_att)
            ds['RTLS_PARAM'] = da

            dims = ["nkernel","nobs","ch"]
            coords=dict(
                    nkernel=np.arange(3),
                    nobs=np.arange(self.nbatch),
                    ch=[ich],
                )

            da = xr.DataArray(kernel_wt.transpose(0,2,1),dims=dims,coords=coords,attrs=kernel_wt_att)
            ds['RTLS_KERNEL_WT'] = da

            # Add 5-D varaibles
            dims = ["lev","ang","npol","nobs","ch"]
            coords=dict(
                    lev=np.arange(self.nlev),
                    nobs=np.arange(self.nbatch),
                    ch=[ich],
                    npol=np.arange(6),
                )

            da = xr.DataArray(pmatrix.transpose(0,3,4,2,1),dims=dims,coords=coords,attrs=pmatrix_att)
            ds['PMATRIX'] = da

            # Write netcdf file
            encoding = dict(
                ZE={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
                TE={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
                PE={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
                RAA={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
                VZA={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
                SZA={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
                DEPOL_RATIO={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
                ROT={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
                TAU={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
                SSA={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
                G={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
                RTLS_PARAM={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
                RTLS_KERNEL_WT={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
                PMATRIX={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
                )

            ds.to_netcdf(path=self.argsFile,format='NETCDF4',engine='netcdf4',encoding=encoding,unlimited_dims=['nobs','ch'])

        else:
            #  Append along channel dimension
            nc = ncDataset(self.argsFile,mode='a')

            # add 1-D variables
            var = nc.variables['DEPOL_RATIO']
            var[ich] = depol_ratio

            var = nc.variables['VZA']
            var[sob:eob] = vza

            var = nc.variables['SZA']
            var[sob:eob] = sza

            var = nc.variables['RAA']
            var[sob:eob] = raa

            # Add 2-D variables
            var = nc.variables['PE']
            var[:,sob:eob] = pe

            var = nc.variables['TE']
            var[:,sob:eob] = te

            var = nc.variables['ZE']
            var[:,sob:eob] = ze

            # Add 3-D varaibles
            var = nc.variables['ROT']
            var[:,sob:eob,ich] = rot

            var = nc.variables['TAU']
            var[:,sob:eob,ich] = tau.transpose(0,2,1)

            var = nc.variables['SSA']
            var[:,sob:eob,ich] = ssa.transpose(0,2,1)

            var = nc.variables['G']
            var[:,sob:eob,ich] = g.transpose(0,2,1)

            var = nc.variables['RTLS_PARAM']
            var[:,sob:eob,ich] = param.transpose(0,2,1)

            var = nc.variables['RTLS_KERNEL_WT']
            var[:,sob:eob,ich] = kernel_wt.transpose(0,2,1)

            # Add 5-D varaibles
            var = nc.variables['PMATRIX']
            var[:,:,:,sob:eob,ich] = pmatrix.transpose(0,2,3,4,1)

            nc.close()

    #---
    def writeNC (self):
        """
        Write a NetCDF file vlidort output
        """

        if not os.path.exists(os.path.dirname(self.outFile)):
            os.makedirs(os.path.dirname(self.outFile))

        # create xarray data arrays of outputs
        # ------------------------------------
        nch = self.nch
        channels = self.channels
        chcoord = np.arange(self.nch)

        # Define variable attributes
        toa_att = dict(
                    standard_name="TOA Reflectance",
                    long_name="reflectance at the top of the atmosphere",
                    units="None",
                    missing_value=MISSING,
                    )

        I_att = dict(
                    standard_name="TOA I",
                    long_name="intensity at the top of the atmosphere",
                    units="W m-2 sr-1 nm-1",
                    missing_value=MISSING,
                    )
        Q_att = dict(
                    standard_name="TOA Q",
                    long_name="Q-component of the stokes vector at the top of the atmopshere",
                    units="W m-2 sr-1 nm-1",
                    missing_value=MISSING,
                    )

        U_att = dict(
                    standard_name="TOA U",
                    long_name="U-component of the stokes vector at the top of the atmopshere",
                    units="W m-2 sr-1 nm-1",
                    missing_value=MISSING,
                    )

        sref_att = dict(
                    standard_name="Surface Reflectance",
                    long_name="Bi-Directional Surface Reflectance",
                    units="None",
                    missing_value=MISSING,
                    )
        srefq_att = dict(
                    standard_name="Surface Reflectance Q",
                    long_name="Bi-Directional Surface Reflectance Q",
                    units="None",
                    missing_value=MISSING,
                    )

        srefu_att = dict(
                    standard_name="Surface Reflectance U",
                    long_name="Bi-Directional Surface Reflectance U",
                    units="None",
                    missing_value=MISSING,
                    )

        ch_att = dict(
                    standard_name="wavelength",
                    long_name="wavelength",
                    units="nm",
                    missing_value=MISSING,
                    )

        ts_I_att = dict(
                    standard_name="TWOSTREAM TOA I",
                    long_name="intensity at the top of the atmosphere from TWOSTREAM",
                    units="W m-2 sr-1 nm-1",
                    missing_value=MISSING,
                    )

        ts_toa_att = dict(
                    standard_name="TWOSTREAM TOA Reflectance",
                    long_name="reflectance at the top of the atmosphere from TWOSTREAM",
                    units="None",
                    missing_value=MISSING,
                    )

        # Create dataset
        ds = xr.Dataset(
            data_vars=dict(
                wavelength=(['ch'],channels,ch_att),
            ),
            coords=dict(
                nobs=np.arange(self.nobs),
                ch=chcoord,
            ),
            attrs=dict(
                title='VLIDORT-GEOS-SBG Simulator',
                institution = 'NASA/Goddard Space Flight Center',
                source = 'Global Model and Assimilation Office',
                history = 'VLIDORT simulation run on sampled GEOS',
                references = 'n/a',
                contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>',
                Conventions = 'CF',
                inFile = self.inFile,
                ),
        )


        # Add 2-D Variables
        dims = ["nobs", "ch"]
        coords = dict(
                nobs=np.arange(self.nobs),
                ch=chcoord,
            )

        da = xr.DataArray(self.reflectance,dims=dims,coords=coords,attrs=toa_att)
        ds['toa_reflectance'] = da

        da = xr.DataArray(self.I,dims=dims,coords=coords,attrs=I_att)
        ds['I'] = da

        da = xr.DataArray(self.Q,dims=dims,coords=coords,attrs=Q_att)
        ds['Q'] = da

        da = xr.DataArray(self.U,dims=dims,coords=coords,attrs=U_att)
        ds['U'] = da

        da = xr.DataArray(self.surf_reflectance,dims=dims,coords=coords,attrs=sref_att)
        ds['surf_reflectance'] = da

        da = xr.DataArray(self.BR_Q,dims=dims,coords=coords,attrs=srefq_att)
        ds['surf_reflectance_Q'] = da

        da = xr.DataArray(self.BR_U,dims=dims,coords=coords,attrs=srefu_att)
        ds['surf_reflectance_U'] = da

        da = xr.DataArray(self.ts_reflectance,dims=dims,coords=coords,attrs=ts_toa_att)
        ds['ts_toa_reflectance'] = da

        da = xr.DataArray(self.ts_I,dims=dims,coords=coords,attrs=ts_I_att)
        ds['ts_I'] = da

        # Write netcdf file
        encoding = dict(
            toa_reflectance={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            I={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            Q={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            U={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            surf_reflectance={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            surf_reflectance_Q={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            surf_reflectance_U={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            ts_toa_reflectance={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            ts_I={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            )

        ds.to_netcdf(path=self.outFile,format='NETCDF4',engine='netcdf4',encoding=encoding,unlimited_dims=['ch'])

        if self.verbose:
            print(" <> wrote %d %s "%(ich,self.outFile))
 
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
    DT_mins   = 1
    mtFile     = 'm2_aop.yaml'
    albedoType = None
    nproc      = 125

#   Parse command line options
#   --------------------------

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("inputs_yaml",
                        help="yaml file with input configuration file names")

    parser.add_argument("-a","--albedotype", default=albedoType,
                        help="albedo type keyword. default is to figure out according to channel")

    parser.add_argument("--mtFile",default=mtFile,
                        help="mtFile (default=%s)"%mtFile)

    parser.add_argument("-D","--DT_mins", default=DT_mins, type=int,
                        help="Timestep in minutes for each file (default=%i)"%DT_mins)

    parser.add_argument("-n","--nproc", default=nproc, type=int,
                        help="Number of processors to use (default=%i)"%nproc)

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")

    parser.add_argument("--novlidort",action="store_true",
                        help="don't do vlidort calc, aops only (default=False).")

    args = parser.parse_args()
    mtFile         = args.mtFile
    albedoType     = args.albedotype
    do_vlidort     = not args.novlidort

    # figure out albedoType keyword
    if albedoType is None:
        albedoType = 'LAMBERTIAN'


    config = yaml.safe_load(open(args.inputs_yaml))
    args.paths_yaml = config['paths_yaml']
    args.inst_yaml  = config['inst_yaml']
    args.orbit_yaml = config['orbit_yaml']
    # Parse prep config
    # -----------------
    cf             = yaml.safe_load(open(args.inst_yaml))
    instname       = cf['instname']

    cf             = yaml.safe_load(open(args.orbit_yaml))
    orbitname      = cf['orbitname']
    ORBITNAME      = orbitname.upper()

    cf             = yaml.safe_load(open(args.paths_yaml))
    inTemplate     = cf['inDir']     + '/' + cf['inFile']
    outTemplate    = cf['outDir']    + '/' + cf['outFile']

    try:
        brdfTemplate = cf['brdfDir'] + '/' + cf['brdfFile']
    except:
        brdfTemplate = None


    # Loop through dates, running VLIDORT
    # ------------------------------------
    date      = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)
    Dt        = timedelta(minutes=args.DT_mins)

    while date < enddate:
        nymd  = str(date.date()).replace('-','')
        year  = str(date.year)
        month = str(date.month).zfill(2)
        day   = str(date.day).zfill(2)
        hour  = str(date.hour).zfill(2)
        minute = str(date.minute).zfill(2)

        inFile     = inTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%minute',minute).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        outFile    = outTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%minute',minute).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME).replace('%instname',instname)
        argsFile    = outTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%minute',minute).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME).replace('%instname',instname).replace('vlidort','vlidort_args')

        if brdfTemplate is None:
            brdfFile = None
        else:
            brdfFile = brdfTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%minute',minute).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)



        # Initialize VLIDORT class getting aerosol optical properties
        # -----------------------------------------------------------
        print('++++Running VLIDORT with the following arguments+++')
        print('>>>inFile:    ',inFile)
        print('>>>outFile:   ',outFile)
        print('>>>argsFile:  ',argsFile)
        print('>>>mtFile:    ',mtFile)
        print('>>>albedoType:',albedoType)
        print('>>>brdfFile:  ',brdfFile)
        print('>>>verbose:   ',args.verbose)
        print('>>>nproc:     ',args.nproc)
        print('++++End of arguments+++')
        
        vlidort = SBG_VLIDORT(inFile,outFile,argsFile,mtFile,
                            albedoType, 
                            instname,
                            args.dryrun,
                            brdfFile=brdfFile,
                            verbose=args.verbose,
                            do_vlidort=do_vlidort,
                            nproc=args.nproc)


        date += Dt
