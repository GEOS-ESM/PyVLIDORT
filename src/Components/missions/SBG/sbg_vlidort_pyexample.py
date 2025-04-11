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
import yaml
from inputs_vlidort import INPUTS_VLIDORT


from py_vlidort.vlidort import MODIS_BRDF_run
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

    Reads in an SBG granule
    Gets VLIDORT inputs
    Runs VLIDORT

    inFile        : string template for inputs
    outFile       : where to write vlidort outputs
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
    def __init__(self,inFile,outFile,mtFile,albedoType,
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

        if self.nobs == 0:
            self.writeNC(empty=True)
        else:

            # Read in model data
            self.readSampledGEOS()

            # Read in surface data
            self.readSampledAMESBRDF()

            # Calculate P,T atmospheric profile properties needed for Rayleigh calc
            self.getEdgeVars()   

            # Land-Sea Mask
            # limit iGood to land pixels
            self.LandSeaMask()       

            # Get pool of processors
            if not dryrun and do_vlidort:
                p = Pool(self.nproc)

            if self.nobs == 0:
                self.writeNC(empty=True)
            else:
                if not dryrun:
                    # Loop through channels
                    for ich,wavelength in enumerate(self.channels[0:1]):
                        # Get Rayleigh optical depth profile
                        self.getROT(wavelength)

                        # Initiate Output Arrays
                        self.initOutputs()
                        if do_vlidort:
                            # Run VLIDORT using multiprocessing
                            self.runVLIDORT(p,ich)

                        # Write outputs
                        self.writeNC(ich)
 
            # close pool of processors
            if not dryrun and do_vlidort:
                p.close()

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


    #---
    def readSampledAMESBRDF(self):
        """
        Read in AMES BRDF kernel weights
        that have already been sampled on swath
        """
        if self.verbose:
            print('opening BRDF file ',self.brdfFile)
        ds = xr.open_dataset(self.brdfFile,chunks="auto")

        # concatenate kernel weights into one data array
        da = xr.concat([ds.Ki,ds.Kg,ds.Kv],dim="nalong")
        da = da.rename("kernel_wt")

        # reshape to [nkernel,nch,nobs]
        self.kernel_wt = da.squeeze().stack(nobs=("time","ncross")).transpose("nalong","nwav","nobs")

        # make an array of the params, just big enough for a batch
        param1 = np.array([2]*self.nbatch)
        param2 = np.array([1]*self.nbatch)

        #[nparam,nch,nobs]
        param1.shape = (1,1,self.nbatch)
        param2.shape = (1,1,self.nbatch)

        # [nparam,nch,nobs]
        self.RTLSparam = np.append(param1,param2,axis=0)
        self.RTLSparam = self.RTLSparam.reshape(2,1,self.nbatch)

        # filter for nans
        iGood = ~np.isnan(self.kernel_wt[0,0,:])
        self.iGood = self.iGood & iGood.values
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
        self.SSA  = np.ones([nlev,nobs,nch])*MISSING
        self.TAU  = np.ones([nlev,nobs,nch])*MISSING
        self.G    = np.ones([nlev,nobs,nch])*MISSING

    #---
    def runVLIDORT(self,p,ich):
        """
        Calls VLIDORT 
        """

        # Get the index of good obs
        iGood  = np.arange(len(self.iGood))[self.iGood]
        nobs   = self.nobs

        # get wavlength
        channel = [self.channels[ich]]

        # loop through nobs in batches
        for sob in range(0,self.nobs,self.nbatch):
            print('sob, nobs',sob, self.nobs)
            eob = min([self.nobs,sob + self.nbatch])
            iobs = iGood[sob:eob]
            npts = eob - sob

            # Subset ROT for good obs only. dims are [nlev,nobs]
            rot  = self.ROT[:,iobs,:]
            depol_ratio = self.depol_ratio

            # Subset aerosol fields for good obs only.
            # calculate AOPs dims are [nlev,nch,nobs]
            self.aer = self.AER.isel(nobs=iobs)    
            self.getpyobsAOP(self.channels[ich])
            tau = self.tau
            ssa = self.ssa
            pmom = self.pmom

            # Subset vertical levels for good obs only. dims are [nlev+1,nobs]
            pe   = self.pe[:,iobs].astype('float64')
            ze   = self.ze[:,iobs].astype('float64')
            te   = self.te[:,iobs].astype('float64')


            # Surbset surface data for good obs only. dims are [nparam,nch,nobs]
            param     = self.RTLSparam[:,:,0:npts].astype('float64')
            kernel_wt = self.kernel_wt[:,ich:ich+1,iobs].astype('float64')

            # subset angles for good obs only. dims are [nobs]
            vza = self.VZA[iobs].astype('float64')
            sza = self.SZA[iobs].astype('float64')
            raa = self.RAA[iobs].astype('float64')
           
            # initialize temporary outputs 
            I = []
            Q = []
            U = []
            reflectance = []
            surf_reflectance = []
            BR_Q = []
            BR_U = []

            # create list of input arguments
            args = [(channel, self.nstreams, self.plane_parallel, rot[:,i:i+1,:], depol_ratio, 
                    tau[:,:,i:i+1], ssa[:,:,i:i+1], pmom[:,:,i:i+1,:,:],
                    pe[:,i:i+1], ze[:,i:i+1], te[:,i:i+1],
                    kernel_wt[:,:,i:i+1], param[:,:,i:i+1],
                    sza[i:i+1], raa[i:i+1], vza[i:i+1],
                    MISSING,
                    self.verbose) for i in range(npts)]

            # run vlidort distributed across processors
            result = p.map(MODIS_BRDF_run,args)
           
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
            self.TAU[:,sob:eob,:] = np.transpose(tau,[0,2,1])
            self.SSA[:,sob:eob,:] = np.transpose(ssa,[0,2,1])
            self.G[:,sob:eob,:] = np.transpose(self.g,[0,2,1])

    #---
    def writeNC (self,ich=None,empty=False):
        """
        Write a NetCDF file vlidort output
        """

        if not os.path.exists(os.path.dirname(self.outFile)):
            os.makedirs(os.path.dirname(self.outFile))

        # create xarray data arrays of outputs
        # ------------------------------------
        if ich is not None:
            nch = 1
            channels = self.channels[ich]
            chcoord = [ich]
        else:
            nch = self.nch
            channels = self.channels 
            chcoord = np.arange(nch)

        # use an input file to set dims and coords
        col = 'aer_Nv'
        dsaer = xr.open_dataset(self.inFile.replace('%col',col))

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

        depol_att = dict(
                    standard_name="depolarization ratio",
                    long_name="Rayleigh depolarization ratio",
                    units="None",
                    missing_value=MISSING,
                    )

        rot_att = dict(
                    standard_name="Rayleigh optical thickness",
                    long_name="Rayleigh optical thickness",
                    units="None",
                    missing_value=MISSING,
                    )

        ch_att = dict(
                    standard_name="wavelength",
                    long_name="wavelength",
                    units="nm",
                    missing_value=MISSING,
                    )

        aot_att = dict(
                    standard_name="Aerosol Optical Thickness",
                    long_name="Aerosol Optical Thickness",
                    units="None",
                    missing_value=MISSING,
                    )

        ssa_att = dict(
                    standard_name="Aerosol Single Scattering Albedo",
                    long_name="Aerosol Single Scattering Albedo",
                    units="None",
                    missing_value=MISSING,
                    )

        g_att = dict(
                    standard_name="Aerosol Assymetry Parameter",
                    long_name="Aerosol Assymetry Parameter",
                    units="None",
                    missing_value=MISSING,
                    )

        # Create dataset
        ds = xr.Dataset(
            data_vars=dict(
                longitude=(['time','across'],dsaer['longitude'].squeeze().values,dsaer['longitude'].attrs),
                latitude=(['time','across'],dsaer['latitude'].squeeze().values,dsaer['latitude'].attrs),
                wavelength=(['ch'],channels,ch_att),
            ),
            coords=dict(
                time=dsaer['time'],
                lev=dsaer['lev'],
                across=np.arange(self.nacross),
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
        dims = ["time", "across", "ch"]
        coords = dict(
                time=dsaer['time'],
                across=np.arange(self.nacross),
                ch=chcoord,
            )

        temp = np.ones([self.ntyme*self.nacross,nch])*MISSING
        if not empty:
            temp[self.iGood,:] = self.reflectance
        da = xr.DataArray(temp.reshape([self.ntyme,self.nacross,nch]),dims=dims,coords=coords,attrs=toa_att)
        ds['toa_reflectance'] = da

        temp = np.ones([self.ntyme*self.nacross,nch])*MISSING
        if not empty:
            temp[self.iGood,:] = self.I
        da = xr.DataArray(temp.reshape([self.ntyme,self.nacross,nch]),dims=dims,coords=coords,attrs=I_att)
        ds['I'] = da

        temp = np.ones([self.ntyme*self.nacross,nch])*MISSING
        if not empty:
            temp[self.iGood,:] = self.Q
        da = xr.DataArray(temp.reshape([self.ntyme,self.nacross,nch]),dims=dims,coords=coords,attrs=Q_att)
        ds['Q'] = da

        temp = np.ones([self.ntyme*self.nacross,nch])*MISSING
        if not empty:
            temp[self.iGood,:] = self.U
        da = xr.DataArray(temp.reshape([self.ntyme,self.nacross,nch]),dims=dims,coords=coords,attrs=U_att)
        ds['U'] = da

        temp = np.ones([self.ntyme*self.nacross,nch])*MISSING
        if not empty:
            temp[self.iGood,:] = self.surf_reflectance
        da = xr.DataArray(temp.reshape([self.ntyme,self.nacross,nch]),dims=dims,coords=coords,attrs=sref_att)
        ds['surf_reflectance'] = da

        temp = np.ones([self.ntyme*self.nacross,nch])*MISSING
        if not empty:
            temp[self.iGood,:] = self.BR_Q
        da = xr.DataArray(temp.reshape([self.ntyme,self.nacross,nch]),dims=dims,coords=coords,attrs=srefq_att)
        ds['surf_reflectance_Q'] = da

        temp = np.ones([self.ntyme*self.nacross,nch])*MISSING
        if not empty:
            temp[self.iGood,:] = self.BR_U
        da = xr.DataArray(temp.reshape([self.ntyme,self.nacross,nch]),dims=dims,coords=coords,attrs=srefu_att)
        ds['surf_reflectance_U'] = da

        temp = np.ones([self.ntyme*self.nacross,nch])*MISSING
        if not empty:
            temp[self.iGood,:] = self.depol_ratio[:]
        da = xr.DataArray(temp.reshape([self.ntyme,self.nacross,nch]),dims=dims,coords=coords,attrs=depol_att)
        ds['depol_ratio'] = da

        # Add 3-D varaibles
        dims = ["lev","time", "across", "ch"]
        coords=dict(
                time=dsaer['time'],
                lev=dsaer['lev'],
                across=np.arange(self.nacross),
                ch=chcoord,
            )

        temp = np.ones([self.nlev,self.ntyme*self.nacross,nch])*MISSING
        if not empty:
            temp[:,self.iGood,:] = self.ROT
        da = xr.DataArray(temp.reshape([self.nlev,self.ntyme,self.nacross,nch]),dims=dims,coords=coords,attrs=rot_att)
        ds['ROT'] = da

        temp = np.ones([self.nlev,self.ntyme*self.nacross,nch])*MISSING
        if not empty:
            temp[:,self.iGood,:] = self.TAU
        da = xr.DataArray(temp.reshape([self.nlev,self.ntyme,self.nacross,nch]),dims=dims,coords=coords,attrs=aot_att)
        ds['TAU'] = da

        temp = np.ones([self.nlev,self.ntyme*self.nacross,nch])*MISSING
        if not empty:
            temp[:,self.iGood,:] = self.SSA
        da = xr.DataArray(temp.reshape([self.nlev,self.ntyme,self.nacross,nch]),dims=dims,coords=coords,attrs=ssa_att)
        ds['SSA'] = da

        temp = np.ones([self.nlev,self.ntyme*self.nacross,nch])*MISSING
        if not empty:
            temp[:,self.iGood,:] = self.G
        da = xr.DataArray(temp.reshape([self.nlev,self.ntyme,self.nacross,nch]),dims=dims,coords=coords,attrs=g_att)
        ds['G'] = da

        # Write netcdf file
        encoding = dict(
            longitude={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            latitude={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            toa_reflectance={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            I={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            Q={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            U={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            surf_reflectance={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            surf_reflectance_Q={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            surf_reflectance_U={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            depol_ratio={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            ROT={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            TAU={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            SSA={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            G={"zlib":True,"_FillValue":MISSING,"dtype":"f4",},
            )

        if (ich == 0) or (ich is None):
            mode='w'
        else:
            mode='a'
        ds.to_netcdf(path=self.outFile,format='NETCDF4',engine='netcdf4',encoding=encoding,unlimited_dims=['ch'],mode=mode)

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
        albedoType = 'AMES_BRDF'


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

        if brdfTemplate is None:
            brdfFile = None
        else:
            brdfFile = brdfTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%minute',minute).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)



        # Initialize VLIDORT class getting aerosol optical properties
        # -----------------------------------------------------------
        print('++++Running VLIDORT with the following arguments+++')
        print('>>>inFile:    ',inFile)
        print('>>>outFile:   ',outFile)
        print('>>>mtFile:    ',mtFile)
        print('>>>albedoType:',albedoType)
        print('>>>brdfFile:  ',brdfFile)
        print('>>>verbose:   ',args.verbose)
        print('>>>nproc:     ',args.nproc)
        print('++++End of arguments+++')
        
        vlidort = SBG_VLIDORT(inFile,outFile,mtFile,
                            albedoType, 
                            instname,
                            args.dryrun,
                            brdfFile=brdfFile,
                            verbose=args.verbose,
                            do_vlidort=do_vlidort,
                            nproc=args.nproc)


        date += Dt
