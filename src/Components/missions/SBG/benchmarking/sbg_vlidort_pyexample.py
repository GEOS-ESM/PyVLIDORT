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
from py_vlidort.readers import READERS

from py_vlidort.vlidort import MODIS_BRDF_PMATRIX_run as vl_run
from py_vlidort.vlidort import vl_pack_args_MODIS_BRDF as vl_pack_args
from py_vlidort.vlidort import vl_unpack_result, vl_initOutputs
from py_vlidort.twostream import MODIS_BRDF_run as ts_run
from py_vlidort.twostream import ts_pack_args_MODIS_BRDF as ts_pack_args
from py_vlidort.twostream import ts_unpack_result, ts_initOutputs
from py_vlidort.constants import MISSING
from py_vlidort.io_atts import ARGS_ATTS, OUT_ATTS
from writers import WRITERS

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


class SBG_VLIDORT(INPUTS_VLIDORT,READERS,WRITERS):
    """
    Everything needed for calling VLIDORT
    GEOS-5 has already been sampled on satellite track

    Reads in an SBG granule
    Gets VLIDORT inputs
    Runs VLIDORT

    inFile        : string template for inputs
    outFile       : where to write vlidort outputs
    argsFile      : where to write vlidort inputs
    mtFile        : aerosol optics tables
    albedoType    : what kind of surface albedo model to use
    instname      : instrument name
    nstreams      : number of vlidort streams
    plane_parallel: use plane_parallel assumption in vlidort
    brdfFile      : string template for file with brdf parameters
    verbose       : write debugging outputs
    debug         : write debug files
    """
    def __init__(self,inFile,outFile,argsFile,mtFile,albedoType,
                instname,
                nstreams=12,
                plane_parallel=True,
                brdfFile=None,
                verbose=False,
                debug=False,
                nproc=125):
        self.SDS_AER     = SDS_AER
        self.SDS_MET     = SDS_MET
        self.SDS_INV     = SDS_INV
        self.SDS_ANG     = SDS_ANG
        self.AERNAMES    = AERNAMES
        self.ncALIAS     = ncALIAS
        self.MISSING     = MISSING

        for name, value in locals().items():
            if name != "self":
                setattr(self,name,value)

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

        # Return if no good obs
        if self.nobs == 0:
            return        

        # Read and filter data
        # make some choices here specific to bencharking

        # Read in surface data
        self.readSampledAMESBRDF()

        # Land-Sea Mask
        # limit iGood to land pixels
        self.LandSeaMask()      
        self.iGood = self.iGood & self.iLand
        self.nobs = np.sum(self.iGood)

        # Do only 1 batch of  pixels
        self.nobs = self.nbatch
        self.iGood[np.where(self.iGood)[0][self.nobs:]] = False

        # Read in model data
        self.readSampledGEOS() 

        # Calculate P,T atmospheric profile properties needed for Rayleigh calc
        self.getEdgeVars()

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

        self.writeArgs(ich,sob,eob,args,'RTLS')

        return args

    #---
    def runTWOSTREAM(self,p,ich,sob,args):
        """
        Calls TWOSTREAM
        """

        # get wavlength
        channel = [self.channels[ich]]

        # get last index and number of pts
        eob = min([self.nobs,sob + self.nbatch])
        npts = eob - sob

        # create list of input arguments
        packed_args = ts_pack_args(self,channel,args,npts)

        # run vlidort distributed across processors
        result = p.map(ts_run,packed_args)

        # unpack result
        I, reflectance = ts_unpack_result(result)

        # store outputs
        self.ts_I[sob:eob,ich] = np.squeeze(I)
        self.ts_reflectance[sob:eob,ich] = np.squeeze(reflectance)

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

        # vector or scalar
        NSTOKES = 3
 
        # create list of input arguments
        packed_args = vl_pack_args(self,channel,args,NSTOKES,npts)

        # run vlidort distributed across processors
        result = p.map(vl_run,packed_args)

        # unpack result
        I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U = vl_unpack_result(result)

        # store outputs 
        self.I[sob:eob,ich] = np.squeeze(I)
        self.reflectance[sob:eob,ich] = np.squeeze(reflectance)
        self.surf_reflectance[sob:eob,ich] = np.squeeze(surf_reflectance)
        self.Q[sob:eob,ich] = np.squeeze(Q)
        self.U[sob:eob,ich] = np.squeeze(U) 
        self.BR_Q[sob:eob,ich] = np.squeeze(BR_Q)
        self.BR_U[sob:eob,ich] = np.squeeze(BR_U)
 
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

    parser.add_argument("--do_sphericity",action="store_true",
                        help="use spherical atmosphere in VLIDORT")

    parser.add_argument("-D","--DT_mins", default=DT_mins, type=int,
                        help="Timestep in minutes for each file (default=%i)"%DT_mins)

    parser.add_argument("-n","--nproc", default=nproc, type=int,
                        help="Number of processors to use (default=%i)"%nproc)

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("--debug",action="store_true",
                        help="Debug mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")

    parser.add_argument("--novlidort",action="store_true",
                        help="don't do vlidort calc, aops only (default=False).")

    args = parser.parse_args()
    mtFile         = args.mtFile
    albedoType     = args.albedotype
    do_vlidort     = not args.novlidort

    if args.do_sphericity:
        plane_parallel = False
    else:
        plane_parallel = True


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
        print('>>>plane_parallel',plane_parallel)
        print('>>>nproc:     ',args.nproc)
        print('++++End of arguments+++')
        print('') 
        vlidort = SBG_VLIDORT(inFile,outFile,argsFile,mtFile,
                            albedoType, 
                            instname,
                            brdfFile=brdfFile,
                            verbose=args.verbose,
                            debug=args.debug,
                            plane_parallel=plane_parallel
                            nproc=args.nproc)

        # Run VLIDORT
        # ------------
        if vlidort.nobs == 0:
            print('No valid pixels found. Nothing to do.')
            return  # or sys.exit(0), depending on script structure

        if args.dryrun:
            print('Dry run enabled. Skipping execution.')
            return


        # Initiate Output Arrays
        vl_initOutputs(vlidort)
        ts_initOutputs(vlidort)

        # Get the index of good obs
        iGood = np.where(vlidort.iGood)[0]

        with Pool(vlidort.nproc) as p:

            # Loop through channels
            # only doing the first channel for the benchmark
            for ich,channel in enumerate(vlidort.channels[0:1]):
                print(f'ich: {ich}  channel: {channel}')

                # Get Rayleigh optical depth profile
                vlidort.getROT(channel)

                # loop through nobs in batches
                # just doing one batch for benchmark
                for sob in range(0,vlidort.nobs,vlidort.nbatch):
                    print(f'sob: {sob}, nobs: {vlidort.nobs}')

                    eob = min([vlidort.nobs, sob + vlidort.nbatch])
                    iobs = iGood[sob:eob]
                    npts = eob - sob

                    # Subset inputs for batch
                    # And get optical property inputs
                    args  = vlidort.getargs(ich,sob,eob,iobs,npts)

                    print('   - run TWOSTREAM')
                    vlidort.runTWOSTREAM(p,ich,sob,args)

                    print('   - run VLIDORT')
                    vlidort.runVLIDORT(p,ich,sob,args)

            # Write outputs
            vlidort.writeNC()


        date += Dt
