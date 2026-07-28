#!/usr/bin/env python3

"""
    Calculates polarized TOA radiance for sbg imager
    Model fields have already been sampled using trj_sampler
    Uses POLAR_VLIDORT as parent class

    Adapted from sbg_vlidort_pyexample.py
    Patricia Castellanos, Jul 2026

    Adapted from accp_polar_vlidort.py
    Patricia Castellanos, Jul 2024

    Adapted from polar_vlidort.py and lidar_vlidort.py
    Patricia Castellanos, Jan 2020

"""
import psutil
import time
import os,sys
import argparse
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
import numpy   as np
import xarray  as xr
import yaml

from py_vlidort.inputs_vlidort import INPUTS_VLIDORT
from py_vlidort.readers import READERS
from py_vlidort.vlidort import vl_initOutputs, VLIDORT
from py_vlidort.twostream import ts_initOutputs, TWOSTREAM
from py_vlidort.constants import MISSING
from py_vlidort.writers import WRITERS

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


class SBG_VLIDORT(INPUTS_VLIDORT,READERS,WRITERS,VLIDORT,TWOSTREAM):
    """
    Everything needed for calling VLIDORT
    GEOS-5 has already been sampled on satellite track

    Reads in an SBG granule
    Gets VLIDORT inputs

    inFile        : string template for inputs
    outFile       : where to write vlidort outputs
    argsFile      : where to write vlidort inputs
    mtFile        : aerosol optics tables
    NSTOKES       : number stokes vector elements
    albedoType    : list of surface model types to use
    instname      : instrument name
    nstreams      : number of vlidort streams
    plane_parallel: use plane_parallel assumption in vlidort
    do_fullrad    : do SS+MS, if false do SS only
    brdfFile      : string template for file with brdf parameters
    verbose       : write debugging outputs
    debug         : write debug files
    do_inputscacl : do input calculation, if False reads from file
    """
    def __init__(self,inFile,outFile,argsFile,mtFile,albedoType,
                instname,
                do_inputscalc=True,
                NSTOKES=3,
                albedo=None,
                nstreams=12,
                plane_parallel=True,
                do_fullrad=True,
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
        # make some choices here specific to SBG

        # Read in surface data            
        # Get the string name of the function from the inherited map
        method_name = self.ALBEDO_READER_MAP.get(albedoType)
            
        if method_name:
            # Get the actual function object using getattr
            reader_func = getattr(self, method_name)
            
            # Call the function
            reader_func() 
        else:
            print(f"Warning: No reader mapped for albedo type '{atype}'")

        # Land-Sea Mask
        # limit iGood to land pixels
        self.LandSeaMask()      
        self.iGood = self.iGood & self.iLand
        self.nobs = np.sum(self.iGood)

        # Read in model data
        self.readSampledGEOS() 

        # Calculate P,T atmospheric profile properties needed for Rayleigh calc
        self.getEdgeVars()

        # Get Rayleigh optical depth profile
        self.getROT(self.channels)

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
    def getargs(self,ich):
        """
        Generic call to subset args
        """
        # Subset ROT for good obs only. dims are [nlev,nobs]
        rot  = self.rayleigh.ROT[:,:,ich:ich+1].values
        depol_ratio = self.rayleigh.depol_ratio[ich:ich+1].values

        # calculate AOPs dims are [nlev,nch,nobs]
        if self.do_inputscalc:
            t1 = time.perf_counter()
            self.getAOP(self.channels[ich])
            tau  = self.tau
            ssa  = self.ssa
            pmatrix = self.pmatrix
            g    = self.g
            t2 = time.perf_counter()
            if self.verbose:
                print(f"      -> getpyobsAOP took {t2-t1:.4f}s")

        else:
            tau = self.aer.tau[:,ich:ich+1,:].values
            ssa = self.aer.ssa[:,ich:ich+1,:].values
            pmatrix = self.aer.pmatrix[:,ich:ich+1,:,:,:].values
            g   = self.aer.g[:,ich:ich+1,:].values

        # get vertical levels dims are [nlev+1,nobs]
        pe   = self.edges.pe.values.astype('float64')
        ze   = self.edges.ze.values.astype('float64')
        te   = self.edges.te.values.astype('float64')

        # get angles. dims are [nobs]
        vza = self.geoms.VZA.values.astype('float64')
        sza = self.geoms.SZA.values.astype('float64')
        raa = self.geoms.RAA.values.astype('float64')

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

        # surface parameters
        param = self.surface.RTLSparam[:,ich:ich+1,:].values.astype('float64')
        kernel_wt = self.surface.kernel_wt[:,ich:ich+1,:].values.astype('float64')

        args = (rot, depol_ratio, alpha, tau, ssa, g, pmatrix, tauI, ssaI, gI, pmatrixI, tauL, ssaL, gL, pmatrixL, pe, te, ze, param, kernel_wt, vza, sza, raa, flux_factor)

        return args
        
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
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

    parser.add_argument("-n","--nproc", default=nproc, type=int,
                        help="Number of processors to use (default=%i)"%nproc)

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("--debug",action="store_true",
                        help="Debug mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")

    args = parser.parse_args()

    config = yaml.safe_load(open(args.inputs_yaml))

    mtFile = config.get('MTFILE','m2_aop.yaml')
    albedoType = config.get('ALBEDOTYPE','AMES_BRDF')
    plane_parallel = not config.get('DO_SPHERICITY',False)
    do_fullrad = config.get('DO_FULLRAD',True)
    NSTOKES = config.get('NSTOKES',3)
    write_inputs = config.get('WRITE_INPUTS',True)
    do_rtcalc = config.get('DO_RTCALC',True)
    do_inputscalc = config.get('DO_INPUTSCALC',True)

    args.paths_yaml = config['paths_yaml']
    args.inst_yaml  = config['inst_yaml']
    args.orbit_yaml = config['orbit_yaml']

    # Parse yamls
    # -----------------
    cf             = yaml.safe_load(open(args.inst_yaml))
    instname       = cf['instname']

    cf             = yaml.safe_load(open(args.orbit_yaml))
    orbitname      = cf['orbitname']
    ORBITNAME      = orbitname.upper()

    cf             = yaml.safe_load(open(args.paths_yaml))
    inTemplate     = cf['inDir']     + '/' + cf['inFile']
    outTemplate    = cf['outDir']    + '/' + cf['outFile']
    DT_mins        = cf.get('DT_MINS',1)

    brdfTemplate = cf['brdfDir'] + '/' + cf['brdfFile']


    # Loop through dates, running VLIDORT
    # ------------------------------------
    date      = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)
    Dt        = timedelta(minutes=DT_mins)

    while date < enddate:
        nymd  = str(date.date()).replace('-','')
        year  = str(date.year)
        month = str(date.month).zfill(2)
        day   = str(date.day).zfill(2)
        hour  = str(date.hour).zfill(2)
        minute = str(date.minute).zfill(2)

        replacements = {
            '%year': year, '%month': month, '%day': day, '%nymd': nymd,
            '%hour': hour, '%minute': minute, '%orbitname': orbitname, 
            '%ORBITNAME': ORBITNAME, '%instname': instname
        }


        inFile = inTemplate
        outFile = outTemplate
        argsFile = outTemplate.replace('vlidort', 'vlidort_args')
        brdfFile = brdfTemplate
        for k, v in replacements.items():
            inFile = inFile.replace(k, v)
            outFile = outFile.replace(k, v)
            argsFile = argsFile.replace(k, v)
            brdfFile = brdfFile.replace(k, v)

        argsFile = argsFile.replace('.nc', '.zarr')

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
        print('>>>do_fullrad',do_fullrad)
        print('>>>nproc:     ',args.nproc)
        print('++++End of arguments+++')
        print('') 
        vlidort = SBG_VLIDORT(inFile,outFile,argsFile,mtFile,
                            albedoType,
                            instname,
                            NSTOKES=NSTOKES,
                            brdfFile=brdfFile,
                            verbose=args.verbose,
                            debug=args.debug,
                            plane_parallel=plane_parallel,
                            do_fullrad=do_fullrad,
                            do_inputscalc=do_inputscalc,
                            nproc=args.nproc)

        # Run VLIDORT
        # ------------
        if vlidort.nobs == 0:
            print('No valid pixels found. Nothing to do.')
            sys.exit(0)

        if args.dryrun:
            print('Dry run enabled. Skipping execution.')
            sys.exit(0)


        # Initiate Output Arrays
        vl_initOutputs(vlidort)
        ts_initOutputs(vlidort)

        # Get the index of good obs
        iGood = np.where(vlidort.iGood)[0]

        if do_rtcalc:
            p = Pool(vlidort.nproc)


        # loop through nobs in batches
        for sob in range(0,vlidort.nobs,vlidort.nbatch):
            mem = psutil.virtual_memory()
            print(f'sob: {sob}, Available: {mem.available/1e9:.1f} GB, Used: {mem.percent}%')

            print(f'sob: {sob}, nobs: {vlidort.nobs}')

            eob = min([vlidort.nobs, sob + vlidort.nbatch])
            iobs = iGood[sob:eob]
           
            # Subset inputs for batch
            vlidort.aer = vlidort.AER.isel(nobs=slice(sob,eob)).load()
            vlidort.edges = vlidort.EDGES.isel(nobs=slice(sob,eob)).load()
            vlidort.rayleigh = vlidort.RAYLEIGH.isel(nobs=slice(sob,eob)).load()
            if vlidort.do_inputscalc:
                vlidort.geoms = vlidort.GEOMS.isel(nobs=iobs).load()
                vlidort.surface = vlidort.SURFACE.isel(nobs=iobs).load()
            else:
                vlidort.geoms = vlidort.GEOMS.isel(nobs=slice(sob,eob)).load()
                vlidort.surface = vlidort.SURFACE.isel(nobs=slice(sob,eob)).load()

            # Loop through channels
            for ich,channel in enumerate(vlidort.channels):
                print(f'ich: {ich}  channel: {channel}')

                # Get optical property inputs
                t_start = time.perf_counter()
                batch_args  = vlidort.getargs(ich)
                t_end = time.perf_counter()
                if vlidort.verbose:
                    print(f'   -> getargs took: {t_end - t_start:.4f} seconds')

                if write_inputs: 
                    t_start = time.perf_counter()
                    vlidort.writeArgs(ich, sob, eob, iobs, batch_args, 'RTLS')
                    t_end = time.perf_counter()
                    if vlidort.verbose:
                        print(f'   -> writeargs took: {t_end - t_start:.4f} seconds')

                if do_rtcalc:
                    print('   - run TWOSTREAM')
                    t_start = time.perf_counter()
                    vlidort.runTWOSTREAM(p,ich,sob,batch_args)
                    t_end = time.perf_counter()
                    if vlidort.verbose:
                        print(f'   -> TWOSTREAM took: {t_end - t_start:.4f} seconds')


                    print('   - run VLIDORT')
                    t_start = time.perf_counter()
                    vlidort.runVLIDORT(p,ich,sob,batch_args)
                    t_end = time.perf_counter()
                    if vlidort.verbose:
                        print(f'   -> VLIDORT took: {t_end - t_start:.4f} seconds')

        # Write outputs
        vlidort.writeNC()
            
#        # rewrite args onto the correct grid
#        vlidort.expand_dimensions() 

        if do_rtcalc:
            pool.close()
            pool.join()


        date += Dt
