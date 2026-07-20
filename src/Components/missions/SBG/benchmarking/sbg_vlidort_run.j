#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
 
#SBATCH -J sbg_vlidort
#SBATCH --nodes=1
#SBATCH --constraint=mil
#SBATCH --time=1:00:00
#SBATCH -A @GROUPID
#SBATCH -o output_sbg_vlidort-%j.log
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --qos=debug
#######################################################################
#  Run vlidort code for SBG
#######################################################################
#           Architecture Specific Environment Variables
#######################################################################

setenv SRC_DIR @SRCDIR
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules


#######################################################################
#          Run Sampler 
#######################################################################
if (! -d ExtData) then
    ln -s /discover/nobackup/pcastell/workspace/GEOSmie_refactor/dustupdate/AerosolOptics/v2.0.0/x  ExtData
endif

python3 -u ./sbg_vlidort_pyexample.py --nproc=1 2006-01-16T17:35 2006-01-16T17:36  sbg_vlidort.yaml >& sbg_vlidort-${SLURM_JOB_ID}.log
