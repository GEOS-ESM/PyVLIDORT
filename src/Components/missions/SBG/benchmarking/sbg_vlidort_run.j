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
    ln -s /home/pcastell/opendap/dasilva_fvinput/ExtData/chemistry/AerosolOptics/v2.x.x/  ExtData
endif

python3 -u ./sbg_vlidort_pyexample.py 2006-01-16T17:35 2006-01-16T17:36  sbg_vlidort.yaml >& aaq_sampler-${SLURM_JOB_ID}.log
