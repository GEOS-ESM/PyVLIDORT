#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
 
#SBATCH -J lamb_vlidort
#SBATCH --nodes=1
#SBATCH --constraint=mil
#SBATCH --time=1:00:00
#SBATCH -A @GROUPID
#SBATCH -o output_lamb_vlidort-%j.log
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --qos=debug
#######################################################################
#  Run vlidort code for Lambertian Test
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
    ln -s /home/pcastell/opendap/dasilva_fvinput/ExtData/chemistry/AerosolOptics/v1.0.0/  ExtData
endif

python3 -u ./lamb_vlidort_pyexample.py --nproc=1 2006-01-16T17:35 2006-01-16T17:36  lamb_vlidort.yaml 0.0 >& lamb_vlidort-${SLURM_JOB_ID}.log
