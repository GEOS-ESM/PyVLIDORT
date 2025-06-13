#!/usr/bin/env python3

"""
    Parent class with utilities needed to 
    setup input data
    and run python wrappers for twostream
"""

from   py_vlidort import TWOSTREAM_ 

MISSING = -1.e+20

#---
def LAMBERTIAN_run(args):

    # Call TWO_STREAM wrapper function
    I, reflectance, rc = TWOSTREAM_.twostream_lambert_driver(*args)

    return I,reflectance

#---
def MODIS_BRDF_run(args):

    # Call TWO_STREAM wrapper function
    I, reflectance, rc = TWOSTREAM_.twostream_brdf_rtls(*args)

    return I,reflectance
