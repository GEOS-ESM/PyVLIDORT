# input and output attributes
from .constants import MISSING

ARGS_ATTS : {
    "ch" : {
                        "standard_name": "wavelength",
                        "long_name": "wavelength",
                        "units"":" "nm",
                        "missing_value": MISSING,
           },
    "ang" : {
                        "standard_name": "angle",
                        "long_name": "scattering angle",
                        "units": "degrees",
                        "missing_value": MISSING,
            },
    "rot" : {
                        "standard_name": "ROT",
                        "long_name": "rayleigh optical thickenss",
                        "units": "None",
                        "missing_value": MISSING,
            },
    "depol" : {
                        "standard_name": "depolarization ratio",
                        "long_name": "Rayleigh depolarization ratio",
                        "units": "None",
                        "missing_value": MISSING,
              },
    "tau" : {
                        "standard_name": "TAU",
                        "long_name": "aerosol optical thickenss",
                        "units": "None",
                        "missing_value": MISSING,
            },
    "ssa" : {
                        "standard_name": "SSA",
                        "long_name": "single scattering albedo",
                        "units": "None",
                        "missing_value": MISSING,
            },
    "g" : {
                        "standard_name": "G",
                        "long_name": "assymetry parameter",
                        "units": "None",
                        "missing_value": MISSING,
          },
    "pmatrix" : {
                        "standard_name": "PMATRIX",
                        "long_name": "aerosol scattering matrix (p11,22,33,44,12,34)",
                        "units": "None",
                        "missing_value": MISSING,
                },
    "pe" : {
                        "standard_name": "PE",
                        "long_name": "pressure at layer edges",
                        "units": "Pa",
                        "missing_value": MISSING,
            },
    "ze" : {
                        "standard_name": "ZE",
                        "long_name": "height above sea level at layer edges",
                        "units": "m",
                        "missing_value": MISSING,
           },
    "te" : {
                        "standard_name": "TE",
                        "long_name": "temperature at layer edges",
                        "units": "K",
                        "missing_value": MISSING,
            },
    "albedo" : {
                        "standard_name": "ALBEDO",
                        "long_name": "lambertian albedo",
                        "units": "None",
                        "missing_value": MISSING,
                },
    "sza" : {
                        "standard_name": "SZA",
                        "long_name": "solar zenith angle",
                        "units": "degrees",
                        "missing_value": MISSING,
            },
    "vza" : {
                        "standard_name": "VZA",
                        "long_name": "sensor zenith angle",
                        "units": "degrees",
                        "missing_value": MISSING,
             },
    "raa" : {
                        "standard_name": "RAA",
                        "long_name": "relative azimuth angle",
                        "units": "degrees",
                        "missing_value": MISSING,
            },
    "flux" : {
                        "standard_name": "FLUX",
                        "long_name": "solar flux (F0)",
                        "units": "Done",
                        "missing_value": MISSING,
            },
    }
