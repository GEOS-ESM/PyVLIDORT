#!/usr/bin/env python3
"""
Concatenate chunked VLIDORT RT output files along ch dimension.
Uses the same YAML configuration as the main sbg_vlidort.py script.
"""
import xarray as xr
import glob
import sys
import os
import argparse
import yaml
from datetime import timedelta
from dateutil.parser import parse as isoparser


def concatenate_chunks(outFile):
    """
    Find all chunk files and concatenate along ch dimension.
    """
    base = outFile.replace('.nc4', '')
    chunk_files = sorted(glob.glob(f'{base}_ch*.nc4'))

    if not chunk_files:
        print(f"No chunk files found matching: {base}_ch*.nc4")
        return False

    print(f"Found {len(chunk_files)} chunk files:")
    for f in chunk_files:
        print(f"  {f}")

    # Open and concatenate along ch dimension
    ds = xr.open_mfdataset(chunk_files, concat_dim='ch', combine='nested')

    # Write combined output
    print(f"Writing concatenated file: {outFile}")
    ds.to_netcdf(outFile)
    print(f"Done: {outFile}")

    return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Concatenate chunked VLIDORT output files")
    parser.add_argument("iso_t1", help="starting iso time")
    parser.add_argument("iso_t2", help="ending iso time")
    parser.add_argument("inputs_yaml", help="yaml file with input configuration file names")
    parser.add_argument("--cleanup", action="store_true",
                        help="Remove chunk files after successful concatenation")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose mode")

    args = parser.parse_args()

    # Parse same YAML files as main script
    config = yaml.safe_load(open(args.inputs_yaml))

    args.paths_yaml = config['paths_yaml']
    args.inst_yaml = config['inst_yaml']
    args.orbit_yaml = config['orbit_yaml']

    cf = yaml.safe_load(open(args.inst_yaml))
    instname = cf['instname']

    cf = yaml.safe_load(open(args.orbit_yaml))
    orbitname = cf['orbitname']
    ORBITNAME = orbitname.upper()

    cf = yaml.safe_load(open(args.paths_yaml))
    outTemplate = cf['outDir'] + '/' + cf['outFile']
    DT_mins = cf.get('DT_MINS', 1)

    # Loop through dates (same as main script)
    date = isoparser(args.iso_t1)
    enddate = isoparser(args.iso_t2)
    Dt = timedelta(minutes=DT_mins)

    while date < enddate:
        nymd = str(date.date()).replace('-', '')
        year = str(date.year)
        month = str(date.month).zfill(2)
        day = str(date.day).zfill(2)
        hour = str(date.hour).zfill(2)
        minute = str(date.minute).zfill(2)

        replacements = {
            '%year': year, '%month': month, '%day': day, '%nymd': nymd,
            '%hour': hour, '%minute': minute, '%orbitname': orbitname,
            '%ORBITNAME': ORBITNAME, '%instname': instname
        }

        outFile = outTemplate
        for k, v in replacements.items():
            outFile = outFile.replace(k, v)

        print(f"\n--- Processing: {outFile} ---")

        if os.path.exists(outFile):
            print(f"  Output already exists, skipping: {outFile}")
            date += Dt
            continue

        success = concatenate_chunks(outFile)

        # Optionally remove chunk files
        if success and args.cleanup:
            base = outFile.replace('.nc4', '')
            chunk_files = sorted(glob.glob(f'{base}_ch*.nc4'))
            for f in chunk_files:
                if args.verbose:
                    print(f"  Removing: {f}")
                os.remove(f)
            print(f"  Cleaned up {len(chunk_files)} chunk files.")

        date += Dt
