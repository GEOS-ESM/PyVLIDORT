#!/usr/bin/env python3
"""
Concatenate chunked VLIDORT RT output files along ch dimension.
"""
import xarray as xr
import glob
import sys

def concatenate_chunks(outFile_template):
    """
    outFile_template: e.g., 'path/to/ssd650-sbg-g5nr.lc.vlidort.20060116_1735z.nc4'
    Finds all chunk files like *_ch000-053.nc4, *_ch053-106.nc4, etc.
    """
    base = outFile_template.replace('.nc4', '')
    chunk_files = sorted(glob.glob(f'{base}_ch*.nc4'))

    if not chunk_files:
        print(f"No chunk files found matching: {base}_ch*.nc4")
        sys.exit(1)

    print(f"Found {len(chunk_files)} chunk files:")
    for f in chunk_files:
        print(f"  {f}")

    # Open and concatenate along ch dimension
    ds = xr.open_mfdataset(chunk_files, concat_dim='ch', combine='nested')

    # Write combined output
    print(f"Writing concatenated file: {outFile_template}")
    ds.to_netcdf(outFile_template)
    print("Done.")

if __name__ == "__main__":
    concatenate_chunks(sys.argv[1])
