#!/usr/bin/python3

import sys
import rawpy
from astropy.io import fits
import numpy as np
import re
import os.path

def printHelp(progName):
    print("")
    print("This program reads a dng (RAW) file and converts it to ROOT format.")
    print("It does a 1 to 1 conversion, no Demosaic nor any postprocessing.")
    print("All the pixels are treated the same regardless of their color.\n")
    print("Usage:")
    print("    ",progName, " <dng file>\n\n")


def isfloat(str):
    try: 
        float(str)
    except ValueError: 
        return False
    return True


def extractRunID(s):
    # Regular expression pattern to match the desired format
    pattern = r'.*_(\d+).dng$'

    match = re.search(pattern, s)
    if match:
        # Extract the integer part
        number = int(match.group(1))
        return number
    return -1


def dng_to_fits(dng_filename, fits_filename):
    # Step 1: Read the DNG file
    raw = rawpy.imread(dng_filename)
    rawImg = raw.raw_image

    # Create a new FITS file
    hdu = fits.PrimaryHDU(rawImg)
    hdulist = fits.HDUList([hdu])

    # Read RunID from the file name and add it to the header of the fits image
    runID = extractRunID(dng_filename)
    if isfloat(runID):
        hdulist[0].header["RUNID"] = float(runID)
    else:
        hdulist[0].header["RUNID"] = runID
    hdulist.writeto(fits_filename, overwrite=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        printHelp(sys.argv[0])
        exit(1)
    inFileName = sys.argv[1]
    outFileName = inFileName[:-3] + "fits"

    if os.path.isfile(inFileName) == False:
        print( "\nThe file: ", inFileName)
        print( "does not exist!\n")
        sys.exit(1) 

    print("Reading: ", inFileName)
    print("Will output: ", outFileName)
    dng_to_fits(inFileName, outFileName)
    print("Done!")
