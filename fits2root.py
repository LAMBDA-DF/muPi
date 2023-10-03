#!/usr/bin/python3

import sys
import rawpy
import uproot3
import numpy as np
from astropy.io import fits
import glob
import os.path


def printHelp(progName):
    print("")
    print("This program reads a FITS file and converts it to ROOT format.")
    print("The root file will contain a TTree named \"imgData\" that has the following")
    print("branches: \"x:y:pix\" where \"pix\" is the value of the (x,y) pixel.\n")
    print("Usage:")
    print("    ",progName, " <FITS file>\n\n")


def fits_to_root(filename, root_filename):

    # Step 0: Prepare data for the TTree
    x   = []
    y   = []
    pix = []
    
    with fits.open(filename) as hdul:
        # Assuming the image data is in the primary HDU
        data = hdul[0].data

        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                    x.append(j)
                    y.append(i)
                    pix.append(data[i,j])
        hdul.close()

    data = {
        "x":    np.array(x, dtype=np.int32),
        "y":    np.array(y, dtype=np.int32),
        "pix":  np.array(pix, dtype=np.float64),
    }
    
    # Step 3: Write data to the ROOT file
    with uproot3.recreate(root_filename) as f:
        f["imgData"] = uproot3.newtree({"x": "int32", "y": "int32", "pix": "float64"})
        f["imgData"].extend(data)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        printHelp(sys.argv[0])
        exit(1)
    inFileName = sys.argv[1]
    outFileName = inFileName + ".root"

    if os.path.isfile(inFileName) == False:
        print( "\nThe file: ", inFileName)
        print( "does not exist!\n")
        sys.exit(1) 

    print("Reading: ", inFileName)
    print("Will output: ", outFileName)
    fits_to_root(inFileName, outFileName)
    print("Done!")
