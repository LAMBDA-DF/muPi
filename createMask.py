#!/usr/bin/python3

import sys
import glob
import numpy as np
from astropy.io import fits
from tqdm import tqdm


def printHelp(progName):
    print("")
    print("This program reads a set of FITS images and:")
    print("\t1. computes the pixel MAD from the first image (assuming they are all equivalent).")
    print("\t2. sets a threshold = 3*sigma, with sigma=1.4826*MAD.")
    print("\t3. counts the number of times each pixel is avobe a threshold.")
    print("\t4. the pixels that are avobe the threshold in more that 30\% of the images are marked with 1 (the others are 0).")
    print("\t5. creates an output mask image (in FITS format) that has 1 for the pixels that frequently avobe the threshold.")
    print("Usage:")
    print("    ",progName, " <\"name patern for the input files\">\n\n")
    print("Example:")
    print("    ",progName, " \"pepe*.fits\" \n\n")

def computeMedianAndMAD(data):
    # Flatten the data for computation
    flattened_data = data.flatten()

    # Compute the median
    median_val = np.median(flattened_data)

    # Compute the MAD
    mad_val = np.median(np.abs(flattened_data - median_val))

    return median_val, mad_val


def createMask(pattern):

    # Use glob to get fits_files matching the pattern
    fits_files = glob.glob(pattern)
    if not fits_files:
        raise ValueError("No files found matching the pattern.")

    print("Will read the following", len(fits_files),"files:")
    print(*fits_files, sep = ", ")


    # Read the first image to get the shape, median and MAD (assuming they are all equivalent)
    first_image  = fits.getdata(fits_files[0])
    image_shape  = first_image.shape
    image_median, image_MAD = computeMedianAndMAD(first_image)

    print("\n===================================================\n")
    print("Will use: median=",image_median," MAD=", image_MAD)
    print("Extracted from: ", fits_files[0])
    print("\n===================================================\n")

    # Threshold value
    threshold = 2*image_MAD*1.4826  # Threshold set to 3*sigma (sigma=MAD*1.4826)

    # Initialize a count array to zero
    count_array = np.zeros(image_shape, dtype=int)

    # Loop through the list of FITS files and count pixels above the threshold
    for file in tqdm(fits_files):
        image_data = fits.getdata(file)
        count_array += (image_data > threshold)

    # Calculate the fraction of images where each pixel is above the threshold
    fraction = count_array / len(fits_files)

    # Create the output image: pixels with fraction > 0.3 are 1, others are 0
    output_image = np.where(fraction > 0.3, 1, 0).astype(np.uint8)

    # Write the output image to a new FITS file
    hdu = fits.PrimaryHDU(output_image.astype(np.uint8))
    hdu.writeto('mask.fits', overwrite=True)

    #Print a summary
    non_zero_count = np.count_nonzero(output_image)
    print("Number of bad pixels:   ", non_zero_count)
    print("Fraction of bad pixels: ", non_zero_count/(image_shape[0]*image_shape[1]))
    print("\n===================================================\n")


if __name__ == "__main__":

    if len(sys.argv) != 2:
        printHelp(sys.argv[0])
        exit(1)
    inFileName = sys.argv[1]

    inFileNamePattern = sys.argv[1]

    createMask(inFileNamePattern)

