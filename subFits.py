#!/usr/bin/python3
import sys
from astropy.io import fits

def printHelp(progName):
    print("")
    print("This program subtracts the data of the second FITS file from the first and saves")
    print("it to a new file. The name of the output file is \"mbs_\"<file1>:\n")
    print("Usage:")
    print("    ",progName, " <file1> <file2>\n")
    print("Example:")
    print("    ",progName, " \"pepe.fits\" \"median.fits\" \"mbs_pepe.fits\"\n\n")

def subtract_fits(file1, file2, output_file):
    """
    Subtract the data of the second FITS file from the first and save to a new file.
    
    Parameters:
        file1 (str): Path to the first FITS file.
        file2 (str): Path to the second FITS file.
        output_file (str): Path to the output FITS file.
    """
    
    with fits.open(file1) as hdul1, fits.open(file2) as hdul2:
        # Assuming the image data is in the primary HDU for both files
        data1 = hdul1[0].data
        data2 = hdul2[0].data
        
        if data1.shape != data2.shape:
            raise ValueError("The two FITS files have different data dimensions!")
        
        result_data = data1 - data2
        
        # Create a new HDU with the subtracted data and copy the header from the first file
        hdu = fits.PrimaryHDU(result_data, header=hdul1[0].header)
        
        # Save to output file
        hdu.writeto(output_file, overwrite=True)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        printHelp(sys.argv[0])
        exit(1)

    inFileName1 = sys.argv[1]
    inFileName2 = sys.argv[2]
    outFileName = sys.argv[3]
    print("Reading: ", inFileName1)
    print("Will output: ", outFileName)
    subtract_fits(inFileName1, inFileName2, outFileName)
