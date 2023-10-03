#!/usr/bin/python3
import os.path  
import sys
import numpy as np
from astropy.io import fits


def isfloat(str):
    try: 
        float(str)
    except ValueError: 
        return False
    return True


# get argument list using sys module
sys.argv

# Get the total number of args passed
total = len(sys.argv)

if total == 5:
	file = str(sys.argv[1])
	if os.path.isfile(file) == False:
		print( "\nThe fits file: ", file)
		print( "does not exist!\n")
		sys.exit(1) 
else:
	print( "\nAdds a field to the header of the FITS file if the field doesn\'t exist.")
	print( "It updates the value of the field if it already exist.")
	print( "\nUsage:")
	print("  %s <input file> <hdu id> <card name> <card value>\n" % (sys.argv[0]))
	sys.exit(0)


file     = str(sys.argv[1])
hduid    = int(sys.argv[2])
cardName = str(sys.argv[3])
value    = str(sys.argv[4])
hdulist = fits.open(file, mode='update')
if hduid>=len(hdulist):
	print( "\nError: the requested HDU header does not exist\n\n")
	sys.exit(1)
if isfloat(value):
	hdulist[hduid].header[cardName] = float(value)
else:
	hdulist[hduid].header[cardName] = value
hdulist.flush()
hdulist.close()
print( "Done!")
sys.exit(0)
