h1. Description

This program computes the median images from N input fit files.
It handles all the available HDUs. The HDUs in the output fit file
will be:
 * float (32bits) for:  int8, int16, int32 and float input images
 * double (64bits) for: double input images

h2. To compile:
1. install cfitsio libraries (In Ubuntu: apt-get install libcfitsio-bin libcfitsio-dev) 
2. make

h2. Usage examples:
 * ./computeMedian.exe <input file 1> .. <input file N> -o <output filename>
 * ./computeMedian.exe -i <input list file > -o <output filename>
 * ./computeMedian.exe <input file 1> .. <input file N> -i <input list file > -o <output filename>

h2. Options:
 * **-q** for quiet (no screen output)
 * **-s** <HDU number> for processing a single HDU
 * **-m** compute MAD instead of Median
 * **-x** compute sigma from MAD instead of Median. sigma = 1.4826\*MAD
 * **-a** to stack (add all images) instead of Median