h1. Description
This program extracts hits and tracks, computes their relevant parameters and saves them in a TTree in ROOT file.

h2. To compile:
1. install cfitsio libraries (In Ubuntu: apt-get install libcfitsio-bin libcfitsio-dev) and [CERN ROOT](https://root.cern/)
2. make

h2. Usage examples:
 * ./extract.exe <input file> <mask file> -o <output filename>
 * ./extract.exe `ls pepe*.fits` -o <output filename>

h2. Options:
 * **-q** for quiet (no screen output)
 * **-b** for computing and subtracting flat baseline from the image 
 * **-w** for computing and subtracting running window baseline from the image 
 * **-l** for computing and subtracting baseline from the median of each row in the image 
 * **-n** for computing noise sigma from the image 
 * **-s** <HDU number> for processing a single HDU