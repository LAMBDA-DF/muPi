#!/bin/bash

# Stop execution and exit on error
set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


# 1. convert RAW data to FITS images
$SCRIPT_DIR/procImages.sh

# 2. compute median image to use as master-bias 
$SCRIPT_DIR/computeMedian/computeMedian.exe `ls fits/*.fits` -o median.fits

# 3. subtract master-bias from FITS images to equalize them
mkdir -p mbs
cd mbs
$SCRIPT_DIR/subImages.sh ../fits ../median.fits

# 4. create a mask file
$SCRIPT_DIR/createMask.py "*.fits"
mv mask.fits ../.

# 5. extract events and create a catalog
cp $SCRIPT_DIR/extractEvents/extractConfig.xml .
$SCRIPT_DIR/extractEvents/extract.exe `ls *.fits` -m ../mask.fits -n -o ../events.root
cd ..

echo
echo =====================================================
echo DONE!
echo The catalog file is events.root
echo =====================================================
echo 