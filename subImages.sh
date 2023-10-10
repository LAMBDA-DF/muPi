#!/bin/bash

# Stop execution and exit on error
set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

baseDir=$1
baseDir="${baseDir:=.}"

masterBias=$2
masterBias="${masterBias:=median.fits}"

for i in `ls $baseDir/*.fits`; do
	outputName=mbs_`echo $i | sed 's#.*/##'`
	echo ========== $i ========== 	
	# echo $SCRIPT_DIR/subFits.py $i $masterBias $outputName
	$SCRIPT_DIR/subFits.py $i $masterBias $outputName
	echo
done
