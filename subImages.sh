#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


for i in `ls d*.fits`; do
	#statements
	echo ========== $i ========== 	
	python3 $SCRIPT_DIR/subFits.py $i median.fits
	echo
done