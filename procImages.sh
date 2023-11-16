#!/bin/bash

# Stop execution and exit on error
#set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

mkdir -p fits
mkdir -p raw

for i in `ls *.dng`; do
	echo =============== $i ===============
	python3 $SCRIPT_DIR/raw2fits.py $i
	echo
	mv $i `echo $i| sed 's#dng$#jpg#'` raw
done

mv *.fits fits

echo =====================================================
echo DONE!
echo All the processed raw files \(and jpg\) were moved to
echo a new folder called \"raw\"
echo =====================================================
echo 