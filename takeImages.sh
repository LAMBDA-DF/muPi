#!/bin/bash

for g in 4; do
for t in 50; do

  for (( i = 0; i < 500; i++ )); do
	echo ################## $i ################### 	
	libcamera-still -n -r -o 55fe_raw_g$g"_"$t"s_"$i.jpg --shutter $t"000000" --gain $g --awbgains 1,1 --immediate
	echo
	echo
  done

done
done
