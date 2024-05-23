#!/bin/bash

for g in 1
do
    for t in 60
    do

        for (( i = 0; i < 30; i++ ))
        do
            echo ################## $i ###################
            libcamera-still -n -r -o ./img/muons/Camera0_raw_g$g"_"$t"s_"$i.jpg --camera 0 --shutter $t"000000" --gain $g --awbgains 1,1 --immediate > /dev/null 2>&1 &
            pid1=$!
            libcamera-still -n -r -o ./img/muons/Camera1_raw_g$g"_"$t"s_"$i.jpg --camera 1 --shutter $t"000000" --gain $g --awbgains 1,1 --immediate > /dev/null 2>&1 &
            pid2=$!

            echo "camera 0 (1) started with PID:$pid1 ($pid2)"
            # wait for all pids
            for pid in $pid1 $pid2
            do
               wait $pid
	       echo "finished $pid"
            done
        done
        done
done

