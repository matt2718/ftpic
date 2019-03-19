#!/bin/bash

export OMP_NUM_THREADS=16

for np in $(seq 1000 1000 10000); do
    for ng in $(seq 16 16 160); do
	echo "Np=$np, Ng=$ng"
	./oldpic2d -q -pt -c 2stream -np $np -ng $ng -t 0.001,0.1 | grep 'ms per'
	./ftpic2d  -q -pt -c 2stream -np $np -ng $ng -t 0.001,0.1 | grep 'ms per'
    done
done
