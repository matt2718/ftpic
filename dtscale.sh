#!/bin/bash

mkdir -p dtscale

for i in {0..10}; do
    instr="0.0001 * 2 ^ $i"
    dt=$(echo "$instr" | bc -l)

    echo "Old: 0${dt} ..."
    oprefix="dtscale/o-0${dt}"
#    ./oldpic -q -t ${dt},10 -m "${oprefix}-mode.csv" > "${oprefix}.csv"
    ./oldpic -q -t ${dt},10 -c 2stream > "${oprefix}.csv"
    
    echo " FT: 0${dt} ..."
    fprefix="dtscale/f-0${dt}"
    ./ftpic -q -t ${dt},10 -c 2stream > "${fprefix}.csv"
done
