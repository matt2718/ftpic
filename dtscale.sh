#!/bin/bash

mkdir -p dtscale

for i in {0..9}; do
    instr="0.0001 * 2 ^ $i"
    dt=$(echo "$instr" | bc -l)

    echo "Old: ${dt} ..."
    oprefix="dtscale/o-0.${dt}"
    ./oldpic -q -t ${dt},10 -m "${oprefix}-mode.csv" > "${oprefix}.csv"

    echo " FT: ${dt} ..."
    fprefix="dtscale/f-0.${dt}"
    ./ftpic -q -t ${dt},10 -m "${fprefix}-mode.csv" > "${fprefix}.csv"
done
