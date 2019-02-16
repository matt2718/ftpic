#!/bin/bash

if [[ -n "$2" ]]; then
    outdir="$2"
else
    outdir="dtscale2d"
fi

# create run directory
if [[ "$1" == "-r" ]]; then
    mkdir -p $outdir
fi

# logging header
if [[ "$1" == "-a" ]]; then
    echo "dt,oldpic,ftpic";
fi

params="-q -c 2stream -np 20000 -ng 64"

# print parameter info to stderr
>&2 echo $params
./oldpic2d $params -t 0.001,0.01 -p /dev/stderr >/dev/null

#for i in {0..10}; do
for i in {0..10}; do
    dt=$(echo "0.00025 * 2 ^ $i" | bc -l)
    oname="$outdir/o-0${dt}.csv"
    fname="$outdir/f-0${dt}.csv"
   
    if [[ "$1" == "-r" ]]; then
        # run
        echo "Old: 0${dt} ..."
        ./oldpic2d $params -t ${dt},5 > "${oname}"

        echo " FT: 0${dt} ..."
        ./ftpic2d $params -t ${dt},5 > "${fname}"

    elif [[ "$1" == "-a" ]]; then
        # analyze
        echo $dt,$(./getmax.py $oname),$(./getmax.py $fname)
    fi
done

exit 0

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
