#!/bin/bash

echo threads,oldpic,ftpic

for th in {1..16}; do
    export OMP_NUM_THREADS=$th

    otime=$(/usr/bin/time -f '%e' ./oldpic -q -c wave 2>&1 >/dev/null)
    ftime=$(/usr/bin/time -f '%e' ./ftpic -q -c wave 2>&1 >/dev/null)

    echo $th,$otime,$ftime
done
