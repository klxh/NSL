#!/bin/bash

lc=$(wc input.dat | awk '{print $1}')
line_count=$(echo "$lc - 1" | bc)

./clean.sh

for i in {0..10};
do
    temp=$(echo "$i * 0.15 + 0.5" | bc)
    echo $temp > tmp.in
    tail -$line_count input.dat >> tmp.in
    mv tmp.in input.dat
    ./Monte_Carlo_ISING_1D.exe
done
