#!/bin/bash

for name in {ene,hea,chi,mag}
do
    for i in {0..10}
    do 
        tail -1 output.$name.$i
    done | awk '{print $2, "\t", $3, "\t", $4}' > tmp_$name.dat
    paste temperatures tmp_$name.dat > gibbs_$name.dat
    rm tmp_$name.dat
done
