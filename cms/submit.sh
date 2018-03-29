#!/bin/sh

for i in {1..10}; do
    echo $i
    condor_submit tsWbatch_Lam.jds
    sleep 5
done
