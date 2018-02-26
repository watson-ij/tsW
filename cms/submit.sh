#!/bin/sh

for i in {1..100}; do
    echo $i
    condor_submit tsWbatch.jds
    sleep 2
done
