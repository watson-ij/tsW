#!/bin/sh

PROCESS=tt012j_bbars_2l_FxFx_AOD

mkdir -p /home/iwatson/tsW/mcATnlo/data/$PROCESS/LOG
mkdir -p /home/iwatson/tsW/mcATnlo/data/$PROCESS/GEN
mkdir -p /home/iwatson/tsW/mcATnlo/data/$PROCESS/MINIAODSIM
mkdir -p /home/iwatson/tsW/mcATnlo/data/$PROCESS/HADAOD
mkdir -p /home/iwatson/tsW/mcATnlo/data/$PROCESS/NANOAOD
mkdir -p /home/iwatson/tsW/mcATnlo/data/$PROCESS/AODSIM

for i in {1..25}; do
    condor_submit ${PROCESS}_batch.jds -append "arguments = $i"
done
