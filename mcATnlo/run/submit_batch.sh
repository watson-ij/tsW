#!/bin/sh

PROCESS=tt01j_bbar_1lp_FxFx_AOD

mkdir -p /home/iwatson/tsW/mcATnlo/data/$PROCESS/LOG
mkdir -p /home/iwatson/tsW/mcATnlo/data/$PROCESS/GEN
mkdir -p /home/iwatson/tsW/mcATnlo/data/$PROCESS/MINIAODSIM
mkdir -p /home/iwatson/tsW/mcATnlo/data/$PROCESS/HADAOD
mkdir -p /home/iwatson/tsW/mcATnlo/data/$PROCESS/HADTRUTHAOD
mkdir -p /home/iwatson/tsW/mcATnlo/data/$PROCESS/NANOAOD
mkdir -p /home/iwatson/tsW/mcATnlo/data/$PROCESS/AODSIM

chmod o+rwx /home/iwatson/tsW/mcATnlo/data/$PROCESS/
chmod o+rwx /home/iwatson/tsW/mcATnlo/data/$PROCESS/*

for i in {1001..1002}; do
    condor_submit ${PROCESS}_batch.jds -append "arguments = $i"
done
