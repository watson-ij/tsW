#!/bin/sh

PROCESS=tt01j_bsbar_2l_FxFx_AOD

mkdir -p /xrootd/store/user/iawatson/$PROCESS/LOG
mkdir -p /xrootd/store/user/iawatson/$PROCESS/GEN
mkdir -p /xrootd/store/user/iawatson/$PROCESS/MINIAODSIM
mkdir -p /xrootd/store/user/iawatson/$PROCESS/HADAOD
mkdir -p /xrootd/store/user/iawatson/$PROCESS/NANOAOD
mkdir -p /xrootd/store/user/iawatson/$PROCESS/AOD

for i in {1001..1010}; do
    condor_submit ${PROCESS}_batch.jds -append "arguments = $i"
done
