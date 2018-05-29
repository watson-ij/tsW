#!/bin/sh

BASEP="/xrootd/store/user/iawatson"
BASE="root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson"
XRDCP=1

#BASEP=/home/iwatson/tsW/mcATnlo/data
#BASE=/home/iwatson/tsW/mcATnlo/data
#XRDCP=0

PROCESS="tt01j_bsbar_2l_FxFx_AOD"

mkdir -p $BASEP/$PROCESS/LOG
mkdir -p $BASEP/$PROCESS/GEN
mkdir -p $BASEP/$PROCESS/MINIAODSIM
mkdir -p $BASEP/$PROCESS/HADAOD
mkdir -p $BASEP/$PROCESS/HADTRUTHAOD
mkdir -p $BASEP/$PROCESS/NANOAOD
mkdir -p $BASEP/$PROCESS/AODSIM

for i in {10001..10002}; do
    condor_submit batch.jds -append "arguments = $BASE $PROCESS $i $XRDCP"
done
