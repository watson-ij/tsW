#!/bin/sh

BASEP="/xrootd/store/user/iawatson"
BASE="root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson"
XRDCP=1

#BASEP=/home/iwatson/tsW/mcATnlo/data
#BASE=/home/iwatson/tsW/mcATnlo/data
#XRDCP=0

PROCESS="tt012j_bbars_2l_FxFx"

mkdir -p $BASEP/$PROCESS/LOG
mkdir -p $BASEP/$PROCESS/GEN
mkdir -p $BASEP/$PROCESS/MINIAODSIM
mkdir -p $BASEP/$PROCESS/HADAOD
mkdir -p $BASEP/$PROCESS/HADTRUTHAOD
mkdir -p $BASEP/$PROCESS/NANOAOD
mkdir -p $BASEP/$PROCESS/AODSIM

if [ "$HOSTNAME" = gate.sscc.uos.ac.kr ]; then
    chmod o+rwx $BASEP/$PROCESS/
    chmod o+rwx $BASEP/$PROCESS/*
fi

for i in {1..25}; do
    condor_submit batch.jds -append "arguments = $BASE $PROCESS $i $XRDCP"
done
