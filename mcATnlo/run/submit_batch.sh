#!/bin/sh

BASEP="/xrootd/store/user/iawatson"
BASEX="/xrd/store/user/iawatson"
XRDFS="root://cms-xrdr.sdfarm.kr:1094/"
BASE=$XRDFS"/"$BASEX
XRDCP=1

#BASEP=/home/iwatson/tsW/mcATnlo/data
#BASE=/home/iwatson/tsW/mcATnlo/data
#XRDCP=0

PROCESS="tt012j_bbar_1lm_FxFx"

xrdfs $XRDFS mkdir -p $BASEX/$PROCESS/LOG
xrdfs $XRDFS mkdir -p $BASEX/$PROCESS/GEN
xrdfs $XRDFS mkdir -p $BASEX/$PROCESS/MINIAODSIM
xrdfs $XRDFS mkdir -p $BASEX/$PROCESS/HADAOD
xrdfs $XRDFS mkdir -p $BASEX/$PROCESS/HADTRUTHAOD
xrdfs $XRDFS mkdir -p $BASEX/$PROCESS/NANOAOD
xrdfs $XRDFS mkdir -p $BASEX/$PROCESS/AODSIM

if [ "$HOSTNAME" = gate.sscc.uos.ac.kr ]; then
    chmod o+rwx $BASEP/$PROCESS/
    chmod o+rwx $BASEP/$PROCESS/*
fi

for i in {51..150}; do
    condor_submit batch.jds -append "arguments = $BASE $PROCESS $i $XRDCP"
done
