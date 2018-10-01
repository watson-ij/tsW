#!/bin/sh

BASEP="/xrootd/store/user/iawatson"
BASE="root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson"
XRDCP=1

#BASEP=/home/iwatson/tsW/mcATnlo/data
#BASE=/home/iwatson/tsW/mcATnlo/data
#XRDCP=0

PROCESS="tt012j_bbars_2l_FxFx"

mkdir -p $BASEP/${PROCESS}_herwigpp/LOG
mkdir -p $BASEP/${PROCESS}_herwigpp/GEN
mkdir -p $BASEP/${PROCESS}_herwigpp/MINIAODSIM
mkdir -p $BASEP/${PROCESS}_herwigpp/HADAOD
mkdir -p $BASEP/${PROCESS}_herwigpp/HADTRUTHAOD
mkdir -p $BASEP/${PROCESS}_herwigpp/NANOAOD
mkdir -p $BASEP/${PROCESS}_herwigpp/AODSIM

if [ "$HOSTNAME" = gate.sscc.uos.ac.kr ]; then
    chmod o+rwx $BASEP/${PROCESS}_herwigpp/
    chmod o+rwx $BASEP/${PROCESS}_herwigpp/*
fi

for i in {1..50}; do
    condor_submit batch_herwigpp.jds -append "arguments = $BASE $PROCESS $i $XRDCP"
done
