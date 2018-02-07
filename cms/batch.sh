#!/bin/bash

date

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cms/scratch/iwatson/tsW/cms/CMSSW_9_4_0/src
eval `scramv1 runtime -sh`
cd ../..

echo "-- Setup CMS environment"

TDIR=`mktemp --tmpdir=result`
cp *py $TDIR
cd $TDIR

echo "-- Running Step 1"
cmsRun ./TTbar_13TeV_TuneCUETP8M1_cfi_GEN_SIM.py > step1.log 2>&1
echo "-- Running Step 2"
cmsRun ./step2_DIGI_L1_DIGI2RAW_HLT.py > step2.log 2>&1
echo "-- Running Step 3"
cmsRun ./step3_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM.py > step3.log 2>&1

echo "-- Done."
date
