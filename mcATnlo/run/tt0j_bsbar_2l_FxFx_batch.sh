#!/bin/bash

date
hostname

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cms/scratch/iwatson/tsW/mcATnlo/CMSSW_8_0_30/src
eval `scramv1 runtime -sh`
cd -

echo "-- Setup CMS environment"

TDIR=`mktemp -d`

mkdir -p $TDIR
cp /cms/scratch/iwatson/tsW/mcATnlo/run/*py $TDIR
cd $TDIR

echo `pwd` @ `hostname` [ `date` ] > test

xrdcp test root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/tt0j_bsbar_2l_FxFx/LOG/$1.start

echo "-- Finished setup"
pwd
echo "-- Running Step 1"
cmsRun ./tt0j_bsbar_2l_FxFx_LHEGS.py > step1.log 2>&1
echo "-- Running Step 2"
cmsRun ./step2_DIGI_L1_DIGI2RAW_HLT_2016.py > step2.log 2>&1
echo "-- Running Step 3"
cmsRun ./step3_RAW2DIGI_L1Reco_RECO_EI_PAT_2016.py > step3.log 2>&1

rm step1.root
rm step2.root

xrdcp step1.log root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/tt0j_bsbar_2l_FxFx/LOG/$1_step1.log
xrdcp step2.log root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/tt0j_bsbar_2l_FxFx/LOG/$1_step2.log
xrdcp step3.log root:///xrootd/store/user/iawatson/tt0j_bsbar_2l_FxFx/LOG/$1_step3.log
xrdcp step3.root root:///xrootd/store/user/iawatson/tt0j_bsbar_2l_FxFx/GEN/$1.root
xrdcp step3_inMINIAODSIM.root root:///xrootd/store/user/iawatson/tt0j_bsbar_2l_FxFx/MINIAODSIM/$1.root

echo "-- Done."
date
