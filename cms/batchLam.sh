#!/bin/bash

date
hostname

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cms/scratch/iwatson/tsW/cms/CMSSW_9_4_6_patch1/src
eval `scramv1 runtime -sh`
cd -

echo "-- Setup CMS environment"

if [[ $1 == "2" ]]; then
    TDIR=./result/`date +%Y%m%d_%H%M%S`_`hostname`_ttsWbW_lam
fi
if [[ $1 == "1" ]]; then
    TDIR=./result/`date +%Y%m%d_%H%M%S`_`hostname`_tsW_lam
fi
if [[ $1 == "0" ]]; then
    TDIR=./result/`date +%Y%m%d_%H%M%S`_`hostname`_tbW_lam
fi
mkdir -p $TDIR
cp /cms/scratch/iwatson/tsW/cms/*py $TDIR
cd $TDIR

SRT_PYTHIA8DATA_SCRAMRTDEL=/cms/scratch/iwatson/tsW/pythia-xml
PYTHIA8DATA=/cms/scratch/iwatson/tsW/pythia-xml

echo "-- Finished setup"
pwd
echo "-- Running Step 1"
cmsRun ./TTbar_13TeV_TuneCUETP8M1_cfi_GEN_SIM_LamJPsiFilter.py $1 > step1.log 2>&1
echo "-- Running Step 2"
cmsRun ./step2_DIGI_L1_DIGI2RAW_HLT.py > step2.log 2>&1
echo "-- Running Step 3"
cmsRun ./step3a_RAW2DIGI_L1Reco_RECO_RECOSIM_EI_PAT_VALIDATION_DQM.py > step3.log 2>&1
echo "-- Running NanoAOD"
cmsRun nanoAOD_NANO.py > nano.log 2>&1
echo "-- Running HadAOD"
cmsRun hadAOD17.py > had.log 2>&1
rm step1.root
rm step2.root
# rm step3.root
echo "-- Done."
date
