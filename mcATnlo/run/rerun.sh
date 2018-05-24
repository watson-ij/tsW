#!/bin/sh

hostname

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home/iwatson/tsW/cms/CMSSW_9_4_4/src
eval `scramv1 runtime -sh`
eval `scramv1 runtime -sh`
cd -

cd `mktemp -d`

echo `pwd`
cmsRun /home/iwatson/tsW/mcATnlo/run/run2_2016MC_NANO.py $1/MINIAODSIM/$2.root > $2.nano.log 2>&1
cp nano*root /home/iwatson/tsW/mcATnlo/data/$1/NANOAOD/$2.root
cp *log /home/iwatson/tsW/mcATnlo/data/$1/LOG/$2.nano.log
rm nano*root

echo `pwd`
cmsRun /home/iwatson/tsW/mcATnlo/run/run2_2016MC_HADAOD.py $1/AODSIM/$2.root > $2.nano.log 2>&1
cp nanoAOD_AOD.root /home/iwatson/tsW/mcATnlo/data/$1/HADAOD/$2.root
cp *log /home/iwatson/tsW/mcATnlo/data/$1/LOG/$2.had.log
rm nanoAOD_AOD.root

echo `pwd`
cmsRun /home/iwatson/tsW/mcATnlo/run/hadAOD.py $1/GEN/$2.root > $2.nano.log 2>&1
cp hadAOD.root /home/iwatson/tsW/mcATnlo/data/$1/HADTRUTHAOD/$2.root
cp *log /home/iwatson/tsW/mcATnlo/data/$1/LOG/$2.had.log
rm hadAOD.root
