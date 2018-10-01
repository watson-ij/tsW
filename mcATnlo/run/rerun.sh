#!/bin/sh

hostname

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cms/scratch/iwatson/tsW/cms/CMSSW_9_4_6_patch1/src
eval `scramv1 runtime -sh`
eval `scramv1 runtime -sh`
cd -

cd `mktemp -d`

echo `pwd`
cmsRun /cms/scratch/iwatson/tsW/mcATnlo/run/run2_2016MC_NANO.py $1/MINIAODSIM/$2.root > $2.nano.log 2>&1
xrdcp nano*root root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/$1/NANOAOD/$2.root
xrdcp *log root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/$1/LOG/$2.nano.log
rm nano*root

# echo `pwd`
# cmsRun /home/iwatson/tsW/mcATnlo/run/run2_2016MC_HADAOD.py $1/AODSIM/$2.root > $2.nano.log 2>&1
# xrdcp nanoAOD_AOD.root root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/$1/HADAOD/$2.root
# xrdcp *log root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/$1/LOG/$2.had.log
# rm nanoAOD_AOD.root

# echo `pwd`
# # echo $1/GEN/$2.root
# cmsRun /cms/scratch/iwatson/tsW/mcATnlo/run/hadAOD.py $1/GEN/$2.root > $2.nano.log 2>&1
# xrdcp hadAOD.root root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/$1/HADTRUTHAOD/$2.root
# xrdcp *log root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/$1/LOG/$2.had.log
# rm hadAOD.root
