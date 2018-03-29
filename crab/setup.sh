#!/bin/sh

cd CMSSW_8_0_30/src
eval `scramv1 runtime -sh`
cd -
source /cvmfs/cms.cern.ch/crab3/crab.sh
# voms-proxy-init --voms cms
