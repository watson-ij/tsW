#!/bin/sh

CMS=/cms/scratch/iwatson/tsW/cms/CMSSW_9_4_6_patch1/src

source /cvmfs/cms.cern.ch/cmsset_default.sh

cd ${CMS}
eval `scramv1 runtime -sh`
cd -

# slCalibAnalysis '/xrootd/store/group/nanoAOD/run2_2016v4/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180430_152649/0000/*root' TTJets.root
# slCalibAnalysis '/xrootd/store/group/nanoAOD/run2_2016v4/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180430_152649/0000/*root' TTJets.root

slCalibAnalysis '/xrootd/store/group/nanoAOD/run2_2016v3/SingleMuon/Run2016C-18Apr2017-v1/180124_072705/0000/*.root' SingleMuon_RunC_2016v3.root
slCalibAnalysis '/xrootd/store/group/nanoAOD/run2_2016v4/SingleMuon/Run2016C-07Aug17-v1/180504_150105/0000/*.root' SingleMuon_RunC_2016v4.root
