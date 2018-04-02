from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'ttbar_lambdab_MC_generation_test_6'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'TTbarLepton_13TeV_TuneCUETP8M1_cfi_py_GEN_SIM_DIGIPREMIX_S2_DATAMIX_L1_DIGI2RAW_HLT.py'

config.Data.outputPrimaryDataset = 'TTLambdaB'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 20
NJOBS = 10  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'ttbar_lambdab_MC_generation_test_6'

config.Site.whitelist = ['T2_CH_CERN']
config.Site.storageSite = 'T3_KR_KISTI'
