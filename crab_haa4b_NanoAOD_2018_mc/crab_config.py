
'''
crab submit -c crab/crabConfigMC.py

crab status -d <dir name>
'''

import CRABClient
from CRABClient.UserUtilities import config
import os

config = config()

config.General.workArea = 'crab/r1'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Nano_Hto4bPlus_2018MC_cfg.py' 
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.inputFiles = ['crab_script.sh', 'Hto4b_postproc.py']
config.JobType.outputFiles = ['PNet_v1.root', 'PNet_v1_Skim.root']
config.JobType.maxMemoryMB = 10000  ## 2500*4
config.JobType.maxJobRuntimeMin = 1315  ## 21.9 hours
config.JobType.numCores = 4
# config.JobType.priority = 100  ## High priority for small test 

#config.Data.inputDataset = '/SUSY_ZH_ZToAll_HToAATo4B_Pt150_M-20_TuneCP5_13TeV_madgraph_pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'
config.Data.inputDataset = 'DUMMY'

config.General.requestName = 'DUMMY'
config.Data.outputDatasetTag = 'r1'
config.Data.publication = False

config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased' # 'Automatic' 'LumiBased' 'FileBased' 'EventAwareLumiBased'
config.Data.unitsPerJob = 100000  ## For most samples, should be < 5 hours
#config.Data.unitsPerJob = 10000  ## Perform a small test
#config.Data.totalUnits = 100000  ## Perform a small test
config.Site.storageSite = 'T2_CH_CERN'  ## Choose your site
config.Data.outLFNDirBase = '/store/group/phys_susy/HToaaTo4b/NanoAOD/2018/MC/PNet_v2_2024_11_22/'
