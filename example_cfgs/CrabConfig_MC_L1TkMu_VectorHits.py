import sys

from CRABClient.UserUtilities import config

config = config()

config.JobType.pluginName   = 'Analysis'
config.JobType.outputFiles  = ['muonNtuple_phase2_MC.root']

config.Data.unitsPerJob     = 100
config.Data.totalUnits      = 100000

config.Data.splitting       = 'EventAwareLumiBased'

config.Data.useParent       = True #!!!!
#config.Data.useParent       = False #!!!!

config.Site.storageSite     = 'T2_US_Purdue'
config.JobType.numCores     = 1
config.JobType.maxMemoryMB  = 2500
config.JobType.allowUndistributedCMSSW = True
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

tag = "muonHLT_phase2_DYToLL_PU140_L1TkMu_default_VHenabled"

config.General.workArea   = tag
config.Data.outLFNDirBase = '/store/user/dkondrat/' + tag

config.JobType.psetName    = 'HLT_test_L1TkMu_VectorHits.py'
config.General.requestName = tag
config.General.transferLogs = True

# no PU
#config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2HLTTDRWinter20RECOMiniAOD-NoPU_pilot_110X_mcRun4_realistic_v3-v2/MINIAODSIM'
#config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-NoPU_pilot_111X_mcRun4_realistic_T15_v1-v1/FEVT'
#config.Data.inputDataset = config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/PhaseIITDRSpring19DR-NoPU_pilot_106X_upgrade2023_realistic_v3-v2/GEN-SIM-DIGI-RAW'
# PU 140
#config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2HLTTDRWinter20RECOMiniAOD-PU140_pilot_110X_mcRun4_realistic_v3-v2/MINIAODSIM'
config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2HLTTDRSummer20ReRECOMiniAOD-PU140_pilot_111X_mcRun4_realistic_T15_v1-v1/FEVT'
# PU 200
#config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2HLTTDRWinter20RECOMiniAOD-PU200_pilot_110X_mcRun4_realistic_v3-v2/MINIAODSIM'



#config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader'

config.Data.outputDatasetTag   = tag

