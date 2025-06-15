import sys


from CRABClient.UserUtilities import config

config = config()

config.JobType.pluginName   = 'Analysis'
config.JobType.outputFiles  = ['muonNtuple_phase2_MC.root']

config.Data.unitsPerJob     = 1000
config.Data.totalUnits      = 100000
#config.Data.splitting       = 'LumiBased'
config.Data.splitting       = 'EventAwareLumiBased'
#config.Data.splitting       = 'Automatic'

#config.Data.useParent       = True #!!!!
config.Data.useParent       = False #!!!!

config.Site.storageSite     = 'T2_US_Purdue'
config.JobType.numCores     = 1
config.JobType.maxMemoryMB  = 5000
config.JobType.allowUndistributedCMSSW = True
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
#from httplib import HTTPException

tag = "SingleMuon_muonHLT_L2Muons"

config.General.workArea   = tag
config.Data.outLFNDirBase = '/store/user/amkaur/' + tag

config.JobType.psetName    = 'SingleMuon_Phase2_HLT_AF_L2.py'
config.General.requestName = tag
config.General.transferLogs = True

config.Data.inputDataset = '/SingleMuon_Pt-0To200_Eta-1p4To3p1-gun/Phase2Fall22DRMiniAOD-PU200_125X_mcRun4_realistic_v2-v1/GEN-SIM-DIGI-RAW-MINIAOD'
config.Data.outputDatasetTag   = tag
