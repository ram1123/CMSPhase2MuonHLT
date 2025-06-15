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

tag = "spring24_DisplacedMuons_L2Muons"

config.General.workArea   = tag
config.Data.outLFNDirBase = '/store/user/amkaur/' + tag

config.JobType.psetName    = 'Phase2_HLT_AF_L2.py'
config.General.requestName = tag
config.General.transferLogs = True

config.Data.inputDataset = '/DisplacedMuons_Pt-30To100_Dxy-0To3000-gun/Phase2Spring24DIGIRECOMiniAOD-PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/GEN-SIM-DIGI-RAW-MINIAOD' 

#config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2Spring24DIGIRECOMiniAOD-PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/GEN-SIM-DIGI-RAW-MINIAOD'
#config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2Spring23DIGIRECOMiniAOD-PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/GEN-SIM-DIGI-RAW-MINIAOD'
config.Data.outputDatasetTag   = tag
