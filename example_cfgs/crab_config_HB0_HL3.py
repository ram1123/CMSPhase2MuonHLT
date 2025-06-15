import sys

import CRABClient
from WMCore.Configuration import Configuration
config = Configuration()

#from CRABClient.UserUtilities import config

#config = config()
config.section_('JobType')
config.JobType.pluginName   = 'Analysis'
config.JobType.outputFiles  = ['muonNtuple_phase2_MC.root']

config.section_('Data')
config.Data.unitsPerJob     = 1000
config.Data.totalUnits      = 100000
#config.Data.splitting       = 'LumiBased'
config.Data.splitting       = 'EventAwareLumiBased'
#config.Data.splitting       = 'Automatic'

#config.Data.useParent       = True #!!!!
config.Data.useParent       = False #!!!!

config.section_('Site')
config.Site.storageSite     = 'T2_US_Purdue'
config.JobType.numCores     = 1
config.JobType.maxMemoryMB  = 5000
config.JobType.allowUndistributedCMSSW = True
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
#from httplib import HTTPException

tag = "muonHLT_seeds_HB_0_HL_3_AF"

config.section_("General")
config.General.workArea   = tag
config.Data.outLFNDirBase = '/store/user/amkaur/' + tag

config.JobType.psetName    = 'python/Phase2_HLT_HB0_HL3.py'
config.General.requestName = tag
config.General.transferLogs = True

config.Data.inputDataset = '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2Spring23DIGIRECOMiniAOD-PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/GEN-SIM-DIGI-RAW-MINIAOD' 
config.Data.outputDatasetTag   = tag
