# submit_crab.py
import copy
from CRABClient.UserUtilities import config
from CRABAPI.RawCommand       import crabCommand


baseTag = "muonHLT_phase2_DYToLL"
datasets = {
    'noPU' : '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2Spring24DIGIRECOMiniAOD-noPU_Trk1GeV_140X_mcRun4_realistic_v4-v1/GEN-SIM-DIGI-RAW-MINIAOD',
    'PU200': '/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2Spring24DIGIRECOMiniAOD-PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/GEN-SIM-DIGI-RAW-MINIAOD'
}

for puLabel, ds in datasets.items():
    cfg = config()

    #–– General ––#
    cfg.section_('General')
    cfg.General.requestName    = f"{baseTag}_{puLabel}"
    cfg.General.workArea       = f"crab_{baseTag}_{puLabel}"
    cfg.General.transferLogs   = True

    #–– JobType ––#
    cfg.section_('JobType')
    cfg.JobType.pluginName     = 'Analysis'
    cfg.JobType.psetName       = 'Phase2Config_15July2025_HLT.py'
    cfg.JobType.outputFiles    = ['muonNtuple_phase2_MC.root']
    cfg.JobType.numCores       = 1
    cfg.JobType.maxMemoryMB    = 5000
    cfg.JobType.allowUndistributedCMSSW = True

    #–– Data ––#
    cfg.section_('Data')
    cfg.Data.inputDataset      = ds
    cfg.Data.splitting         = 'EventAwareLumiBased'
    cfg.Data.unitsPerJob       = 100
    cfg.Data.totalUnits        = 100000
    cfg.Data.useParent         = True
    cfg.Data.outputDatasetTag  = f"{baseTag}_{puLabel}"
    cfg.Data.outLFNDirBase     = f"/store/user/rasharma/MuonHLT/{baseTag}_{puLabel}"

    #–– Site ––#
    cfg.section_('Site')
    cfg.Site.storageSite       = 'T2_US_Purdue'

    print(f"Submitting {puLabel} → requestName = {cfg.General.requestName}")
    crabCommand('submit', config=cfg)
