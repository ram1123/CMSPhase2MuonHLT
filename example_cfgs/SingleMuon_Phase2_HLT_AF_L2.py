# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Phase2 -s HLT:75e33 --processName=HLTX --conditions auto:phase2_realistic_T21 --geometry Extended2026D95 --era Phase2C17I13M9 --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000 --eventcontent FEVTDEBUGHLT --filein=/store/mc/Phase2Spring23DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/01607282-0427-4687-a122-ef0a41220590.root --inputCommands=keep *, drop *_hlt*_*_HLT, drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT -n 100 --nThreads 1 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('MYHLT',Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D95Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('HLTrigger.Configuration.HLT_75e33_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2000),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

SkipEvent = cms.untracked.vstring('ProductNotFound')

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    #fileNames = cms.untracked.vstring('file:/tmp/009a10c2-f642-4e1d-9dd5-4fbf9c923035.root'),
    fileNames = cms.untracked.vstring('file:/tmp/004ac631-8090-4eac-b87f-224497a31651d.root'),
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/Phase2Fall22DRMiniAOD/SingleMuon_Pt-0To200_Eta-1p4To3p1-gun/GEN-SIM-DIGI-RAW-MINIAOD/PU200_125X_mcRun4_realistic_v2-v1/30000/004ac631-8090-4eac-b87f-24497a31651d.root'),
    #fileNames = cms.untracked.vstring('root://cms-xcache.rcac.purdue.edu///store/mc/Phase2Spring23DIGIRECOMiniAOD/QCD_Pt-15To3000_TuneCP5_Flat_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v2/50000/01297651-d031-41dd-9689-d71b5d1ee0ed.root'),
    #fileNames = cms.untracked.vstring('root://cms-xcache.rcac.purdue.edu//store/mc/Phase2Spring23DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/009a10c2-f642-4e1d-9dd5-4fbf9c923035.root'),
    #fileNames = cms.untracked.vstring('root://cms-xcache.rcac.purdue.edu//store/mc/Phase2Spring23DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/009a10c2-f642-4e1d-9dd5-4fbf9c923035.root'),
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/Phase2Spring23DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/009a10c2-f642-4e1d-9dd5-4fbf9c923035.root'),
    #fileNames = cms.untracked.vstring('/store/mc/Phase2Spring23DIGIRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/30000/01607282-0427-4687-a122-ef0a41220590.root'),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop *_hlt*_*_HLT',
        'drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

# -- L1 emulation -- #
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
#process.L1TkMuons.L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks", "RECO")
#process.L1TkMuons.applyQualityCuts = cms.bool(True)
# -- #

# -- HLTriggerFinalPath -- #
process.hltTriggerSummaryAOD = cms.EDProducer( "TriggerSummaryProducerAOD",
    moduleLabelPatternsToSkip = cms.vstring(  ),
    processName = cms.string( "@" ),
    moduleLabelPatternsToMatch = cms.vstring( 'hlt*', 'L1Tk*' ),
    throw = cms.bool( False )
)
process.hltTriggerSummaryRAW = cms.EDProducer( "TriggerSummaryProducerRAW",
    processName = cms.string( "@" )
)
process.hltBoolFalse = cms.EDFilter( "HLTBool",
    result = cms.bool( False )
)
process.HLTriggerFinalPath = cms.Path(
    process.hltTriggerSummaryAOD+
    process.hltTriggerSummaryRAW+
    process.hltBoolFalse
)

# -- HLT paths -- #
from HLTrigger.PhaseII.Muon.Customizers.loadPhase2MuonHLTPaths_OIFromL2_cfi import loadPhase2MuonHLTPaths
process = loadPhase2MuonHLTPaths(process)

process.dump = cms.EDAnalyzer('EventContentAnalyzer')
#other statements
#from HLTrigger.Configuration.CustomConfigs import ProcessName
#process = ProcessName(process)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *

process.muonNtuples = cms.EDAnalyzer("MuonNtuples",
                   MuonServiceProxy,
                   offlineVtx             = cms.InputTag("offlineSlimmedPrimaryVertices"),
                   offlineMuons             = cms.InputTag("slimmedMuons"),
                   triggerResult            = cms.untracked.InputTag("TriggerResults::MYHLT"),
                   triggerSummary           = cms.untracked.InputTag("hltTriggerSummaryAOD::MYHLT"),
                   tagTriggerResult         = cms.untracked.InputTag("TriggerResults::HLT"),
                   tagTriggerSummary        = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
                   #triggerProcess   = cms.string("MYHLT"),
                   triggerProcess   = cms.string("TEST"),
                   #L3Candidates             = cms.untracked.InputTag("hltPhase2L3Muons"),
                   L3Candidates             = cms.untracked.InputTag("hltPhase2L3MuonCandidates"),
                   L3CandidatesNoID         = cms.untracked.InputTag("hltPhase2MuonsNoID"),
                   L2Candidates             = cms.untracked.InputTag("hltL2MuonFromL1TkMuonCandidates"),
                   L1Candidates             = cms.untracked.InputTag('hltGtStage2Digis','Muon'),
                   #L1TkCandidates             = cms.untracked.InputTag("l1tTrackerMuons_l1tTkMuonsGmt__HLT", "", "HLT"),
                   #L1TkCandidates             = cms.untracked.InputTag("L1TkMuons", "", "MYHLT"),
                   L1TkCandidates             = cms.untracked.InputTag("l1tTkMuonsGmt", "", "HLT"),
                   TkMuCandidates           = cms.untracked.InputTag("hltPhase2L3OIL3MuonCandidates"),
                   L3OIMuCandidates         = cms.untracked.InputTag("hltPhase2L3OIL3MuonCandidates"),
                   L3IOMuCandidates         = cms.untracked.InputTag("hltPhase2IOFromL2MuonCandidates"),
                   MuonLinksTag = cms.untracked.InputTag("hltPhase2MuonsFromL2LinksCombination"),
                   globalMuons = cms.InputTag("globalMuons"),
                   theTrackOI               = cms.untracked.InputTag("hltPhase2L3OIMuonTrackSelectionHighPurity"),
                   theTrackIOL2             = cms.untracked.InputTag("hltIter3Phase2MuonMerged"),
                   theTrackIOL1             = cms.untracked.InputTag("hltIter3Phase2FromL1TkMuonMerged"),
                   l3filterLabel    = cms.string("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q"),
                   lumiScalerTag            = cms.untracked.InputTag("scalersRawToDigi"),
                   puInfoTag                = cms.untracked.InputTag("addPileupInfo"),
                   genParticlesTag          = cms.untracked.InputTag("genParticles"),
                   doOffline                = cms.untracked.bool(True),
                   seedsForOIFromL2         = cms.InputTag("hltPhase2OISeedsFromL2Muons"),
                   theTrajOI                = cms.untracked.InputTag("hltPhase2OITrackCandidates"),
                   simTracks            = cms.untracked.InputTag("mix","MergedTrackTruth", "HLT"),
                   propagatorName       = cms.string('PropagatorWithMaterialParabolicMf'),
)

process.HLTValidation = cms.EndPath(
    process.muonNtuples
)

#Amandeep

# Path and EndPath definitions
#process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step)
#process.schedule.extend(process.HLTSchedule)
#process.schedule.extend([process.endjob_step])

#process.schedule = cms.Schedule(process.HLTSchedule)

#process.schedule.extend([process.endjob_step])
#process.schedule.extend([process.HLTriggerFinalPath])
#process.schedule.extend([process.HLTValidation])

#Amandeep
process.schedule = cms.Schedule(
   #process.L1simulation_step,
   #process.L1_SingleTkMuon_22,
   #process.L1_DoubleTkMuon_17_8,
   #process.L1_TripleTkMuon_5_3_3,
   #process.HLT_Mu50_FromL1TkMuon_Open,
   #process.HLT_Mu50_FromL1TkMuon,
   #process.HLT_IsoMu24_FromL1TkMuon,
   #process.HLT_Mu37_Mu27_FromL1TkMuon,
   #process.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_FromL1TkMuon,
   #process.HLT_TriMu_10_5_5_DZ_FromL1TkMuon,
   process.HLTriggerFinalPath,
   process.HLTValidation
)



process.TFileService = cms.Service("TFileService",
                               fileName = cms.string("muonNtuple_phase2_MC.root"),
                               closeFileFast = cms.untracked.bool(False)
)


process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Phase2 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

#process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
#    dataset = cms.untracked.PSet(
#        dataTier = cms.untracked.string(''),
#        filterName = cms.untracked.string('')
#    ),
#    fileName = cms.untracked.string('Phase2_HLT.root'),
#    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
#    splitLevel = cms.untracked.int32(0)
#)

# Additional output definition

import FWCore.ParameterSet.Config as cms

#
# attention: default is changed to work on unsuppressed digis!! ##############
#
simEcalTriggerPrimitiveDigis = cms.EDProducer("EcalTrigPrimProducer",
    BarrelOnly = cms.bool(False),
    InstanceEB = cms.string(''),
    InstanceEE = cms.string(''),
    binOfMaximum = cms.int32(6), ## optional from release 200 on, from 1-10
    Famos = cms.bool(False),
    TcpOutput = cms.bool(False),
    Debug = cms.bool(False),
    Label = cms.string('simEcalUnsuppressedDigis')
)


# Path and EndPath definitions
#process.L1T_DoubleNNTau52 = cms.Path(process.HLTL1Sequence+process.hltL1DoubleNNTau52)
#process.L1T_SingleNNTau150 = cms.Path(process.HLTL1Sequence+process.hltL1SingleNNTau150)
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
# process.schedule imported from cff in HLTrigger.Configuration
#process.schedule.extend([process.endjob_step,process.FEVTDEBUGHLToutput_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseReEmul
#from L1Trigger.Configuration.customiseReEmul import L1TReEmulMCFromRAW

#call to customisation function L1TReEmulMCFromRAW imported from L1Trigger.Configuration.customiseReEmul
#process = L1TReEmulMCFromRAW(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseSettings
#from L1Trigger.Configuration.customiseSettings import L1TSettingsToCaloParams_2018_v1_4

#call to customisation function L1TSettingsToCaloParams_2018_v1_4 imported from L1Trigger.Configuration.customiseSettings
#process = L1TSettingsToCaloParams_2018_v1_4(process)

# Automatic addition of the customisation function from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3
#from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3 import customizeMuonHLTForDoubletRemoval,customizeMuonHLTForCscSegment,customizeMuonHLTForGEM
#Aman
#from HLTrigger.PhaseII.Muon.Customizers.customizeMuonHLTForRun3 import customizeMuonHLTForDoubletRemoval,customizeMuonHLTForCscSegment,customizeMuonHLTForGEM

#call to customisation function customizeMuonHLTForDoubletRemoval imported from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3
#process = customizeMuonHLTForDoubletRemoval(process)

#call to customisation function customizeMuonHLTForCscSegment imported from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3
#process = customizeMuonHLTForCscSegment(process)

#call to customisation function customizeMuonHLTForGEM imported from HLTrigger.Configuration.MuonHLTForRun3.customizeMuonHLTForRun3
#process = customizeMuonHLTForGEM(process)

# Automatic addition of the customisation function from HLTrigger.Configuration.MuonHLTForRun3.customizeOIseedingForRun3
#from HLTrigger.AutoSubmit.customizer_muonHLTtest_Run3_DYJets_xx3 import customizeOIseeding

#call to customisation function customizeOIseeding imported from HLTrigger.Configuration.MuonHLTForRun3.customizeOIseedingForRun3
#process = customizeOIseeding(process)

#from HLTrigger.PhaseII.Muon.autosubmit.customizer_muonHLTtest_Run3_DYJets_xx3 import customizeOIseeding
#process = customizeOIseeding(process)


# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC


# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.aging
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000 

#call to customisation function customise_aging_1000 imported from SLHCUpgradeSimulations.Configuration.aging
process = customise_aging_1000(process)

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
