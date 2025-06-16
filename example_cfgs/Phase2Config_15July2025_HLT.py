# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: Phase2Config_15July2025 --step HLT:75e33 --processName=MYHLT --conditions auto:phase2_realistic_T33 --geometry Extended2026D110 --era Phase2C17I13M9 --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000 --eventcontent FEVTDEBUGHLT --filein root://cms-xcache.rcac.purdue.edu//store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/.../d47130cf-22dd-4d2b-a458-93981dbf480b.root --inputCommands keep *,drop *_hlt*_*_HLT,drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT -n 100 --nThreads 1 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('MYHLT',Phase2C17I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('HLTrigger.Configuration.HLT_75e33_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring('file:/tmp/rasharma/d47130cf-22dd-4d2b-a458-93981dbf480b.root'),
    # fileNames = cms.untracked.vstring('root://cms-xcache.rcac.purdue.edu//store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/2810000/d47130cf-22dd-4d2b-a458-93981dbf480b.root'),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop *_hlt*_*_HLT',
        'drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

from RecoMuon.TrackingTools.MuonServiceProxy_cff import MuonServiceProxy

# a debugging tool to print out every branch in every event to the log.
# process.dump = cms.EDAnalyzer('EventContentAnalyzer')
process.muonNtuples = cms.EDAnalyzer("MuonNtuples",
                   MuonServiceProxy,
                   offlineVtx             = cms.InputTag("offlineSlimmedPrimaryVertices"),
                   offlineMuons             = cms.InputTag("slimmedMuons"),
                   triggerResult            = cms.untracked.InputTag("TriggerResults::MYHLT"),
                   triggerSummary           = cms.untracked.InputTag("hltTriggerSummaryAOD::MYHLT"),
                   tagTriggerResult         = cms.untracked.InputTag("TriggerResults::MYHLT"),
                   tagTriggerSummary        = cms.untracked.InputTag("hltTriggerSummaryAOD::MYHLT"),
                   triggerProcess   = cms.string("MYHLT"),
                   # triggerProcess   = cms.string("TEST"),
                   # L3Candidates             = cms.untracked.InputTag("hltPhase2L3Muons"),
                   L3Candidates             = cms.untracked.InputTag("hltPhase2L3MuonCandidates"),
                   L3CandidatesNoID         = cms.untracked.InputTag("hltPhase2MuonsNoID"),
                   L2Candidates             = cms.untracked.InputTag("hltL2MuonFromL1TkMuonCandidates"),
                   L1Candidates             = cms.untracked.InputTag('hltGtStage2Digis','Muon'),
                   L1TkCandidates             = cms.untracked.InputTag("l1tTkMuonsGmt", "", "MYHLT"),
                   #L1TkCandidates             = cms.untracked.InputTag("l1tTkMuonsGmt", "", "MYHLT"),
                   #L1TkCandidates             = cms.untracked.InputTag("L1TkMuons", "", "MYHLT"),
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
                   #seedsForOIFromL2         = cms.InputTag("hltPhase2OISeedsFromL2Muons"),
                   seedsForOIFromL2         = cms.InputTag("hltPhase2L3OISeedsFromL2Muons"),
                   theTrajOI                = cms.untracked.InputTag("hltPhase2OITrackCandidates"),
                   simTracks            = cms.untracked.InputTag("mix","MergedTrackTruth", "MYHLT"),
                   propagatorName       = cms.string('PropagatorWithMaterialParabolicMf'),
)

process.HLTValidation = cms.EndPath(
    process.muonNtuples
)

from L1Trigger.Phase2L1GT.l1tGTAlgoBlockProducer_cff import collectAlgorithmPaths
process.schedule = cms.Schedule(
    # First execute my GlobalTrigger emulation (GTemulation_step)
    process.GTemulation_step,
    *collectAlgorithmPaths(process),
    # Then execute every L1T and HLT selection path in the order they were collected
    process.L1T_SingleTkMuon_22,
   process.L1T_DoubleTkMuon_15_7,
   process.L1T_TripleTkMuon_5_3_3,
   process.HLT_Mu50_FromL1TkMuon,
   process.HLT_IsoMu24_FromL1TkMuon,
   process.HLT_Mu37_Mu27_FromL1TkMuon,
   process.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_FromL1TkMuon,
   process.HLT_TriMu_10_5_5_DZ_FromL1TkMuon,
    # Finally run my HLTValidation endpath (which contains your muonNtuples EDAnalyzer)

   process.HLTriggerFinalPath,
   process.HLTValidation
)


process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
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
    modulesToCallForTryToContinue = cms.untracked.vstring(),
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
    annotation = cms.untracked.string('Phase2Config_15July2025 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
#     dataset = cms.untracked.PSet(
#         dataTier = cms.untracked.string(''),
#         filterName = cms.untracked.string('')
#     ),
#     fileName = cms.untracked.string('Phase2Config_15July2025_HLT.root'),
#     outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
#     splitLevel = cms.untracked.int32(0)
# )

process.TFileService = cms.Service("TFileService",
                               fileName = cms.string("muonNtuple_phase2_MC.root"),
                               closeFileFast = cms.untracked.bool(False)
)

# Additional output definition

# Other statements
from HLTrigger.Configuration.CustomConfigs import ProcessName
process = ProcessName(process)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# Path and EndPath definitions
# process.L1T_DoubleNNTau52 = cms.Path(process.HLTL1Sequence+process.hltL1DoubleNNTau52)
# process.L1T_SingleNNTau150 = cms.Path(process.HLTL1Sequence+process.hltL1SingleNNTau150)
process.endjob_step = cms.EndPath(process.endOfProcess)
# process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
# process.schedule imported from cff in HLTrigger.Configuration
# process.schedule.extend([process.endjob_step,process.FEVTDEBUGHLToutput_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

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
