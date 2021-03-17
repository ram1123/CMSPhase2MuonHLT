
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process("MYHLT", eras.Phase2C9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
#process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.MessageLogger.cerr.threshold = "DEBUG"
#process.MessageLogger.debugModules = ["hltPhase2L3OISeedsFromL1TkMu"]
process.MessageLogger = cms.Service("MessageLogger",
                                    destinations   = cms.untracked.vstring('myDebugOutputFile'),
                                    #cerr           = cms.untracked.PSet(
                                    #    threshold      = cms.untracked.string('ERROR'),
                                    #),
                                    myDebugOutputFile = cms.untracked.PSet(
                                        threshold = cms.untracked.string('DEBUG'),
                                        default = cms.untracked.PSet(limit = cms.untracked.int32(-1)),
                                    ),
                                    debugModules  = cms.untracked.vstring('*')
                                    )

#process.maxEvents = cms.untracked.PSet(
#    input  = cms.untracked.int32(-1),
#    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
#)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(''),
    secondaryFileNames = cms.untracked.vstring()
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

# -- L1 emulation -- #
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.L1TkMuons.L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks", "RECO")
process.L1TkMuons.applyQualityCuts = cms.bool(True)
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
# -- #

# -- HLT paths -- #
from HLTrigger.PhaseII.Muon.Customizers.loadPhase2MuonHLTPaths_OIFromL1TkMu_cfi import loadPhase2MuonHLTPaths
process = loadPhase2MuonHLTPaths(process)

process.dump = cms.EDAnalyzer('EventContentAnalyzer')
#process.MessageLogger = cms.Service("MessageLogger",
#                                    debugModules = cms.untracked.vstring('MultiTrackValidator')
#                                    )

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *

process.muonNtuples = cms.EDAnalyzer("MuonNtuples",
                                     MuonServiceProxy,
                                     #offlineVtx               = cms.InputTag("offlinePrimaryVertices"),
                                     offlineVtx             = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                     #offlineMuons             = cms.InputTag("muons"),
                                     offlineMuons             = cms.InputTag("slimmedMuons"),
                                     triggerResult            = cms.untracked.InputTag("TriggerResults::MYHLT"),
                                     triggerSummary           = cms.untracked.InputTag("hltTriggerSummaryAOD::MYHLT"),
                                     tagTriggerResult         = cms.untracked.InputTag("TriggerResults::HLT"),
                                     tagTriggerSummary        = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
                                     triggerProcess   = cms.string("MYHLT"),
                                     L3Candidates             = cms.untracked.InputTag("hltPhase2L3MuonCandidates"),
                                     L3CandidatesNoID         = cms.untracked.InputTag("hltPhase2L3MuonsNoID"),
                                     #L2Candidates             = cms.untracked.InputTag("hltL2MuonCandidates"),
                                     L2Candidates             = cms.untracked.InputTag("hltL2MuonFromL1TkMuonCandidates"),
                                     L1Candidates             = cms.untracked.InputTag('hltGtStage2Digis','Muon'),
                                     L1TkCandidates             = cms.untracked.InputTag("L1TkMuons", "", "MYHLT"),
                                     #L1TkCandidates             = cms.untracked.InputTag("TTTracksFromTrackletEmulation","Level1TTTracks","RECO"),
                                     TkMuCandidates           = cms.untracked.InputTag("hltPhase2L3OIL3MuonCandidates"),
                                     L3OIMuCandidates         = cms.untracked.InputTag("hltPhase2L3OIL3MuonCandidates"),
                                     L3IOMuCandidates         = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonCkfTrackCandidates"),
                                     MuonLinksTag = cms.untracked.InputTag("hltPhase2L3OIL3MuonsLinksCombination"),
                                     globalMuons = cms.InputTag("globalMuons"),
                                     theTrackOI               = cms.untracked.InputTag("hltPhase2L3OIMuonTrackSelectionHighPurity"),
                                     #theTrackOI               = cms.untracked.InputTag("hltPhase2L3OIMuCtfWithMaterialTracks"),
                                     theTrackIOL2             = cms.untracked.InputTag("hltIter3IterL3MuonMerged"),
                                     theTrackIOL1             = cms.untracked.InputTag("hltIter2Phase2L3FromL1TkMuonMerged"),
                                     #l3filterLabel    = cms.string("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q"),
                                     l3filterLabel    = cms.string("hltDiMuon178Mass3p8Filtered"),
                                     lumiScalerTag            = cms.untracked.InputTag("scalersRawToDigi"),
                                     puInfoTag                = cms.untracked.InputTag("addPileupInfo"),
                                     genParticlesTag          = cms.untracked.InputTag("genParticles"),
                                     doOffline                = cms.untracked.bool(True),
                                     #seedsForOIFromL2         = cms.InputTag("hltPhase2L3OISeedsFromL2Muons"),
                                     seedsForOIFromL2         = cms.InputTag("hltPhase2L3OISeedsFromL1TkMuons"),
                                     theTrajOI                = cms.untracked.InputTag("hltPhase2L3OITrackCandidates"),
                                     simTracks            = cms.untracked.InputTag("mix","MergedTrackTruth", "HLT"), 
                                 )

process.TFileService = cms.Service("TFileService",
                               fileName = cms.string("muonNtuple_phase2_MC.root"),
                               closeFileFast = cms.untracked.bool(False)
)
process.HLTValidation = cms.EndPath(
    process.muonNtuples
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False), SkipEvent = cms.untracked.vstring('ProductNotFound') )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.schedule = cms.Schedule(
    process.L1simulation_step,
    process.L1_SingleTkMuon_22,
    process.L1_DoubleTkMuon_17_8,
    process.L1_TripleTkMuon_5_3_3,
    process.HLT_Mu50_FromL1TkMuon_Open,
    process.HLT_Mu50_FromL1TkMuon,
    process.HLT_IsoMu24_FromL1TkMuon,
    process.HLT_Mu37_Mu27_FromL1TkMuon,
    process.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_FromL1TkMuon,
    process.HLT_TriMu_10_5_5_DZ_FromL1TkMuon,
    process.HLTValidation,
    process.HLTriggerFinalPath
)
# -- #

# -- Test Setup -- #
process.load( "DQMServices.Core.DQMStore_cfi" )
process.DQMStore.enableMultiThread = True

process.GlobalTag.globaltag = "113X_mcRun4_realistic_v1"

process.source.fileNames = cms.untracked.vstring(
    "/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/FEVT/NoPU_pilot_111X_mcRun4_realistic_T15_v1-v1/100000/0018FD1C-F75F-9E4A-A329-1E671A1CA267.root",
    #"file:output.root"
    #"/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/FEVT/PU200_pilot_111X_mcRun4_realistic_T15_v1-v1/270000/FC1C5501-17FF-AD4E-B0C2-78B114D94AD6.root",
    #"/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/FEVT/PU200_pilot_111X_mcRun4_realistic_T15_v1-v1/100000/00695E54-EAD4-3444-A833-3FE1C2BC8880.root",
    #"/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/FEVT/PU200_pilot_111X_mcRun4_realistic_T15_v1-v1/100000/08B5765F-A152-F043-AF98-F84241563451.root",
    #"/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/FEVT/PU200_pilot_111X_mcRun4_realistic_T15_v1-v1/100000/08D8F8B9-6A14-0940-9444-E43ECC3E7D27.root",
    #"/store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/FEVT/PU200_pilot_111X_mcRun4_realistic_T15_v1-v1/100000/0955D65C-03AE-9C4F-91F1-F7979C242772.root",
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( -1 )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True ),
    numberOfThreads = cms.untracked.uint32( 1 ),
    numberOfStreams = cms.untracked.uint32( 0 ),
    sizeOfStackForThreadsInKB = cms.untracked.uint32( 10*1024 )
)

"""
if 'MessageLogger' in process.__dict__:
    process.MessageLogger.categories.append('TriggerSummaryProducerAOD')
    process.MessageLogger.categories.append('L1GtTrigReport')
    process.MessageLogger.categories.append('L1TGlobalSummary')
    process.MessageLogger.categories.append('HLTrigReport')
    process.MessageLogger.categories.append('FastReport')
    process.MessageLogger.cerr.FwkReport.reportEvery = 1  # 1000
# -- #
"""

from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000
process = customise_aging_1000(process)

from L1Trigger.Configuration.customisePhase2TTNoMC import customisePhase2TTNoMC
process = customisePhase2TTNoMC(process)

from HLTrigger.Configuration.Eras import modifyHLTforEras
modifyHLTforEras(process)


# process.Timing = cms.Service("Timing",
#     summaryOnly = cms.untracked.bool(True),
#     useJobReport = cms.untracked.bool(True)
# )

# process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#     ignoreTotal = cms.untracked.int32(1)
# )

# print process.dumpPython()

