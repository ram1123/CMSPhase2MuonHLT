# from: [Configuration/Geometry/README.md](https://github.com/cms-sw/cmssw/blob/master/Configuration/Geometry/README.md)

- T32: Phase2 tilted tracker. The tracker description is identical to T25. The outer radius of the tracker volume is reduced to avoid a clash with the BTL geometry (same as T31). The positions of the tracker components are not affected. This geometry is intended as a transition step towards a realistic configuration with 3D sensors in TBPX layer1.
- T33: Phase2 tilted tracker. Identical to T32 apart from a more realistic description of the 3D sensors in TBPX layer1.
- Module where these geometries are defined: [Configuration/AlCa/python/autoCondPhase2.py](https://github.com/cms-sw/cmssw/blob/CMSSW_14_0_9_patchX/Configuration/AlCa/python/autoCondPhase2.py)

# ERA

```python
from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
```


# cmsDriverCommand

```bash
cmsDriver.py Phase2Config_15July2025 \
  --step HLT:75e33 \
  --processName=MYHLT \
  --conditions auto:phase2_realistic_T33 \
  --geometry Extended2026D110 \
  --era Phase2C17I13M9 \
  --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000 \
  --eventcontent FEVTDEBUGHLT \
  --filein 'root://cms-xcache.rcac.purdue.edu//store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/.../d47130cf-22dd-4d2b-a458-93981dbf480b.root' \
  --inputCommands 'keep *','drop *_hlt*_*_HLT','drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT' \
  -n 100 \
  --nThreads 1 \
  --no_exec
```

Using the above command, we will obtain a configuration file named `Phase2Config_15July2025.py`. Then connect it with `MuonNtuples.cc` to extract the required information.

To add it to the configuration file append the following lines, after the standard process.load() statement:

```python
from RecoMuon.TrackingTools.MuonServiceProxy_cff import MuonServiceProxy

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

process.HLTValidation = cms.EndPath( process.muonNtuples )
```

Then schedule the path in the configuration file:

```python
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
```

- Use of `FEVTDEBUGHLToutput`: Writes out a full‐event EDM file (raw→digi→reco→HLT products) in CMSSW format. If we need to re-run HLT inspection downstream, or to inspect individual reco objects in a standalone job) → then add  the FEVTDEBUGHLToutput PoolOutputModule on an EndPath.
- We need simple root file so comment out the `FEVTDEBUGHLToutput` and utilize the `TFileService` to write the output to a ROOT file.

```python
process.TFileService = cms.Service("TFileService",
                               fileName = cms.string("muonNtuple_phase2_MC.root"),
                               closeFileFast = cms.untracked.bool(False)
)
```
