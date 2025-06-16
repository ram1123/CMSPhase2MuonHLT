# CMS Phase-2 Muon HLT

## Running L3 Muon reco (OI + IO) in CMSSW_11_3_0_pre3

### Setup Instruction

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_14_0_9
cd  CMSSW_14_0_9/src
cmsenv
git cms-init
git cms-addpkg HLTrigger/Configuration
git cms-addpkg Configuration/Geometry
# git cms-addpkg RecoMuon/TrackerSeedGenerator
# git cms-addpkg Configuration/Eras
git clone git@github.com:ram1123/CMSPhase2MuonHLT.git HLTrigger/PhaseII/python/Muon
cp -r HLTrigger/PhaseII/python/Muon/Analyzers HLTrigger/
# cp -r HLTrigger/PhaseII/python/Muon/TSG/* RecoMuon/TrackerSeedGenerator/plugins/
scramv1 b clean
scram b -j 10
```

-  Main code to modify or update hist or variables: `HLTrigger/PhaseII/python/Muon/TSG/*.cc`
    - TSGForOIDNN.cc
    - TSGForOIFromL1TkMu.cc
    - TSGForOIFromL2.cc

- Config files that we need to run: `HLTrigger/PhaseII/python/Muon/example_cfgs/`
    - The configs contain full Phase2 HLT setup, but you can choose whether to create OI seeds from L2 muons or L1TkMuons (skip L2).
    - VectorHits are enabled for both tracking and seeding in setups with suffix `_VectorHits`.


```bash
cd HLTrigger/PhaseII/python/Muon/example_cfgs/
cmsRun Phase2_HLT_AF_L2.py
# CrabConfig_Phase2Config_15July2025_HLT.py
# crab_config.py
# crab_config_L1TkMuons.py
# SingleMuon_Phase2_HLT_AF_L2.py
#single_muon_crab_config.py
# Phase2_HLT_AF_L1TkMuon.py
# L1TkMuon_QCD_crab_config_L1TkMuons.py
# L2_QCD_crab_config.py
# Phase2_HLT_AF.py
# crab_config_HB0_HL3.py
```

#### Crab Job Submission

```bash
cd HLTrigger/PhaseII/python/Muon/example_cfgs/
voms-proxy-init --voms cms --valid 168:00 --out $(pwd)/voms_proxy.txt
export X509_USER_PROXY=$(pwd)/voms_proxy.txt
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit -c CrabConfig_Phase2Config_15July2025_HLT.py
```

## Dataset:
#### Reference: Amandeep presentation: 26 Nov 2024

```bash
dasgoclient --query="dataset=/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2Spring24DIGIRECOMiniAOD*/GEN-SIM-DIGI-RAW-MINIAOD"
/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2Spring24DIGIRECOMiniAOD-PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/GEN-SIM-DIGI-RAW-MINIAOD
/DYToLL_M-50_TuneCP5_14TeV-pythia8/Phase2Spring24DIGIRECOMiniAOD-noPU_Trk1GeV_140X_mcRun4_realistic_v4-v1/GEN-SIM-DIGI-RAW-MINIAOD
```

```
/store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/2810000/d47130cf-22dd-4d2b-a458-93981dbf480b.root

/store/mc/Phase2Spring24DIGIRECOMiniAOD/DYToLL_M-50_TuneCP5_14TeV-pythia8/GEN-SIM-DIGI-RAW-MINIAOD/noPU_Trk1GeV_140X_mcRun4_realistic_v4-v1/2810000/cf329170-79f1-42ea-afca-9a419eae27ba.root
```


# Questions/Doubts

- The geometry: `GeometryExtended2026D110Reco_cff`. Is it the correct one for Phase2 HLT?
- How should I relate geometry and ERA?

# Some references

1. [Geometry version as well as GT info](https://github.com/cms-sw/cmssw/blob/master/Configuration/Geometry/README.md)

