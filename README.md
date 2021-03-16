# CMS Phase-2 Muon HLT

## Running L3 Muon reco (OI + IO) in CMSSW_11_3_0_pre3

### Setup
```shell
cmsrel CMSSW_11_3_0_pre3
cd CMSSW_11_3_0_pre3/src
cmsenv

git cms-init
git cms-addpkg RecoMuon/TrackerSeedGenerator
git cms-addpkg HLTrigger/HLTfilters

git clone https://github.com/kondratyevd/CMSPhase2MuonHLT.git HLTrigger/PhaseII/python/Muon

cp -r HLTrigger/PhaseII/python/Muon/Analyzers HLTrigger/
cp -r HLTrigger/PhaseII/python/Muon/TSG/* RecoMuon/TrackerSeedGenerator/plugins/

scram b -j 10
```

### Configs
Located under `HLTrigger/PhaseII/python/Muon/example_cfgs/`
```shell
HLT_test.py # default
HLT_test_VectorHits.py # hit-based seeds in OI will be created using VectorHits
HLT_test_L1TkMu.py # use L1TkMuons to create seeds (skip L2 step)
HLT_test_L1TkMu_VectorHits.py # L1TkMu-based setup with VectorHits
```


