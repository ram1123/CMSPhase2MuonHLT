import FWCore.ParameterSet.Config as cms

def loadPhase2MuonHLTPaths(process, processName = "MYHLT"):
    # N.B. this customizer overwrites previously defined filters and paths

    print """
########################################################
Imported muon paths:
    - L1_SingleTkMuon_22
    - L1_DoubleTkMuon_17_8
    - L1_TripleTkMuon_5_3_3
    - HLT_Mu50_FromL1TkMuon
    - HLT_Mu50_FromL1TkMuon_Open
    - HLT_IsoMu24_FromL1TkMuon
    - HLT_Mu37_Mu27_FromL1TkMuon
    - HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_FromL1TkMuon
    - HLT_TriMu_10_5_5_DZ_FromL1TkMuon

L1TkMuon scalings (ScalingsV11p1.txt):
    L1TkMuonScalings = cms.PSet(
        barrel=cms.vdouble(0.820128, 1.04124, 0.0),
        overlap=cms.vdouble(0.920897, 1.03712, 0.0),
        endcap=cms.vdouble(0.864715, 1.03215, 0.0),
    )
########################################################
"""

    # -- L1TkMuon scalings and allowed qualities
    L1TkMuonScalings = cms.PSet(
        barrel=cms.vdouble(0.820128, 1.04124, 0.0),
        overlap=cms.vdouble(0.920897, 1.03712, 0.0),
        endcap=cms.vdouble(0.864715, 1.03215, 0.0),
    )

    # -- L3 Reconstruction and Isolation
    from HLTrigger.PhaseII.Muon.Customizers.customizerForPhase2MuonHLT_OIFromL2_VectorHits import customizePhase2MuonHLTReconstruction
    process = customizePhase2MuonHLTReconstruction(process)

    from HLTrigger.PhaseII.Muon.Customizers.customizerForPhase2MuonHLT_OIFromL2_VectorHits import customizePhase2MuonHLTIsolation
    process = customizePhase2MuonHLTIsolation(process)

    for moduleType in [process.producers_(), process.filters_(), process.analyzers_()]:
        for name, module in moduleType.iteritems():
            if hasattr(module, "mightGet") and module.mightGet:
                module.mightGet = cms.optional.untracked.vstring

    # -- Single muon filters
    process.hltL1TkSingleMuFiltered22 = cms.EDFilter("L1TTkMuonFilter",
        saveTags = cms.bool( True ),
        MinPt = cms.double( 22.0 ),
        MinN=cms.int32(1),
        MinEta=cms.double(-2.4),
        MaxEta=cms.double(2.4),
        inputTag = cms.InputTag("L1TkMuons", "", processName),
        Scalings = L1TkMuonScalings
    )
    process.hltL3fL1TkSingleMu22L3Filtered50Q = cms.EDFilter( "HLTMuonTrkL1TkMuFilter",
        saveTags = cms.bool( True ),
        maxNormalizedChi2 = cms.double( 1.0E99 ),
        maxAbsEta = cms.double( 1.0E99 ),
        minPt = cms.double( 50.0 ),
        minN = cms.uint32( 1 ),
        minMuonStations = cms.int32( 1 ),
        minMuonHits = cms.int32( -1 ),
        minTrkHits = cms.int32( -1 ),
        previousCandTag = cms.InputTag( "hltL1TkSingleMuFiltered22" ),
        inputMuonCollection = cms.InputTag( "hltPhase2L3Muons" ),
        inputCandCollection = cms.InputTag( "hltPhase2L3MuonCandidates" )
    )
    process.hltL3fL1TkSingleMu22L3Filtered24Q = cms.EDFilter( "HLTMuonTrkL1TkMuFilter",
        saveTags = cms.bool( True ),
        maxNormalizedChi2 = cms.double( 1.0E99 ),
        maxAbsEta = cms.double( 1.0E99 ),
        minPt = cms.double( 24.0 ),
        minN = cms.uint32( 1 ),
        minMuonStations = cms.int32( 1 ),
        minMuonHits = cms.int32( -1 ),
        minTrkHits = cms.int32( -1 ),
        previousCandTag = cms.InputTag( "hltL1TkSingleMuFiltered22" ),
        inputMuonCollection = cms.InputTag( "hltPhase2L3Muons" ),
        inputCandCollection = cms.InputTag( "hltPhase2L3MuonCandidates" )
    )
    process.hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41 = cms.EDFilter( "HLTMuonGenericFilter",
        thrOverE2EE = cms.vdouble( -1.0 ),
        effectiveAreas = cms.vdouble( 0.0, 0.0 ),
        energyLowEdges = cms.vdouble( 0.0 ),
        doRhoCorrection = cms.bool( False ),
        saveTags = cms.bool( True ),
        thrOverE2EB = cms.vdouble( -1.0 ),
        thrRegularEE = cms.vdouble( -1.0 ),
        thrOverEEE = cms.vdouble( 0.41 ),
        varTag = cms.InputTag( "hltPhase2L3MuonsEcalIsodR0p3dRVeto0p000" ),
        thrOverEEB = cms.vdouble( 0.41 ),
        thrRegularEB = cms.vdouble( -1.0 ),
        lessThan = cms.bool( True ),
        l1EGCand = cms.InputTag( "hltPhase2L3MuonCandidates" ),
        ncandcut = cms.int32( 1 ),
        absEtaLowEdges = cms.vdouble( 0.0, 1.479 ),
        candTag = cms.InputTag( "hltL3fL1TkSingleMu22L3Filtered24Q" ),
        rhoTag = cms.InputTag( "" ),
        rhoMax = cms.double( 9.9999999E7 ),
        useEt = cms.bool( True ),
        rhoScale = cms.double( 1.0 )
    )
    process.hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40 = cms.EDFilter( "HLTMuonGenericFilter",
        thrOverE2EE = cms.vdouble( -1.0 ),
        effectiveAreas = cms.vdouble( 0.0, 0.0 ),
        energyLowEdges = cms.vdouble( 0.0 ),
        doRhoCorrection = cms.bool( False ),
        saveTags = cms.bool( True ),
        thrOverE2EB = cms.vdouble( -1.0 ),
        thrRegularEE = cms.vdouble( -1.0 ),
        thrOverEEE = cms.vdouble( 0.4 ),
        varTag = cms.InputTag( "hltPhase2L3MuonsHcalIsodR0p3dRVeto0p000" ),
        thrOverEEB = cms.vdouble( 0.4 ),
        thrRegularEB = cms.vdouble( -1.0 ),
        lessThan = cms.bool( True ),
        l1EGCand = cms.InputTag( "hltPhase2L3MuonCandidates" ),
        ncandcut = cms.int32( 1 ),
        absEtaLowEdges = cms.vdouble( 0.0, 1.479 ),
        candTag = cms.InputTag( "hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41" ),
        rhoTag = cms.InputTag( "" ),
        rhoMax = cms.double( 9.9999999E7 ),
        useEt = cms.bool( True ),
        rhoScale = cms.double( 1.0 )
    )
    process.hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70 = cms.EDFilter("HLTMuonGenericFilter",
        absEtaLowEdges = cms.vdouble(0.0, 1.479),
        candTag = cms.InputTag("hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40"),
        doRhoCorrection = cms.bool(False),
        effectiveAreas = cms.vdouble(0.0, 0.0),
        energyLowEdges = cms.vdouble(0.0),
        l1EGCand = cms.InputTag("hltPhase2L3MuonCandidates"),
        lessThan = cms.bool(True),
        ncandcut = cms.int32(1),
        rhoMax = cms.double(99999999.0),
        rhoScale = cms.double(1.0),
        rhoTag = cms.InputTag(""),
        saveTags = cms.bool(True),
        thrOverE2EB = cms.vdouble(-1.0),
        thrOverE2EE = cms.vdouble(-1.0),
        thrOverEEB = cms.vdouble(4.70),
        thrOverEEE = cms.vdouble(4.70),
        thrRegularEB = cms.vdouble(-1.0),
        thrRegularEE = cms.vdouble(-1.0),
        useEt = cms.bool(True),
        varTag = cms.InputTag("hltPhase2L3MuonsHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00")
    )
    process.hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk = cms.EDFilter( "HLTMuonIsoFilter",
        saveTags = cms.bool( True ),
        PreviousCandTag = cms.InputTag( "hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70" ),
        MinN = cms.int32( 1 ),
        IsolatorPSet = cms.PSet(  ),
        CandTag = cms.InputTag( "hltPhase2L3MuonCandidates" ),
        DepTag = cms.VInputTag( 'hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p07' )
    )

    # -- Double muon filters
    process.hltL1TkDoubleMuFiltered7 = cms.EDFilter("L1TTkMuonFilter",
        saveTags = cms.bool( True ),
        MinPt=cms.double(7.0),
        MinN=cms.int32(2),
        MinEta=cms.double(-2.4),
        MaxEta=cms.double(2.4),
        inputTag = cms.InputTag("L1TkMuons", "", processName),
        Scalings = L1TkMuonScalings
    )
    process.hltL1TkSingleMuFiltered15 = cms.EDFilter("L1TTkMuonFilter",
        saveTags = cms.bool( True ),
        MinPt=cms.double(15.0),
        MinN=cms.int32(1),
        MinEta=cms.double(-2.4),
        MaxEta=cms.double(2.4),
        inputTag = cms.InputTag("L1TkMuons", "", processName),
        Scalings = L1TkMuonScalings
    )
    process.hltDoubleMuon7DZ1p0 = cms.EDFilter("HLT2L1TkMuonL1TkMuonDZ",
        saveTags = cms.bool( True ),
        originTag1 = cms.VInputTag( 'L1TkMuons::%s' % processName ),
        originTag2 = cms.VInputTag( 'L1TkMuons::%s' % processName ),
        MinPixHitsForDZ = cms.int32( 0 ),
        MinN = cms.int32( 1 ),
        triggerType1 = cms.int32( -114 ),
        triggerType2 = cms.int32( -114 ),
        MinDR = cms.double( -1 ),
        MaxDZ = cms.double( 1.0 ),
        inputTag1 = cms.InputTag( "hltL1TkDoubleMuFiltered7" ),
        checkSC = cms.bool( False ),
        inputTag2 = cms.InputTag( "hltL1TkDoubleMuFiltered7" )
    )
    process.hltDoubleMuon7DR0 = cms.EDFilter("HLT2L1TkMuonL1TkMuonMuRefDR",
        saveTags = cms.bool( True ),
        originTag1 = cms.VInputTag( 'L1TkMuons::%s' % processName ),
        originTag2 = cms.VInputTag( 'L1TkMuons::%s' % processName ),
        MinN = cms.int32( 1 ),
        MinDR = cms.double( 0 ),
        inputTag1 = cms.InputTag( "hltL1TkDoubleMuFiltered7" ),
        inputTag2 = cms.InputTag( "hltL1TkDoubleMuFiltered7" )
    )
    process.hltL3fL1DoubleMu155fPreFiltered8 = cms.EDFilter( "HLTMuonTrkL1TkMuFilter",
        saveTags = cms.bool( True ),
        maxNormalizedChi2 = cms.double( 1.0E99 ),
        maxAbsEta = cms.double( 1.0E99 ),
        minPt = cms.double( 8.0 ),
        minN = cms.uint32( 2 ),
        minMuonStations = cms.int32( 1 ),
        minMuonHits = cms.int32( -1 ),
        minTrkHits = cms.int32( -1 ),
        previousCandTag = cms.InputTag( "hltL1TkDoubleMuFiltered7" ),
        inputMuonCollection = cms.InputTag( "hltPhase2L3Muons" ),
        inputCandCollection = cms.InputTag( "hltPhase2L3MuonCandidates" )
    )
    process.hltL3fL1DoubleMu155fPreFiltered27 = cms.EDFilter( "HLTMuonTrkL1TkMuFilter",
        saveTags = cms.bool( True ),
        maxNormalizedChi2 = cms.double( 1.0E99 ),
        maxAbsEta = cms.double( 1.0E99 ),
        minPt = cms.double( 27.0 ),
        minN = cms.uint32( 2 ),
        minMuonStations = cms.int32( 1 ),
        minMuonHits = cms.int32( -1 ),
        minTrkHits = cms.int32( -1 ),
        previousCandTag = cms.InputTag( "hltL1TkDoubleMuFiltered7" ),
        inputMuonCollection = cms.InputTag( "hltPhase2L3Muons" ),
        inputCandCollection = cms.InputTag( "hltPhase2L3MuonCandidates" )
    )
    process.hltL3fL1DoubleMu155fFiltered17 = cms.EDFilter( "HLTMuonTrkL1TkMuFilter",
        saveTags = cms.bool( True ),
        maxNormalizedChi2 = cms.double( 1.0E99 ),
        maxAbsEta = cms.double( 1.0E99 ),
        minPt = cms.double( 17.0 ),
        minN = cms.uint32( 1 ),
        minMuonStations = cms.int32( 1 ),
        minMuonHits = cms.int32( -1 ),
        minTrkHits = cms.int32( -1 ),
        previousCandTag = cms.InputTag( "hltL1TkSingleMuFiltered15" ),
        inputMuonCollection = cms.InputTag( "hltPhase2L3Muons" ),
        inputCandCollection = cms.InputTag( "hltPhase2L3MuonCandidates" )
    )
    process.hltL3fL1DoubleMu155fFiltered37 = cms.EDFilter( "HLTMuonTrkL1TkMuFilter",
        saveTags = cms.bool( True ),
        maxNormalizedChi2 = cms.double( 1.0E99 ),
        maxAbsEta = cms.double( 1.0E99 ),
        minPt = cms.double( 37.0 ),
        minN = cms.uint32( 1 ),
        minMuonStations = cms.int32( 1 ),
        minMuonHits = cms.int32( -1 ),
        minTrkHits = cms.int32( -1 ),
        previousCandTag = cms.InputTag( "hltL1TkSingleMuFiltered15" ),
        inputMuonCollection = cms.InputTag( "hltPhase2L3Muons" ),
        inputCandCollection = cms.InputTag( "hltPhase2L3MuonCandidates" )
    )
    process.hltDiMuon178RelTrkIsoFiltered0p4 = cms.EDFilter( "HLTMuonIsoFilter",
        saveTags = cms.bool( True ),
        PreviousCandTag = cms.InputTag( "hltL3fL1DoubleMu155fPreFiltered8" ),
        MinN = cms.int32( 2 ),
        IsolatorPSet = cms.PSet(  ),
        CandTag = cms.InputTag( "hltPhase2L3MuonCandidates" ),
        DepTag = cms.VInputTag( 'hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p4' )
    )
    process.hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2 = cms.EDFilter( "HLT2MuonMuonDZ",
        saveTags = cms.bool( True ),
        originTag1 = cms.VInputTag( 'hltPhase2L3MuonCandidates' ),
        originTag2 = cms.VInputTag( 'hltPhase2L3MuonCandidates' ),
        MinPixHitsForDZ = cms.int32( 0 ),
        MinN = cms.int32( 1 ),
        triggerType1 = cms.int32( 83 ),
        triggerType2 = cms.int32( 83 ),
        MinDR = cms.double( 0.001 ),
        MaxDZ = cms.double( 0.2 ),
        inputTag1 = cms.InputTag( "hltDiMuon178RelTrkIsoFiltered0p4" ),
        checkSC = cms.bool( False ),
        inputTag2 = cms.InputTag( "hltDiMuon178RelTrkIsoFiltered0p4" )
    )
    # Not working yet
        # process.hltDiMuon178Mass3p8Filtered = cms.EDFilter( "HLTMuonDimuonL3Filter",
        #     saveTags = cms.bool( True ),
        #     ChargeOpt = cms.int32( 0 ),
        #     MaxPtMin = cms.vdouble( 1.0E125 ),
        #     FastAccept = cms.bool( False ),
        #     MatchToPreviousCand = cms.bool( True ),
        #     CandTag = cms.InputTag( "hltPhase2L3MuonCandidates" ),
        #     L1CandTag = cms.InputTag( "" ),
        #     inputMuonCollection = cms.InputTag( "hltPhase2L3Muons" ),
        #     InputLinks = cms.InputTag( "hltL3MuonsPhase2L3Links" ),
        #     PreviousCandIsL2 = cms.bool( False ),
        #     PreviousCandTag = cms.InputTag( "hltDiMuon178RelTrkIsoFiltered0p4" ),
        #     MaxPtBalance = cms.double( 999999.0 ),
        #     MaxPtPair = cms.vdouble( 1.0E125 ),
        #     MaxAcop = cms.double( 9999.0 ),
        #     MinPtMin = cms.vdouble( 8.0 ),
        #     MaxInvMass = cms.vdouble( 9999.0 ),
        #     MinPtMax = cms.vdouble( 0.0 ),
        #     BeamSpotTag = cms.InputTag( "hltOnlineBeamSpot" ),
        #     MinN = cms.int32( 1 ),
        #     MaxDz = cms.double( 9999.0 ),
        #     MinPtPair = cms.vdouble( 0.0 ),
        #     MaxDr = cms.double( 2.0 ),
        #     MinAcop = cms.double( -1.0 ),
        #     MaxDCAMuMu = cms.double( 99999.9 ),
        #     MinNhits = cms.int32( 0 ),
        #     NSigmaPt = cms.double( 0.0 ),
        #     MinPtBalance = cms.double( -1.0 ),
        #     MaxEta = cms.double( 2.5 ),
        #     L1MatchingdR = cms.double( 0.3 ),
        #     MaxRapidityPair = cms.double( 999999.0 ),
        #     CutCowboys = cms.bool( False ),
        #     MinInvMass = cms.vdouble( 3.8 )
        # )

    # -- Triple muon filters
    process.hltL1TripleMuFiltered3 = cms.EDFilter("L1TTkMuonFilter",
        saveTags = cms.bool( True ),
        MinPt=cms.double(3.0),
        MinN=cms.int32(3),
        MinEta=cms.double(-2.4),
        MaxEta=cms.double(2.4),
        inputTag = cms.InputTag("L1TkMuons", "", processName),
        Scalings = L1TkMuonScalings
    )
    process.hltL1SingleMuFiltered5 = cms.EDFilter("L1TTkMuonFilter",
        saveTags = cms.bool( True ),
        MinPt=cms.double(5.0),
        MinN=cms.int32(1),
        MinEta=cms.double(-2.4),
        MaxEta=cms.double(2.4),
        inputTag = cms.InputTag("L1TkMuons", "", processName),
        Scalings = L1TkMuonScalings
    )
    process.hltTripleMuon3DZ1p0 = cms.EDFilter("HLT2L1TkMuonL1TkMuonDZ",
        saveTags = cms.bool( True ),
        originTag1 = cms.VInputTag( 'L1TkMuons::%s' % processName ),
        originTag2 = cms.VInputTag( 'L1TkMuons::%s' % processName ),
        MinPixHitsForDZ = cms.int32( 0 ),
        MinN = cms.int32( 3 ),
        triggerType1 = cms.int32( -114 ),
        triggerType2 = cms.int32( -114 ),
        MinDR = cms.double( -1 ),
        MaxDZ = cms.double( 1.0 ),
        inputTag1 = cms.InputTag( "hltL1TripleMuFiltered3" ),
        checkSC = cms.bool( False ),
        inputTag2 = cms.InputTag( "hltL1TripleMuFiltered3" )
    )
    process.hltTripleMuon3DR0 = cms.EDFilter("HLT2L1TkMuonL1TkMuonMuRefDR",
        saveTags = cms.bool( True ),
        originTag1 = cms.VInputTag( 'L1TkMuons::%s' % processName ),
        originTag2 = cms.VInputTag( 'L1TkMuons::%s' % processName ),
        MinN = cms.int32( 3 ),
        MinDR = cms.double( 0 ),
        inputTag1 = cms.InputTag( "hltL1TripleMuFiltered3" ),
        inputTag2 = cms.InputTag( "hltL1TripleMuFiltered3" )
    )
    process.hltL3fL1TkTripleMu533PreFiltered555 = cms.EDFilter( "HLTMuonTrkL1TkMuFilter",
        saveTags = cms.bool( True ),
        maxNormalizedChi2 = cms.double( 1.0E99 ),
        maxAbsEta = cms.double( 2.5 ),
        minPt = cms.double( 5.0 ),
        minN = cms.uint32( 3 ),
        minMuonStations = cms.int32( 1 ),
        minMuonHits = cms.int32( -1 ),
        minTrkHits = cms.int32( -1 ),
        previousCandTag = cms.InputTag( "hltL1TripleMuFiltered3" ),
        inputMuonCollection = cms.InputTag( "hltPhase2L3Muons" ),
        inputCandCollection = cms.InputTag( "hltPhase2L3MuonCandidates" )
    )
    process.hltL3fL1TkTripleMu533L3Filtered1055 = cms.EDFilter( "HLTMuonTrkL1TkMuFilter",
        saveTags = cms.bool( True ),
        maxNormalizedChi2 = cms.double( 1.0E99 ),
        maxAbsEta = cms.double( 2.5 ),
        minPt = cms.double( 10.0 ),
        minN = cms.uint32( 1 ),
        minMuonStations = cms.int32( 1 ),
        minMuonHits = cms.int32( -1 ),
        minTrkHits = cms.int32( -1 ),
        previousCandTag = cms.InputTag( "hltL1TripleMuFiltered3" ),  # hltL1SingleMuFiltered5
        inputMuonCollection = cms.InputTag( "hltPhase2L3Muons" ),
        inputCandCollection = cms.InputTag( "hltPhase2L3MuonCandidates" )
    )
    process.hltL3fL1TkTripleMu533L31055DZFiltered0p2 = cms.EDFilter( "HLT2MuonMuonDZ",
        saveTags = cms.bool( True ),
        originTag1 = cms.VInputTag( 'hltPhase2L3MuonCandidates' ),
        originTag2 = cms.VInputTag( 'hltPhase2L3MuonCandidates' ),
        MinPixHitsForDZ = cms.int32( 1 ),
        MinN = cms.int32( 3 ),
        triggerType1 = cms.int32( 83 ),
        triggerType2 = cms.int32( 83 ),
        MinDR = cms.double( 0.001 ),
        MaxDZ = cms.double( 0.2 ),
        inputTag1 = cms.InputTag( "hltL3fL1TkTripleMu533PreFiltered555" ),
        checkSC = cms.bool( False ),
        inputTag2 = cms.InputTag( "hltL3fL1TkTripleMu533PreFiltered555" )
    )

    # -- Paths
    process.L1_SingleTkMuon_22 = cms.Path(
        process.HLTBeginSequence+
        process.hltL1TkSingleMuFiltered22+
        process.HLTEndSequence
    )

    process.L1_DoubleTkMuon_17_8 = cms.Path(
        process.HLTBeginSequence+
        process.hltL1TkDoubleMuFiltered7+
        process.hltL1TkSingleMuFiltered15+
        process.hltDoubleMuon7DZ1p0+
        process.HLTEndSequence
    )

    process.L1_TripleTkMuon_5_3_3 = cms.Path(
        process.HLTBeginSequence+
        process.hltL1TripleMuFiltered3+
        process.hltL1SingleMuFiltered5+
        process.hltTripleMuon3DZ1p0+
        process.hltTripleMuon3DR0+
        process.HLTEndSequence
    )

    process.HLT_Mu50_FromL1TkMuon = cms.Path(
        process.HLTBeginSequence+

        # L1TkMuon filter
        process.hltL1TkSingleMuFiltered22+

        # local reco
        process.HLTMuonLocalRecoSequence+
        process.HLTDoLocalPixelSequence+
        process.HLTDoLocalStripSequence+

        # L2 + L3 reco
        process.HLTL2muonrecoSequence+
        process.HLTPhase2L3MuonRecoSequence+

        process.hltL3fL1TkSingleMu22L3Filtered50Q+

        process.HLTEndSequence
    )

    process.HLT_Mu50_FromL1TkMuon_Open = cms.Path(
        process.HLTBeginSequence+

        # L1TkMuon filter
        cms.ignore( process.hltL1TkSingleMuFiltered22 )+

        # local reco
        process.HLTMuonLocalRecoSequence+
        process.HLTDoLocalPixelSequence+
        process.HLTDoLocalStripSequence+

        # L2 + L3 reco
        process.HLTL2muonrecoSequence+
        process.HLTPhase2L3MuonRecoSequence+

        cms.ignore( process.hltL3fL1TkSingleMu22L3Filtered50Q )+

        process.HLTEndSequence
    )

    process.HLT_IsoMu24_FromL1TkMuon = cms.Path(
        process.HLTBeginSequence+

        # L1TkMuon filter
        process.hltL1TkSingleMuFiltered22+

        # local reco
        process.HLTMuonLocalRecoSequence+
        process.HLTDoLocalPixelSequence+
        process.HLTDoLocalStripSequence+

        # L2 + L3 reco
        process.HLTL2muonrecoSequence+
        process.HLTPhase2L3MuonRecoSequence+

        process.hltL3fL1TkSingleMu22L3Filtered24Q+

        # Isolation
        process.HLTDoFullUnpackingEgammaEcalSequence+
        process.HLTDoLocalHcalSequence+
        process.HLTFastJetForMuons+

        process.HLTPFClusteringForMuonsUnseeded+
        process.hltPhase2L3MuonsEcalIsodR0p3dRVeto0p000+
        process.hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41+

        process.HLTPFHcalClusteringForMuons+
        process.hltPhase2L3MuonsHcalIsodR0p3dRVeto0p000+
        process.hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40+

        process.HLTHgcalTiclLayerClusteringForMuons + 
        process.hltPhase2L3MuonsHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00+
        process.hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70+

        process.tracking_v6_1_L3Muon_NoVertexReco+
        process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p07+
        process.hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk+

        process.HLTEndSequence
    )

    process.HLT_Mu37_Mu27_FromL1TkMuon = cms.Path( 
        process.HLTBeginSequence+

        # L1TkMuon filter
        process.hltL1TkDoubleMuFiltered7+
        process.hltL1TkSingleMuFiltered15+
        process.hltDoubleMuon7DZ1p0+

        # local reco
        process.HLTMuonLocalRecoSequence+
        process.HLTDoLocalPixelSequence+
        process.HLTDoLocalStripSequence+

        # L2 + L3 reco
        process.HLTL2muonrecoSequence+
        process.HLTPhase2L3MuonRecoSequence+

        process.hltL3fL1DoubleMu155fPreFiltered27+
        process.hltL3fL1DoubleMu155fFiltered37+

        process.HLTEndSequence
    )

    process.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_FromL1TkMuon = cms.Path( 
        process.HLTBeginSequence+

        # L1TkMuon filter
        process.hltL1TkDoubleMuFiltered7+
        process.hltL1TkSingleMuFiltered15+
        process.hltDoubleMuon7DZ1p0+

        # local reco
        process.HLTMuonLocalRecoSequence+
        process.HLTDoLocalPixelSequence+
        process.HLTDoLocalStripSequence+

        # L2 + L3 reco
        process.HLTL2muonrecoSequence+
        process.HLTPhase2L3MuonRecoSequence+

        process.hltL3fL1DoubleMu155fPreFiltered8+
        process.hltL3fL1DoubleMu155fFiltered17+

        # Isolation
        process.tracking_v6_1_L3Muon_NoVertexReco+
        process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p4+
        process.hltDiMuon178RelTrkIsoFiltered0p4+

        # Dimuon dz and mass
        process.hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2+
        # process.hltDiMuon178Mass3p8Filtered+

        process.HLTEndSequence
    )

    process.HLT_TriMu_10_5_5_DZ_FromL1TkMuon = cms.Path(
        process.HLTBeginSequence+

        # L1TkMuon filters
        process.hltL1TripleMuFiltered3+
        process.hltL1SingleMuFiltered5+
        process.hltTripleMuon3DZ1p0+
        process.hltTripleMuon3DR0+

        # local reco
        process.HLTMuonLocalRecoSequence+
        process.HLTDoLocalPixelSequence+
        process.HLTDoLocalStripSequence+

        # L2 + L3 reco
        process.HLTL2muonrecoSequence+
        process.HLTPhase2L3MuonRecoSequence+

        # L3 filters
        process.hltL3fL1TkTripleMu533PreFiltered555+
        process.hltL3fL1TkTripleMu533L3Filtered1055+
        process.hltL3fL1TkTripleMu533L31055DZFiltered0p2+

        process.HLTEndSequence
    )

    return process

