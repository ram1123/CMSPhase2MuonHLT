import FWCore.ParameterSet.Config as cms

def loadPhase2MuonHLTPSets(process, processName = "MYHLT"):

    process.hltPhase2PSetPvClusterComparerForIT = cms.PSet(
        track_chi2_max = cms.double( 20.0 ),
        track_pt_max = cms.double( 100.0 ),
        track_prob_min = cms.double( -1.0 ),
        track_pt_min = cms.double( 1.0 )
    )

    process.HLTSiStripClusterChargeCutNone = cms.PSet(
        value = cms.double(-1.0)
    )

    process.HLTPSetMuonCkfTrajectoryFilter = cms.PSet(
        ComponentType = cms.string('CkfBaseTrajectoryFilter'),
        chargeSignificance = cms.double(-1.0),
        constantValueForLostHitsFractionFilter = cms.double(1.0),
        extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
        maxCCCLostHits = cms.int32(9999),
        maxConsecLostHits = cms.int32(1),
        #maxConsecLostHits = cms.int32(3),#AP
        maxLostHits = cms.int32(1),
        maxLostHitsFraction = cms.double(999.0),
        maxNumberOfHits = cms.int32(-1),
        minGoodStripCharge = cms.PSet(
            refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
        ),
        minHitsMinPt = cms.int32(3),
        minNumberOfHitsForLoopers = cms.int32(13),
        minNumberOfHitsPerLoop = cms.int32(4),
        minPt = cms.double(0.9),
        minimumNumberOfHits = cms.int32(5),
        nSigmaMinPt = cms.double(5.0),
        pixelSeedExtension = cms.bool(False),
        seedExtension = cms.int32(0),
        seedPairPenalty = cms.int32(0),
        strictSeedExtension = cms.bool(False)
    )

    process.HLTPSetMuonCkfTrajectoryBuilder = cms.PSet(
        ComponentType = cms.string('MuonCkfTrajectoryBuilder'),
        MeasurementTrackerName = cms.string('hltESPMeasurementTracker'),
        #TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        alwaysUseInvalidHits = cms.bool(True),
        deltaEta = cms.double(-1.0),
        deltaPhi = cms.double(-1.0),
        estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
        intermediateCleaning = cms.bool(False),
        lostHitPenalty = cms.double(30.0),
        maxCand = cms.int32(5),
        propagatorAlong = cms.string('PropagatorWithMaterial'),
        propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
        propagatorProximity = cms.string('SteppingHelixPropagatorAny'),
        rescaleErrorIfFail = cms.double(1.0),
        trajectoryFilter = cms.PSet(
            refToPSet_ = cms.string('HLTPSetMuonCkfTrajectoryFilter')
        ),
        updator = cms.string('hltESPKFUpdator'),
        useSeedLayer = cms.bool(False),
        seedAs5DHit = cms.bool(False),
    )

    process.HLTIter2Phase2L3MuonPSetGroupedCkfTrajectoryBuilderIT = cms.PSet(
        ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
        MeasurementTrackerName = cms.string('hltIter2HighPtTkMuESPMeasurementTracker'),
        #TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        alwaysUseInvalidHits = cms.bool(False),
        bestHitOnly = cms.bool(True),
        estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
        foundHitBonus = cms.double(1000.0),
        intermediateCleaning = cms.bool(True),
        keepOriginalIfRebuildFails = cms.bool(False),
        lockHits = cms.bool(True),
        lostHitPenalty = cms.double(30.0),
        maxCand = cms.int32(2),
        minNrOfHitsForRebuild = cms.int32(5),
        propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
        propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
        requireSeedHitsInRebuild = cms.bool(False),
        trajectoryFilter = cms.PSet(
            refToPSet_ = cms.string('HLTIter2Phase2L3MuonPSetTrajectoryFilterIT')
        ),
        updator = cms.string('hltESPKFUpdator'),
        useSameTrajFilter = cms.bool(True),
        seedAs5DHit = cms.bool(False),
    )

    process.HLTIter2Phase2L3MuonPSetTrajectoryFilterIT = cms.PSet(
        ComponentType = cms.string('CkfBaseTrajectoryFilter'),
        chargeSignificance = cms.double(-1.0),
        constantValueForLostHitsFractionFilter = cms.double(1.0),
        extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
        maxCCCLostHits = cms.int32(9999),
        maxConsecLostHits = cms.int32(3),
        maxLostHits = cms.int32(1),
        maxLostHitsFraction = cms.double(999.0),
        maxNumberOfHits = cms.int32(100),
        minGoodStripCharge = cms.PSet(
            refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
        ),
        minHitsMinPt = cms.int32(3),
        minNumberOfHitsForLoopers = cms.int32(13),
        minNumberOfHitsPerLoop = cms.int32(4),
        minPt = cms.double(0.3),
        minimumNumberOfHits = cms.int32(5),
        nSigmaMinPt = cms.double(5.0),
        pixelSeedExtension = cms.bool(False),
        seedExtension = cms.int32(0),
        seedPairPenalty = cms.int32(0),
        strictSeedExtension = cms.bool(False)
    )

    process.HLTIter0Phase2L3FromL1TkMuonGroupedCkfTrajectoryFilterIT = cms.PSet(
        ComponentType = cms.string('CkfBaseTrajectoryFilter'),
        chargeSignificance = cms.double(-1.0),
        constantValueForLostHitsFractionFilter = cms.double(10.0),
        extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
        maxCCCLostHits = cms.int32(9999),
        maxConsecLostHits = cms.int32(1),
        maxLostHits = cms.int32(999),
        maxLostHitsFraction = cms.double(0.1),
        maxNumberOfHits = cms.int32(100),
        minGoodStripCharge = cms.PSet(
            refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
        ),
        minHitsMinPt = cms.int32(3),
        minNumberOfHitsForLoopers = cms.int32(13),
        minNumberOfHitsPerLoop = cms.int32(4),
        minPt = cms.double(0.9),
        minimumNumberOfHits = cms.int32(3),
        nSigmaMinPt = cms.double(5.0),
        pixelSeedExtension = cms.bool(False),
        seedExtension = cms.int32(0),
        seedPairPenalty = cms.int32(0),
        strictSeedExtension = cms.bool(False)
    )

    process.HLTIter0Phase2L3FromL1TkMuonPSetGroupedCkfTrajectoryBuilderIT = cms.PSet(
        ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
        MeasurementTrackerName = cms.string(''),
        #TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        alwaysUseInvalidHits = cms.bool(True),
        bestHitOnly = cms.bool(True),
        estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
        foundHitBonus = cms.double(1000.0),
        inOutTrajectoryFilter = cms.PSet(
            refToPSet_ = cms.string('HLTIter0Phase2L3FromL1TkMuonGroupedCkfTrajectoryFilterIT')
        ),
        intermediateCleaning = cms.bool(True),
        keepOriginalIfRebuildFails = cms.bool(True),
        lockHits = cms.bool(True),
        lostHitPenalty = cms.double(1.0),
        maxCand = cms.int32(5),
        minNrOfHitsForRebuild = cms.int32(2),
        propagatorAlong = cms.string('PropagatorWithMaterial'),
        propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
        requireSeedHitsInRebuild = cms.bool(True),
        trajectoryFilter = cms.PSet(
            refToPSet_ = cms.string('HLTIter0Phase2L3FromL1TkMuonGroupedCkfTrajectoryFilterIT')
        ),
        updator = cms.string('hltESPKFUpdator'),
        useSameTrajFilter = cms.bool(True),
        seedAs5DHit = cms.bool(False),
    )

    process.HLTIter2Phase2L3FromL1TkMuonPSetGroupedCkfTrajectoryBuilderIT = cms.PSet(
        ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
        MeasurementTrackerName = cms.string('hltIter2HighPtTkMuESPMeasurementTracker'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        #TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
        alwaysUseInvalidHits = cms.bool(False),
        bestHitOnly = cms.bool(True),
        estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
        foundHitBonus = cms.double(1000.0),
        intermediateCleaning = cms.bool(True),
        keepOriginalIfRebuildFails = cms.bool(False),
        lockHits = cms.bool(True),
        lostHitPenalty = cms.double(30.0),
        maxCand = cms.int32(2),
        minNrOfHitsForRebuild = cms.int32(5),
        propagatorAlong = cms.string('PropagatorWithMaterialParabolicMf'),
        propagatorOpposite = cms.string('PropagatorWithMaterialParabolicMfOpposite'),
        requireSeedHitsInRebuild = cms.bool(False),
        trajectoryFilter = cms.PSet(
            refToPSet_ = cms.string('HLTIter2Phase2L3FromL1TkMuonPSetTrajectoryFilterIT')
        ),
        updator = cms.string('hltESPKFUpdator'),
        useSameTrajFilter = cms.bool(True),
        seedAs5DHit = cms.bool(False),
    )

    process.HLTIter2Phase2L3FromL1TkMuonPSetTrajectoryFilterIT = cms.PSet(
        ComponentType = cms.string('CkfBaseTrajectoryFilter'),
        chargeSignificance = cms.double(-1.0),
        constantValueForLostHitsFractionFilter = cms.double(1.0),
        extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
        maxCCCLostHits = cms.int32(9999),
        maxConsecLostHits = cms.int32(3),
        maxLostHits = cms.int32(1),
        maxLostHitsFraction = cms.double(999.0),
        maxNumberOfHits = cms.int32(100),
        minGoodStripCharge = cms.PSet(
            refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
        ),
        minHitsMinPt = cms.int32(3),
        minNumberOfHitsForLoopers = cms.int32(13),
        minNumberOfHitsPerLoop = cms.int32(4),
        minPt = cms.double(0.3),
        minimumNumberOfHits = cms.int32(5),
        nSigmaMinPt = cms.double(5.0),
        pixelSeedExtension = cms.bool(False),
        seedExtension = cms.int32(0),
        seedPairPenalty = cms.int32(0),
        strictSeedExtension = cms.bool(False)
    )

    process.HLTSeedFromProtoTracks = cms.PSet(
        ComponentName = cms.string('SeedFromConsecutiveHitsCreator'),
        MinOneOverPtError = cms.double(1.0),
        OriginTransverseErrorMultiplier = cms.double(1.0),
        SeedMomentumForBOFF = cms.double(5.0),
        TTRHBuilder = cms.string('hltESPTTRHBuilderPixelOnly'),
        forceKinematicWithRegionDirection = cms.bool(False),
        magneticField = cms.string('ParabolicMf'),
        propagator = cms.string('PropagatorWithMaterialParabolicMf')
    )

    process.hltPhase2SeedFromProtoTracks = cms.PSet(
        ComponentName = cms.string("SeedFromConsecutiveHitsCreator" ),
        MinOneOverPtError = cms.double( 1.0 ),
        OriginTransverseErrorMultiplier = cms.double( 1.0 ),
        SeedMomentumForBOFF = cms.double( 5.0 ),
        TTRHBuilder = cms.string("WithTrackAngle"), 
        forceKinematicWithRegionDirection = cms.bool( False ),
        magneticField = cms.string(""),    
        propagator = cms.string("PropagatorWithMaterial")   
    )

    process.HLTIter0Phase2L3MuonGroupedCkfTrajectoryFilterIT = cms.PSet(
        ComponentType = cms.string('CkfBaseTrajectoryFilter'),
        chargeSignificance = cms.double(-1.0),
        constantValueForLostHitsFractionFilter = cms.double(10.0),
        extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
        maxCCCLostHits = cms.int32(9999),
        maxConsecLostHits = cms.int32(1),
        maxLostHits = cms.int32(999),
        maxLostHitsFraction = cms.double(0.1),
        maxNumberOfHits = cms.int32(100),
        minGoodStripCharge = cms.PSet(
            refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
        ),
        minHitsMinPt = cms.int32(3),
        minNumberOfHitsForLoopers = cms.int32(13),
        minNumberOfHitsPerLoop = cms.int32(4),
        minPt = cms.double(0.9),
        minimumNumberOfHits = cms.int32(3),
        nSigmaMinPt = cms.double(5.0),
        pixelSeedExtension = cms.bool(False),
        seedExtension = cms.int32(0),
        seedPairPenalty = cms.int32(0),
        strictSeedExtension = cms.bool(False)
    )

    process.HLTIter0Phase2L3MuonPSetGroupedCkfTrajectoryBuilderIT = cms.PSet(
        ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
        MeasurementTrackerName = cms.string(''),
        #TTRHBuilder = cms.string('hltESPTTRHBWithTrackAngle'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        alwaysUseInvalidHits = cms.bool(True),
        bestHitOnly = cms.bool(True),
        estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
        foundHitBonus = cms.double(1000.0),
        inOutTrajectoryFilter = cms.PSet(
            refToPSet_ = cms.string('HLTIter0Phase2L3MuonGroupedCkfTrajectoryFilterIT')
        ),
        intermediateCleaning = cms.bool(True),
        keepOriginalIfRebuildFails = cms.bool(True),
        lockHits = cms.bool(True),
        lostHitPenalty = cms.double(1.0),
        maxCand = cms.int32(5),
        minNrOfHitsForRebuild = cms.int32(2),
        propagatorAlong = cms.string('PropagatorWithMaterial'),
        propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
        requireSeedHitsInRebuild = cms.bool(True),
        trajectoryFilter = cms.PSet(
            refToPSet_ = cms.string('HLTIter0Phase2L3MuonGroupedCkfTrajectoryFilterIT')
        ),
        updator = cms.string('hltESPKFUpdator'),
        useSameTrajFilter = cms.bool(True),
        seedAs5DHit = cms.bool(False),
    )

    process.HLTIter0Phase2L3MuonPSetGroupedCkfTrajectoryBuilderIT = cms.PSet(
        useSameTrajFilter = cms.bool(True),    
        ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
        MeasurementTrackerName = cms.string(''),
        keepOriginalIfRebuildFails = cms.bool(True),
        lockHits = cms.bool(True),
        lostHitPenalty = cms.double(1.0),
        requireSeedHitsInRebuild = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'),
        propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
        
        trajectoryFilter = cms.PSet(
            refToPSet_ = cms.string('HLTIter0Phase2L3MuonGroupedCkfTrajectoryFilterIT')
        ),
        propagatorAlong = cms.string('PropagatorWithMaterial'),
        minNrOfHitsForRebuild = cms.int32(2),
        maxCand = cms.int32(5),
        alwaysUseInvalidHits = cms.bool(True),
        estimator = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
        foundHitBonus = cms.double(1000.0),
        inOutTrajectoryFilter = cms.PSet(
            refToPSet_ = cms.string('HLTIter0Phase2L3MuonGroupedCkfTrajectoryFilterIT')
        ),
        intermediateCleaning = cms.bool(True),
        updator = cms.string('hltESPKFUpdator'),
        bestHitOnly = cms.bool(True),
        seedAs5DHit = cms.bool(False)
    )

    process.HLTIter0GroupedCkfTrajectoryBuilderIT = cms.PSet( 
        keepOriginalIfRebuildFails = cms.bool( False ),
        lockHits = cms.bool( True ),
        maxDPhiForLooperReconstruction = cms.double( 2.0 ),
        propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
        trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter0PSetTrajectoryFilterIT" ) ),
        doSeedingRegionRebuilding = cms.bool( False ),
        useHitsSplitting = cms.bool( False ),
        maxCand = cms.int32( 2 ),
        estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator9" ),
        intermediateCleaning = cms.bool( True ),
        bestHitOnly = cms.bool( True ),
        useSameTrajFilter = cms.bool( True ),
        MeasurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
        ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
        lostHitPenalty = cms.double( 30.0 ),
        requireSeedHitsInRebuild = cms.bool( True ),
        TTRHBuilder = cms.string( "WithTrackAngle" ),
        maxPtForLooperReconstruction = cms.double( 0.7 ),
        cleanTrajectoryAfterInOut = cms.bool( False ),
        propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
        minNrOfHitsForRebuild = cms.int32( 5 ),
        alwaysUseInvalidHits = cms.bool( False ),
        inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter0PSetTrajectoryFilterIT" ) ),
        foundHitBonus = cms.double( 5.0 ),
        updator = cms.string( "hltESPKFUpdator" ),
        seedAs5DHit = cms.bool(False)
    )

    process.HLTIter0PSetTrajectoryFilterIT = cms.PSet( 
        minimumNumberOfHits = cms.int32( 4 ),
        ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
        seedExtension = cms.int32( 0 ),
        chargeSignificance = cms.double( -1.0 ),
        pixelSeedExtension = cms.bool( False ),
        strictSeedExtension = cms.bool( False ),
        nSigmaMinPt = cms.double( 5.0 ),
        maxCCCLostHits = cms.int32( 0 ),
        minPt = cms.double( 0.3 ),
        maxConsecLostHits = cms.int32( 1 ),
        extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
        constantValueForLostHitsFractionFilter = cms.double( 1.0 ),
        seedPairPenalty = cms.int32( 0 ),
        maxNumberOfHits = cms.int32( 100 ),
        minNumberOfHitsForLoopers = cms.int32( 13 ),
        minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
        minNumberOfHitsPerLoop = cms.int32( 4 ),
        minHitsMinPt = cms.int32( 4 ),
        maxLostHitsFraction = cms.double( 999.0 ),
        maxLostHits = cms.int32( 1 )
    )

    process.HLTSiStripClusterChargeCutLoose = cms.PSet(
        value = cms.double( 1620.0 )
    )

    process.HLTIter1GroupedCkfTrajectoryBuilderIT = cms.PSet( 
        useSameTrajFilter = cms.bool( True ),
        ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
        MeasurementTrackerName = cms.string( "hltIter1ESPMeasurementTracker" ),
        keepOriginalIfRebuildFails = cms.bool( False ),
        lostHitPenalty = cms.double( 30.0 ),
        lockHits = cms.bool( True ),
        requireSeedHitsInRebuild = cms.bool( True ),
        TTRHBuilder = cms.string( "WithTrackAngle" ),
        propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
        trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter1PSetTrajectoryFilterIT" ) ),
        propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
        minNrOfHitsForRebuild = cms.int32( 5 ),
        maxCand = cms.int32( 2 ),
        alwaysUseInvalidHits = cms.bool( False ),
        estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator16" ),
        intermediateCleaning = cms.bool( True ),
        foundHitBonus = cms.double( 5.0 ),
        updator = cms.string( "hltESPKFUpdator" ),
        bestHitOnly = cms.bool( True ),
        seedAs5DHit = cms.bool(False)
    )

    process.HLTIter1PSetTrajectoryFilterIT = cms.PSet( 
        minimumNumberOfHits = cms.int32( 3 ),
        ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
        seedExtension = cms.int32( 0 ),
        chargeSignificance = cms.double( -1.0 ),
        pixelSeedExtension = cms.bool( False ),
        strictSeedExtension = cms.bool( False ),
        nSigmaMinPt = cms.double( 5.0 ),
        maxCCCLostHits = cms.int32( 0 ),
        minPt = cms.double( 0.2 ),
        maxConsecLostHits = cms.int32( 1 ),
        extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
        constantValueForLostHitsFractionFilter = cms.double( 1.0 ),
        seedPairPenalty = cms.int32( 0 ),
        maxNumberOfHits = cms.int32( 100 ),
        minNumberOfHitsForLoopers = cms.int32( 13 ),
        minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
        minNumberOfHitsPerLoop = cms.int32( 4 ),
        minHitsMinPt = cms.int32( 3 ),
        maxLostHitsFraction = cms.double( 999.0 ),
        maxLostHits = cms.int32( 1 )
    )

    process.HLTIter2GroupedCkfTrajectoryBuilderIT = cms.PSet( 
        keepOriginalIfRebuildFails = cms.bool( False ),
        lockHits = cms.bool( True ),
        maxDPhiForLooperReconstruction = cms.double( 2.0 ),
        propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
        trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter2PSetTrajectoryFilterIT" ) ),
        doSeedingRegionRebuilding = cms.bool( False ),
        useHitsSplitting = cms.bool( False ),
        maxCand = cms.int32( 2 ),
        estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator16" ),
        intermediateCleaning = cms.bool( True ),
        bestHitOnly = cms.bool( True ),
        useSameTrajFilter = cms.bool( True ),
        MeasurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
        ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
        lostHitPenalty = cms.double( 30.0 ),
        requireSeedHitsInRebuild = cms.bool( True ),
        TTRHBuilder = cms.string( "WithTrackAngle" ),
        maxPtForLooperReconstruction = cms.double( 0.7 ),
        cleanTrajectoryAfterInOut = cms.bool( False ),
        propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
        minNrOfHitsForRebuild = cms.int32( 5 ),
        alwaysUseInvalidHits = cms.bool( False ),
        inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter2PSetTrajectoryFilterIT" ) ),
        foundHitBonus = cms.double( 5.0 ),
        updator = cms.string( "hltESPKFUpdator" ),
        seedAs5DHit = cms.bool(False)
    )

    process.HLTIter2PSetTrajectoryFilterIT = cms.PSet( 
        minimumNumberOfHits = cms.int32( 3 ),
        ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
        seedExtension = cms.int32( 1 ),
        chargeSignificance = cms.double( -1.0 ),
        pixelSeedExtension = cms.bool( False ),
        strictSeedExtension = cms.bool( False ),
        nSigmaMinPt = cms.double( 5.0 ),
        maxCCCLostHits = cms.int32( 0 ),
        minPt = cms.double( 0.3 ),
        maxConsecLostHits = cms.int32( 1 ),
        extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
        constantValueForLostHitsFractionFilter = cms.double( 1.0 ),
        seedPairPenalty = cms.int32( 0 ),
        maxNumberOfHits = cms.int32( 100 ),
        minNumberOfHitsForLoopers = cms.int32( 13 ),
        minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
        minNumberOfHitsPerLoop = cms.int32( 4 ),
        minHitsMinPt = cms.int32( 3 ),
        maxLostHitsFraction = cms.double( 999.0 ),
        maxLostHits = cms.int32( 1 )
    )

    return process

def loadPhase2MuonHLTESProducers(process, processName = "MYHLT"):

    process.hltESPFastSteppingHelixPropagatorAny = cms.ESProducer("SteppingHelixPropagatorESProducer",
        ApplyRadX0Correction = cms.bool(True),
        AssumeNoMaterial = cms.bool(False),
        ComponentName = cms.string('hltESPFastSteppingHelixPropagatorAny'),
        NoErrorPropagation = cms.bool(False),
        PropagationDirection = cms.string('anyDirection'),
        SetVBFPointer = cms.bool(False),
        VBFName = cms.string('VolumeBasedMagneticField'),
        debug = cms.bool(False),
        endcapShiftInZNeg = cms.double(0.0),
        endcapShiftInZPos = cms.double(0.0),
        returnTangentPlane = cms.bool(True),
        sendLogWarning = cms.bool(False),
        useEndcapShiftsInZ = cms.bool(False),
        useInTeslaFromMagField = cms.bool(False),
        useIsYokeFlag = cms.bool(True),
        useMagVolumes = cms.bool(True),
        useMatVolumes = cms.bool(True),
        useTuningForL2Speed = cms.bool(True)
    )

    process.hltESPFastSteppingHelixPropagatorOpposite = cms.ESProducer("SteppingHelixPropagatorESProducer",
        ApplyRadX0Correction = cms.bool(True),
        AssumeNoMaterial = cms.bool(False),
        ComponentName = cms.string('hltESPFastSteppingHelixPropagatorOpposite'),
        NoErrorPropagation = cms.bool(False),
        PropagationDirection = cms.string('oppositeToMomentum'),
        SetVBFPointer = cms.bool(False),
        VBFName = cms.string('VolumeBasedMagneticField'),
        debug = cms.bool(False),
        endcapShiftInZNeg = cms.double(0.0),
        endcapShiftInZPos = cms.double(0.0),
        returnTangentPlane = cms.bool(True),
        sendLogWarning = cms.bool(False),
        useEndcapShiftsInZ = cms.bool(False),
        useInTeslaFromMagField = cms.bool(False),
        useIsYokeFlag = cms.bool(True),
        useMagVolumes = cms.bool(True),
        useMatVolumes = cms.bool(True),
        useTuningForL2Speed = cms.bool(True)
    )

    process.hltESPMuonTransientTrackingRecHitBuilder = cms.ESProducer("MuonTransientTrackingRecHitBuilderESProducer",
        ComponentName = cms.string('hltESPMuonTransientTrackingRecHitBuilder')
    )

    process.hltESPKFUpdator = cms.ESProducer("KFUpdatorESProducer",
        ComponentName = cms.string('hltESPKFUpdator')
    )

    process.hltESPDummyDetLayerGeometry = cms.ESProducer("DetLayerGeometryESProducer",
        ComponentName = cms.string('hltESPDummyDetLayerGeometry')
    )

    process.hltESPChi2MeasurementEstimator30 = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
        ComponentName = cms.string('hltESPChi2MeasurementEstimator30'),
        MaxChi2 = cms.double(30.0),
        MaxDisplacement = cms.double(100.0),
        MaxSagitta = cms.double(-1.0),
        MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
        MinimalTolerance = cms.double(10.0),
        appendToDataLabel = cms.string(''),
        nSigma = cms.double(3.0)
    )

    process.hltESPKFTrajectorySmootherForMuonTrackLoader = cms.ESProducer("KFTrajectorySmootherESProducer",
        ComponentName = cms.string('hltESPKFTrajectorySmootherForMuonTrackLoader'),
        Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
        Propagator = cms.string('hltESPSmartPropagatorAnyOpposite'),
        RecoGeometry = cms.string('hltESPDummyDetLayerGeometry'),
        Updator = cms.string('hltESPKFUpdator'),
        appendToDataLabel = cms.string(''),
        errorRescaling = cms.double(10.0),
        minHits = cms.int32(3)
    )

    process.hltESPSteppingHelixPropagatorOpposite = cms.ESProducer("SteppingHelixPropagatorESProducer",
        ApplyRadX0Correction = cms.bool(True),
        AssumeNoMaterial = cms.bool(False),
        ComponentName = cms.string('hltESPSteppingHelixPropagatorOpposite'),
        NoErrorPropagation = cms.bool(False),
        PropagationDirection = cms.string('oppositeToMomentum'),
        SetVBFPointer = cms.bool(False),
        VBFName = cms.string('VolumeBasedMagneticField'),
        debug = cms.bool(False),
        endcapShiftInZNeg = cms.double(0.0),
        endcapShiftInZPos = cms.double(0.0),
        returnTangentPlane = cms.bool(True),
        sendLogWarning = cms.bool(False),
        useEndcapShiftsInZ = cms.bool(False),
        useInTeslaFromMagField = cms.bool(False),
        useIsYokeFlag = cms.bool(True),
        useMagVolumes = cms.bool(True),
        useMatVolumes = cms.bool(True),
        useTuningForL2Speed = cms.bool(False)
    )

    process.hltESPChi2MeasurementEstimator100 = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
        ComponentName = cms.string('hltESPChi2MeasurementEstimator100'),
        MaxChi2 = cms.double(40.0),
        MaxDisplacement = cms.double(0.5),
        MaxSagitta = cms.double(2.0),
        MinPtForHitRecoveryInGluedDet = cms.double(1e+12),
        MinimalTolerance = cms.double(0.5),
        appendToDataLabel = cms.string(''),
        nSigma = cms.double(4.0)
    )

    process.hltESPChi2ChargeMeasurementEstimator30 = cms.ESProducer("Chi2ChargeMeasurementEstimatorESProducer",
        ComponentName = cms.string('hltESPChi2ChargeMeasurementEstimator30'),
        MaxChi2 = cms.double(30.0),
        MaxDisplacement = cms.double(100.0),
        MaxSagitta = cms.double(-1.0),
        MinPtForHitRecoveryInGluedDet = cms.double(1000000.0),
        MinimalTolerance = cms.double(10.0),
        appendToDataLabel = cms.string(''),
        clusterChargeCut = cms.PSet(
            refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone')
        ),
        nSigma = cms.double(3.0),
        pTChargeCutThreshold = cms.double(-1.0)
    )

    process.hltESPMeasurementTracker = cms.ESProducer("MeasurementTrackerESProducer",
        ComponentName = cms.string('hltESPMeasurementTracker'),
        DebugPixelModuleQualityDB = cms.untracked.bool(False),
        DebugPixelROCQualityDB = cms.untracked.bool(False),
        DebugStripAPVFiberQualityDB = cms.untracked.bool(False),
        DebugStripModuleQualityDB = cms.untracked.bool(False),
        DebugStripStripQualityDB = cms.untracked.bool(False),
        HitMatcher = cms.string('StandardMatcher'),
        MaskBadAPVFibers = cms.bool(True),
        #PixelCPE = cms.string('hltESPPixelCPEGeneric'),
        PixelCPE = cms.string('PixelCPEGeneric'),
        SiStripQualityLabel = cms.string(''),
        StripCPE = cms.string('hltESPStripCPEfromTrackAngle'),
        UsePixelModuleQualityDB = cms.bool(True),
        UsePixelROCQualityDB = cms.bool(True),
        UseStripAPVFiberQualityDB = cms.bool(True),
        UseStripModuleQualityDB = cms.bool(True),
        UseStripStripQualityDB = cms.bool(True),
        badStripCuts = cms.PSet(
            TEC = cms.PSet(
                maxBad = cms.uint32(4),
                maxConsecutiveBad = cms.uint32(2)
            ),
            TIB = cms.PSet(
                maxBad = cms.uint32(4),
                maxConsecutiveBad = cms.uint32(2)
            ),
            TID = cms.PSet(
                maxBad = cms.uint32(4),
                maxConsecutiveBad = cms.uint32(2)
            ),
            TOB = cms.PSet(
                maxBad = cms.uint32(4),
                maxConsecutiveBad = cms.uint32(2)
            )
        )
    )

    process.hltESPRKTrajectorySmoother = cms.ESProducer("KFTrajectorySmootherESProducer",
        ComponentName = cms.string('hltESPRKTrajectorySmoother'),
        Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
        Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
        RecoGeometry = cms.string('hltESPGlobalDetLayerGeometry'),
        Updator = cms.string('hltESPKFUpdator'),
        appendToDataLabel = cms.string(''),
        errorRescaling = cms.double(100.0),
        minHits = cms.int32(3)
    )

    process.hltESPGlobalDetLayerGeometry = cms.ESProducer("GlobalDetLayerGeometryESProducer",
        ComponentName = cms.string('hltESPGlobalDetLayerGeometry')
    )

    process.hltESPRKTrajectoryFitter = cms.ESProducer("KFTrajectoryFitterESProducer",
        ComponentName = cms.string('hltESPRKTrajectoryFitter'),
        Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
        Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
        RecoGeometry = cms.string('hltESPGlobalDetLayerGeometry'),
        Updator = cms.string('hltESPKFUpdator'),
        appendToDataLabel = cms.string(''),
        minHits = cms.int32(3)
    )

    process.hltESPRungeKuttaTrackerPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
        ComponentName = cms.string('hltESPRungeKuttaTrackerPropagator'),
        Mass = cms.double(0.105),
        MaxDPhi = cms.double(1.6),
        PropagationDirection = cms.string('alongMomentum'),
        SimpleMagneticField = cms.string(''),
        ptMin = cms.double(-1.0),
        useRungeKutta = cms.bool(True)
    )

    process.hltESPKFFittingSmootherWithOutliersRejectionAndRK = cms.ESProducer("KFFittingSmootherESProducer",
        BreakTrajWith2ConsecutiveMissing = cms.bool(True),
        ComponentName = cms.string('hltESPKFFittingSmootherWithOutliersRejectionAndRK'),
        EstimateCut = cms.double(20.0),
        Fitter = cms.string('hltESPRKTrajectoryFitter'),
        LogPixelProbabilityCut = cms.double(-14.0),
        MaxFractionOutliers = cms.double(0.3),
        MaxNumberOfOutliers = cms.int32(3),
        MinDof = cms.int32(2),
        MinNumberOfHits = cms.int32(3),
        NoInvalidHitsBeginEnd = cms.bool(True),
        NoOutliersBeginEnd = cms.bool(False),
        RejectTracks = cms.bool(True),
        Smoother = cms.string('hltESPRKTrajectorySmoother'),
        appendToDataLabel = cms.string('')
    )

    process.hltESPSmartPropagatorAnyOpposite = cms.ESProducer("SmartPropagatorESProducer",
        ComponentName = cms.string('hltESPSmartPropagatorAnyOpposite'),
        Epsilon = cms.double(5.0),
        MuonPropagator = cms.string('SteppingHelixPropagatorAny'),
        PropagationDirection = cms.string('oppositeToMomentum'),
        TrackerPropagator = cms.string('PropagatorWithMaterialOpposite')
    )

    process.hltESPSmartPropagatorAny = cms.ESProducer("SmartPropagatorESProducer",
        ComponentName = cms.string('hltESPSmartPropagatorAny'),
        Epsilon = cms.double(5.0),
        MuonPropagator = cms.string('SteppingHelixPropagatorAny'),
        PropagationDirection = cms.string('alongMomentum'),
        TrackerPropagator = cms.string('PropagatorWithMaterial')
    )

    process.hltESPSmartPropagator = cms.ESProducer("SmartPropagatorESProducer",
        ComponentName = cms.string('hltESPSmartPropagator'),
        Epsilon = cms.double(5.0),
        MuonPropagator = cms.string('hltESPSteppingHelixPropagatorAlong'),
        PropagationDirection = cms.string('alongMomentum'),
        TrackerPropagator = cms.string('PropagatorWithMaterial')
    )

    process.hltESPL3MuKFTrajectoryFitter = cms.ESProducer("KFTrajectoryFitterESProducer",
        ComponentName = cms.string('hltESPL3MuKFTrajectoryFitter'),
        Estimator = cms.string('hltESPChi2MeasurementEstimator30'),
        Propagator = cms.string('hltESPSmartPropagatorAny'),
        RecoGeometry = cms.string('hltESPDummyDetLayerGeometry'),
        Updator = cms.string('hltESPKFUpdator'),
        appendToDataLabel = cms.string(''),
        minHits = cms.int32(3)
    )

    process.hltESPSteppingHelixPropagatorAlong = cms.ESProducer("SteppingHelixPropagatorESProducer",
        ApplyRadX0Correction = cms.bool(True),
        AssumeNoMaterial = cms.bool(False),
        ComponentName = cms.string('hltESPSteppingHelixPropagatorAlong'),
        NoErrorPropagation = cms.bool(False),
        PropagationDirection = cms.string('alongMomentum'),
        SetVBFPointer = cms.bool(False),
        VBFName = cms.string('VolumeBasedMagneticField'),
        debug = cms.bool(False),
        endcapShiftInZNeg = cms.double(0.0),
        endcapShiftInZPos = cms.double(0.0),
        returnTangentPlane = cms.bool(True),
        sendLogWarning = cms.bool(False),
        useEndcapShiftsInZ = cms.bool(False),
        useInTeslaFromMagField = cms.bool(False),
        useIsYokeFlag = cms.bool(True),
        useMagVolumes = cms.bool(True),
        useMatVolumes = cms.bool(True),
        useTuningForL2Speed = cms.bool(False)
    )

    process.hltESPTTRHBuilderPixelOnly = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
        ComponentName = cms.string('hltESPTTRHBuilderPixelOnly'),
        ComputeCoarseLocalPositionFromDisk = cms.bool(False),
        Matcher = cms.string('StandardMatcher'),
        #Phase2StripCPE = cms.string('Phase2StripCPE'),
        Phase2StripCPE = cms.string(''),
        PixelCPE = cms.string('PixelCPEGeneric'),
        StripCPE = cms.string('Fake')
    )

    process.hltPixelTracksCleanerBySharedHits = cms.ESProducer("PixelTrackCleanerBySharedHitsESProducer",
        ComponentName = cms.string('hltPixelTracksCleanerBySharedHits'),
        appendToDataLabel = cms.string(''),
        useQuadrupletAlgo = cms.bool(False)
    )

    process.hltESPTrackAlgoPriorityOrder = cms.ESProducer("TrackAlgoPriorityOrderESProducer",
        ComponentName = cms.string('hltESPTrackAlgoPriorityOrder'),
        algoOrder = cms.vstring(),
        appendToDataLabel = cms.string('')
    )

    process.hltESPTrajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
        ComponentName = cms.string('hltESPTrajectoryCleanerBySharedHits'),
        ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
        MissingHitPenalty = cms.double(0.0),
        ValidHitBonus = cms.double(100.0),
        allowSharedFirstHit = cms.bool(False),
        fractionShared = cms.double(0.5)
    )

    process.hltESPChi2ChargeMeasurementEstimator9 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
        appendToDataLabel = cms.string( "" ),
        clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
        MinimalTolerance = cms.double( 0.5 ),
        MaxDisplacement = cms.double( 0.5 ),
        ComponentName = cms.string( "hltESPChi2ChargeMeasurementEstimator9" ),
        pTChargeCutThreshold = cms.double( 15.0 ),
        nSigma = cms.double( 3.0 ),
        MaxSagitta = cms.double( 2.0 ),
        MaxChi2 = cms.double( 9.0 ),
        MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 )
    )

    process.hltESPChi2ChargeMeasurementEstimator16 = cms.ESProducer( "Chi2ChargeMeasurementEstimatorESProducer",
        appendToDataLabel = cms.string( "" ),
        clusterChargeCut = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutLoose" ) ),
        MinimalTolerance = cms.double( 0.5 ),
        MaxDisplacement = cms.double( 0.5 ),
        ComponentName = cms.string( "hltESPChi2ChargeMeasurementEstimator16" ),
        pTChargeCutThreshold = cms.double( -1.0 ),
        nSigma = cms.double( 3.0 ),
        MaxSagitta = cms.double( 2.0 ),
        MaxChi2 = cms.double( 16.0 ),
        MinPtForHitRecoveryInGluedDet = cms.double( 1000000.0 )
    )

    return process

def customizeMuonLocalReco(process, processName = "MYHLT"):
    # -- using offline muon local reconstruction
    from Configuration.StandardSequences.Reconstruction_cff import \
    dt1DRecHits, dt4DSegments, csc2DRecHits, cscSegments, \
    rpcRecHits, gemRecHits, gemSegments, me0RecHits, me0Segments

    process.hltDt1DRecHits = dt1DRecHits.clone(
        dtDigiLabel = cms.InputTag("simMuonDTDigis")
    )

    process.hltDt4DSegments = dt4DSegments.clone(
        recHits1DLabel = cms.InputTag("hltDt1DRecHits")
    )

    process.hltCsc2DRecHits = csc2DRecHits.clone(
        stripDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigi"),
        wireDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")
    )

    process.hltCscSegments = cscSegments.clone(
        inputObjects = cms.InputTag("hltCsc2DRecHits")
    )

    process.hltRpcRecHits = rpcRecHits.clone(
        rpcDigiLabel = cms.InputTag("simMuonRPCDigis")
    )

    process.hltGemRecHits = gemRecHits.clone(
        gemDigiLabel = cms.InputTag("simMuonGEMDigis")
    )

    process.hltGemSegments = gemSegments.clone(
        gemRecHitLabel = cms.InputTag("hltGemRecHits")
    )

    process.hltMe0RecHits = me0RecHits.clone(
        me0DigiLabel = cms.InputTag("simMuonME0PseudoReDigis")
    )

    process.hltMe0Segments = me0Segments.clone(
        me0RecHitLabel = cms.InputTag("hltMe0RecHits")
    )

    process.HLTMuonGemLocalRecoSequence = cms.Sequence(
        process.hltGemRecHits+
        process.hltGemSegments+
        process.hltMe0RecHits+
        process.hltMe0Segments
    )

    process.HLTMuonLocalRecoSequence = cms.Sequence(
        process.hltDt1DRecHits+
        process.hltDt4DSegments+
        process.hltCsc2DRecHits+
        process.hltCscSegments+
        process.hltRpcRecHits+
        process.HLTMuonGemLocalRecoSequence
    )

    return process

def customizeTrackerLocalReco(process, processName = "MYHLT"):

    # process.siPixelClusters = process.siPixelClusters.clone(
    #     src = cms.InputTag("simSiPixelDigis","Pixel")
    # )

    # process.hltSiPixelClustersCache = cms.EDProducer("SiPixelClusterShapeCacheProducer",
    #     onDemand = cms.bool(False),
    #     src = cms.InputTag("siPixelClusters")
    # )

    # process.siPixelRecHits = process.siPixelRecHits.clone(
    #     src = cms.InputTag("siPixelClusters")
    # )

    # process.hltSiStripClusters = cms.EDProducer("MeasurementTrackerEventProducer",
    #     Phase2TrackerCluster1DProducer = cms.string('siPhase2Clusters'),
    #     badPixelFEDChannelCollectionLabels = cms.VInputTag("siPixelDigis"),
    #     inactivePixelDetectorLabels = cms.VInputTag(),
    #     inactiveStripDetectorLabels = cms.VInputTag("siStripDigis"),
    #     measurementTracker = cms.string(''),
    #     mightGet = cms.optional.untracked.vstring,
    #     pixelCablingMapLabel = cms.string(''),
    #     pixelClusterProducer = cms.string('siPixelClusters'),
    #     skipClusters = cms.InputTag(""),
    #     stripClusterProducer = cms.string(''),
    #     switchOffPixelsIfEmpty = cms.bool(True)
    # )

    process.HLTDoLocalPixelSequence = cms.Sequence(
        process.siPixelClusters+
        process.siPixelClusterShapeCache+
        process.siPixelRecHits
    )

    process.HLTDoLocalStripSequence = cms.Sequence(
        process.siPhase2Clusters+
        process.siPhase2VectorHits+
        process.MeasurementTrackerEvent
    )

    return process

def customizeL2MuonReco(process, processName = "MYHLT"):

    process.hltL2OfflineMuonSeeds = cms.EDProducer("MuonSeedGenerator",
        CSCRecSegmentLabel = cms.InputTag("hltCscSegments"),
        CSC_01 = cms.vdouble(
            0.166, 0.0, 0.0, 0.031, 0.0, 
            0.0
        ),
        CSC_01_1_scale = cms.vdouble(-1.915329, 0.0),
        CSC_02 = cms.vdouble(
            0.612, -0.207, 0.0, 0.067, -0.001, 
            0.0
        ),
        CSC_03 = cms.vdouble(
            0.787, -0.338, 0.029, 0.101, -0.008, 
            0.0
        ),
        CSC_12 = cms.vdouble(
            -0.161, 0.254, -0.047, 0.042, -0.007, 
            0.0
        ),
        CSC_12_1_scale = cms.vdouble(-6.434242, 0.0),
        CSC_12_2_scale = cms.vdouble(-1.63622, 0.0),
        CSC_12_3_scale = cms.vdouble(-1.63622, 0.0),
        CSC_13 = cms.vdouble(
            0.901, -1.302, 0.533, 0.045, 0.005, 
            0.0
        ),
        CSC_13_2_scale = cms.vdouble(-6.077936, 0.0),
        CSC_13_3_scale = cms.vdouble(-1.701268, 0.0),
        CSC_14 = cms.vdouble(
            0.606, -0.181, -0.002, 0.111, -0.003, 
            0.0
        ),
        CSC_14_3_scale = cms.vdouble(-1.969563, 0.0),
        CSC_23 = cms.vdouble(
            -0.081, 0.113, -0.029, 0.015, 0.008, 
            0.0
        ),
        CSC_23_1_scale = cms.vdouble(-19.084285, 0.0),
        CSC_23_2_scale = cms.vdouble(-6.079917, 0.0),
        CSC_24 = cms.vdouble(
            0.004, 0.021, -0.002, 0.053, 0.0, 
            0.0
        ),
        CSC_24_1_scale = cms.vdouble(-6.055701, 0.0),
        CSC_34 = cms.vdouble(
            0.062, -0.067, 0.019, 0.021, 0.003, 
            0.0
        ),
        CSC_34_1_scale = cms.vdouble(-11.520507, 0.0),
        DTRecSegmentLabel = cms.InputTag("hltDt4DSegments"),
        DT_12 = cms.vdouble(
            0.183, 0.054, -0.087, 0.028, 0.002, 
            0.0
        ),
        DT_12_1_scale = cms.vdouble(-3.692398, 0.0),
        DT_12_2_scale = cms.vdouble(-3.518165, 0.0),
        DT_13 = cms.vdouble(
            0.315, 0.068, -0.127, 0.051, -0.002, 
            0.0
        ),
        DT_13_1_scale = cms.vdouble(-4.520923, 0.0),
        DT_13_2_scale = cms.vdouble(-4.257687, 0.0),
        DT_14 = cms.vdouble(
            0.359, 0.052, -0.107, 0.072, -0.004, 
            0.0
        ),
        DT_14_1_scale = cms.vdouble(-5.644816, 0.0),
        DT_14_2_scale = cms.vdouble(-4.808546, 0.0),
        DT_23 = cms.vdouble(
            0.13, 0.023, -0.057, 0.028, 0.004, 
            0.0
        ),
        DT_23_1_scale = cms.vdouble(-5.320346, 0.0),
        DT_23_2_scale = cms.vdouble(-5.117625, 0.0),
        DT_24 = cms.vdouble(
            0.176, 0.014, -0.051, 0.051, 0.003, 
            0.0
        ),
        DT_24_1_scale = cms.vdouble(-7.490909, 0.0),
        DT_24_2_scale = cms.vdouble(-6.63094, 0.0),
        DT_34 = cms.vdouble(
            0.044, 0.004, -0.013, 0.029, 0.003, 
            0.0
        ),
        DT_34_1_scale = cms.vdouble(-13.783765, 0.0),
        DT_34_2_scale = cms.vdouble(-11.901897, 0.0),
        EnableCSCMeasurement = cms.bool(True),
        EnableDTMeasurement = cms.bool(True),
        EnableME0Measurement = cms.bool(False),
        ME0RecSegmentLabel = cms.InputTag("hltMe0Segments"),
        OL_1213 = cms.vdouble(
            0.96, -0.737, 0.0, 0.052, 0.0, 
            0.0
        ),
        OL_1213_0_scale = cms.vdouble(-4.488158, 0.0),
        OL_1222 = cms.vdouble(
            0.848, -0.591, 0.0, 0.062, 0.0, 
            0.0
        ),
        OL_1222_0_scale = cms.vdouble(-5.810449, 0.0),
        OL_1232 = cms.vdouble(
            0.184, 0.0, 0.0, 0.066, 0.0, 
            0.0
        ),
        OL_1232_0_scale = cms.vdouble(-5.964634, 0.0),
        OL_2213 = cms.vdouble(
            0.117, 0.0, 0.0, 0.044, 0.0, 
            0.0
        ),
        OL_2213_0_scale = cms.vdouble(-7.239789, 0.0),
        OL_2222 = cms.vdouble(
            0.107, 0.0, 0.0, 0.04, 0.0, 
            0.0
        ),
        OL_2222_0_scale = cms.vdouble(-7.667231, 0.0),
        SMB_10 = cms.vdouble(
            1.387, -0.038, 0.0, 0.19, 0.0, 
            0.0
        ),
        SMB_10_0_scale = cms.vdouble(2.448566, 0.0),
        SMB_11 = cms.vdouble(
            1.247, 0.72, -0.802, 0.229, -0.075, 
            0.0
        ),
        SMB_11_0_scale = cms.vdouble(2.56363, 0.0),
        SMB_12 = cms.vdouble(
            2.128, -0.956, 0.0, 0.199, 0.0, 
            0.0
        ),
        SMB_12_0_scale = cms.vdouble(2.283221, 0.0),
        SMB_20 = cms.vdouble(
            1.011, -0.052, 0.0, 0.188, 0.0, 
            0.0
        ),
        SMB_20_0_scale = cms.vdouble(1.486168, 0.0),
        SMB_21 = cms.vdouble(
            1.043, -0.124, 0.0, 0.183, 0.0, 
            0.0
        ),
        SMB_21_0_scale = cms.vdouble(1.58384, 0.0),
        SMB_22 = cms.vdouble(
            1.474, -0.758, 0.0, 0.185, 0.0, 
            0.0
        ),
        SMB_22_0_scale = cms.vdouble(1.346681, 0.0),
        SMB_30 = cms.vdouble(
            0.505, -0.022, 0.0, 0.215, 0.0, 
            0.0
        ),
        SMB_30_0_scale = cms.vdouble(-3.629838, 0.0),
        SMB_31 = cms.vdouble(
            0.549, -0.145, 0.0, 0.207, 0.0, 
            0.0
        ),
        SMB_31_0_scale = cms.vdouble(-3.323768, 0.0),
        SMB_32 = cms.vdouble(
            0.67, -0.327, 0.0, 0.22, 0.0, 
            0.0
        ),
        SMB_32_0_scale = cms.vdouble(-3.054156, 0.0),
        SME_11 = cms.vdouble(
            3.295, -1.527, 0.112, 0.378, 0.02, 
            0.0
        ),
        SME_11_0_scale = cms.vdouble(1.325085, 0.0),
        SME_12 = cms.vdouble(
            0.102, 0.599, 0.0, 0.38, 0.0, 
            0.0
        ),
        SME_12_0_scale = cms.vdouble(2.279181, 0.0),
        SME_13 = cms.vdouble(
            -1.286, 1.711, 0.0, 0.356, 0.0, 
            0.0
        ),
        SME_13_0_scale = cms.vdouble(0.104905, 0.0),
        SME_21 = cms.vdouble(
            -0.529, 1.194, -0.358, 0.472, 0.086, 
            0.0
        ),
        SME_21_0_scale = cms.vdouble(-0.040862, 0.0),
        SME_22 = cms.vdouble(
            -1.207, 1.491, -0.251, 0.189, 0.243, 
            0.0
        ),
        SME_22_0_scale = cms.vdouble(-3.457901, 0.0),
        SME_31 = cms.vdouble(
            -1.594, 1.482, -0.317, 0.487, 0.097, 
            0.0
        ),
        SME_32 = cms.vdouble(
            -0.901, 1.333, -0.47, 0.41, 0.073, 
            0.0
        ),
        SME_41 = cms.vdouble(
            -0.003, 0.005, 0.005, 0.608, 0.076, 
            0.0
        ),
        SME_42 = cms.vdouble(
            -0.003, 0.005, 0.005, 0.608, 0.076, 
            0.0
        ),
        beamSpotTag = cms.InputTag("offlineBeamSpot"),
        crackEtas = cms.vdouble(0.2, 1.6, 1.7),
        crackWindow = cms.double(0.04),
        deltaEtaCrackSearchWindow = cms.double(0.25),
        deltaEtaSearchWindow = cms.double(0.2),
        deltaPhiSearchWindow = cms.double(0.25),
        scaleDT = cms.bool(True)
    )

    process.hltL2MuonSeedsFromL1TkMuon = cms.EDProducer("L2MuonSeedGeneratorFromL1TkMu",
        EtaMatchingBins = cms.vdouble(0.0, 2.5),
        InputObjects = cms.InputTag("L1TkMuons", "", processName),
        L1MaxEta = cms.double(2.5),
        L1MinPt = cms.double(0.0),
        MatchDR = cms.vdouble(0.3),
        OfflineSeedLabel = cms.untracked.InputTag("hltL2OfflineMuonSeeds"),
        Propagator = cms.string('SteppingHelixPropagatorAny'),
        ServiceParameters = cms.PSet(
            Propagators = cms.untracked.vstring('SteppingHelixPropagatorAny'),
            RPCLayers = cms.bool(True),
            UseMuonNavigation = cms.untracked.bool(True)
        ),
        SetMinPtBarrelTo = cms.double(3.5),
        SetMinPtEndcapTo = cms.double(1.0),
        MinPL1Tk = cms.double(3.5),
        MinPtL1TkBarrel = cms.double(3.5),
        UseOfflineSeed = cms.untracked.bool(True),
        UseUnassociatedL1 = cms.bool(False)
    )

    process.hltL2MuonsFromL1TkMuon = cms.EDProducer("L2MuonProducer",
        DoSeedRefit = cms.bool(False),
        InputObjects = cms.InputTag("hltL2MuonSeedsFromL1TkMuon"),
        L2TrajBuilderParameters = cms.PSet(
            BWFilterParameters = cms.PSet(
                BWSeedType = cms.string('fromGenerator'),
                CSCRecSegmentLabel = cms.InputTag("hltCscSegments"),
                DTRecSegmentLabel = cms.InputTag("hltDt4DSegments"),
                EnableCSCMeasurement = cms.bool(True),
                EnableDTMeasurement = cms.bool(True),
                EnableRPCMeasurement = cms.bool(True),
                FitDirection = cms.string('outsideIn'),
                MaxChi2 = cms.double(100.0),
                MuonTrajectoryUpdatorParameters = cms.PSet(
                    ExcludeRPCFromFit = cms.bool(False),
                    Granularity = cms.int32(0),
                    MaxChi2 = cms.double(25.0),
                    RescaleError = cms.bool(False),
                    RescaleErrorFactor = cms.double(100.0),
                    UseInvalidHits = cms.bool(True)
                ),
                NumberOfSigma = cms.double(3.0),
                Propagator = cms.string('hltESPFastSteppingHelixPropagatorAny'),
                RPCRecSegmentLabel = cms.InputTag("hltRpcRecHits")
            ),
            DoBackwardFilter = cms.bool(True),
            DoRefit = cms.bool(False),
            DoSeedRefit = cms.bool(False),
            FilterParameters = cms.PSet(
                CSCRecSegmentLabel = cms.InputTag("hltCscSegments"),
                DTRecSegmentLabel = cms.InputTag("hltDt4DSegments"),
                GEMRecSegmentLabel = cms.InputTag("hltGemRecHits"),
                ME0RecSegmentLabel = cms.InputTag("hltMe0Segments"),
                EnableCSCMeasurement = cms.bool(True),
                EnableDTMeasurement = cms.bool(True),
                EnableRPCMeasurement = cms.bool(True),
                EnableGEMMeasurement = cms.bool(True),
                EnableME0Measurement = cms.bool(True),
                FitDirection = cms.string('insideOut'),
                MaxChi2 = cms.double(1000.0),
                MuonTrajectoryUpdatorParameters = cms.PSet(
                    ExcludeRPCFromFit = cms.bool(False),
                    Granularity = cms.int32(0),
                    MaxChi2 = cms.double(25.0),
                    RescaleError = cms.bool(False),
                    RescaleErrorFactor = cms.double(100.0),
                    UseInvalidHits = cms.bool(True)
                ),
                NumberOfSigma = cms.double(3.0),
                Propagator = cms.string('hltESPFastSteppingHelixPropagatorAny'),
                RPCRecSegmentLabel = cms.InputTag("hltRpcRecHits")
            ),
            NavigationType = cms.string('Standard'),
            SeedPosition = cms.string('in'),
            SeedPropagator = cms.string('hltESPFastSteppingHelixPropagatorAny'),
            SeedTransformerParameters = cms.PSet(
                Fitter = cms.string('hltESPKFFittingSmootherForL2Muon'),
                MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
                NMinRecHits = cms.uint32(2),
                Propagator = cms.string('hltESPFastSteppingHelixPropagatorAny'),
                RescaleError = cms.double(100.0),
                UseSubRecHits = cms.bool(False)
            )
        ),
        MuonTrajectoryBuilder = cms.string('Exhaustive'),
        SeedTransformerParameters = cms.PSet(
            Fitter = cms.string('hltESPKFFittingSmootherForL2Muon'),
            MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
            NMinRecHits = cms.uint32(2),
            Propagator = cms.string('hltESPFastSteppingHelixPropagatorAny'),
            RescaleError = cms.double(100.0),
            UseSubRecHits = cms.bool(False)
        ),
        ServiceParameters = cms.PSet(
            Propagators = cms.untracked.vstring(
                'hltESPFastSteppingHelixPropagatorAny', 
                'hltESPFastSteppingHelixPropagatorOpposite'
            ),
            RPCLayers = cms.bool(True),
            UseMuonNavigation = cms.untracked.bool(True)
        ),
        TrackLoaderParameters = cms.PSet(
            DoSmoothing = cms.bool(False),
            MuonUpdatorAtVertexParameters = cms.PSet(
                BeamSpotPosition = cms.vdouble(0.0, 0.0, 0.0),
                BeamSpotPositionErrors = cms.vdouble(0.1, 0.1, 5.3),
                MaxChi2 = cms.double(1000000.0),
                Propagator = cms.string('hltESPFastSteppingHelixPropagatorOpposite')
            ),
            Smoother = cms.string('hltESPKFTrajectorySmootherForMuonTrackLoader'),
            TTRHBuilder = cms.string('WithTrackAngle'),
            VertexConstraint = cms.bool(True),
            beamSpot = cms.InputTag("offlineBeamSpot")
        )
    )

    process.HLTL2muonrecoFromL1TkMuonNocandSequence = cms.Sequence(
        process.HLTMuonLocalRecoSequence+
        process.hltL2OfflineMuonSeeds+
        process.hltL2MuonSeedsFromL1TkMuon+
        process.hltL2MuonsFromL1TkMuon
    )

    process.hltL2MuonFromL1TkMuonCandidates = cms.EDProducer("L2MuonCandidateProducer",
        InputObjects = cms.InputTag("hltL2MuonsFromL1TkMuon","UpdatedAtVtx")
    )

    process.HLTL2muonrecoSequence = cms.Sequence(
        process.HLTL2muonrecoFromL1TkMuonNocandSequence+
        process.hltL2MuonFromL1TkMuonCandidates
    )

    return process

def customizeOI(process, processName = "MYHLT"):
    process.hltPhase2L3OISeedsFromL1TkMuons = cms.EDProducer("TSGForOIFromL1TkMu",
        MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
        SF1 = cms.double(3.0),
        SF2 = cms.double(4.0),
        SF3 = cms.double(5.0),
        SF4 = cms.double(7.0),
        SF5 = cms.double(10.0),
        SF6 = cms.double(2.0),
        UseHitLessSeeds = cms.bool(True),
        adjustErrorsDynamicallyForHitless = cms.bool(True),
        adjustErrorsDynamicallyForHits = cms.bool(False),
        debug = cms.untracked.bool(True),
        estimator = cms.string('hltESPChi2MeasurementEstimator100'),
        eta1 = cms.double(0.2),
        eta2 = cms.double(0.3),
        eta3 = cms.double(1.0),
        eta4 = cms.double(1.2),
        eta5 = cms.double(1.6),
        eta6 = cms.double(1.4),
        eta7 = cms.double(2.1),
        fixedErrorRescaleFactorForHitless = cms.double(2.0),
        fixedErrorRescaleFactorForHits = cms.double(1.0),
        hitsToTry = cms.int32(1),
        layersToTry = cms.int32(2),
        maxEtaForTOB = cms.double(1.8),
        maxHitSeeds = cms.uint32(1),
        maxHitlessSeeds = cms.uint32(5),
        maxSeeds = cms.uint32(20),
        minEtaForTEC = cms.double(0.7),
        pT1 = cms.double(13.0),
        pT2 = cms.double(30.0),
        pT3 = cms.double(70.0),
        propagatorName = cms.string('PropagatorWithMaterialParabolicMf'),
        src = cms.InputTag("L1TkMuons", "", processName),
        tsosDiff1 = cms.double(0.2),
        tsosDiff2 = cms.double(0.02)
    )

    process.hltPhase2L3OITrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
        MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
        NavigationSchool = cms.string('SimpleNavigationSchool'),
        RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
        #RedundantSeedCleaner = cms.string('none'),
        SimpleMagneticField = cms.string(''),
        TrajectoryBuilder = cms.string('CkfTrajectoryBuilder'),
        TrajectoryBuilderPSet = cms.PSet(
           refToPSet_ = cms.string('HLTPSetMuonCkfTrajectoryBuilder')
        ),
        TrajectoryCleaner = cms.string('muonSeededTrajectoryCleanerBySharedHits'),
        #TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
        TransientInitialStateEstimatorParameters = cms.PSet(
            numberMeasurementsForFit = cms.int32(4),
            propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
            propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
        ),
        cleanTrajectoryAfterInOut = cms.bool(False),
        doSeedingRegionRebuilding = cms.bool(False),
        maxNSeeds = cms.uint32(500000),
        maxSeedsBeforeCleaning = cms.uint32(5000),
        src = cms.InputTag("hltPhase2L3OISeedsFromL1TkMuons"),
        useHitsSplitting = cms.bool(False)
    )

    process.hltPhase2L3OIMuCtfWithMaterialTracks = cms.EDProducer("TrackProducer",
        AlgorithmName = cms.string('iter10'),
        Fitter = cms.string('FlexibleKFFittingSmoother'),
        GeometricInnerState = cms.bool(True),
        MeasurementTracker = cms.string(''),
        MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
        NavigationSchool = cms.string(''),
        Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
        SimpleMagneticField = cms.string(''),
        TTRHBuilder = cms.string('WithTrackAngle'),
        TrajectoryInEvent = cms.bool(False),
        alias = cms.untracked.string('ctfWithMaterialTracks'),
        beamSpot = cms.InputTag("offlineBeamSpot"),
        clusterRemovalInfo = cms.InputTag(""),
        src = cms.InputTag("hltPhase2L3OITrackCandidates"),
        useHitsSplitting = cms.bool(False),
        useSimpleMF = cms.bool(False)
    )

    process.hltPhase2L3OIMuonTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
        beamspot = cms.InputTag("offlineBeamSpot"),
        ignoreVertices = cms.bool(True),
        mva = cms.PSet(
            dr_par = cms.PSet(
                d0err = cms.vdouble(0.003, 0.003, 3.40282346639e+38),
                d0err_par = cms.vdouble(0.001, 0.001, 3.40282346639e+38),
                dr_exp = cms.vint32(4, 4, 2147483647),
                dr_par1 = cms.vdouble(0.4, 0.4, 3.40282346639e+38),
                dr_par2 = cms.vdouble(0.3, 0.3, 3.40282346639e+38)
            ),
            dz_par = cms.PSet(
                dz_exp = cms.vint32(4, 4, 2147483647),
                dz_par1 = cms.vdouble(0.4, 0.4, 3.40282346639e+38),
                dz_par2 = cms.vdouble(0.35, 0.35, 3.40282346639e+38)
            ),
            maxChi2 = cms.vdouble(3.40282346639e+38, 3.40282346639e+38, 3.40282346639e+38),
            maxChi2n = cms.vdouble(10.0, 1.0, 0.4),
            maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
            maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
            maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 100.0),
            maxLostLayers = cms.vint32(4, 3, 2),
            min3DLayers = cms.vint32(1, 2, 1),
            minLayers = cms.vint32(3, 5, 5),
            minNVtxTrk = cms.int32(3),
            minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
            minPixelHits = cms.vint32(0, 0, 1)
        ),
        qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
        src = cms.InputTag("hltPhase2L3OIMuCtfWithMaterialTracks"),
        vertices = cms.InputTag("Notused")
    )

    process.hltPhase2L3OIMuonTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
        copyExtras = cms.untracked.bool(True),
        copyTrajectories = cms.untracked.bool(False),
        minQuality = cms.string('highPurity'),
        originalMVAVals = cms.InputTag("hltPhase2L3OIMuonTrackCutClassifier","MVAValues"),
        originalQualVals = cms.InputTag("hltPhase2L3OIMuonTrackCutClassifier","QualityMasks"),
        originalSource = cms.InputTag("hltPhase2L3OIMuCtfWithMaterialTracks")
    )

    process.hltL3MuonsPhase2L3OI = cms.EDProducer("L3MuonProducer",
        L3TrajBuilderParameters = cms.PSet(
            GlbRefitterParameters = cms.PSet(
                CSCRecSegmentLabel = cms.InputTag("hltCscSegments"),
                Chi2CutCSC = cms.double(150.0),
                Chi2CutDT = cms.double(10.0),
                Chi2CutRPC = cms.double(1.0),
                DTRecSegmentLabel = cms.InputTag("hltDt4DSegments"),
                DYTthrs = cms.vint32(30, 15),
                DoPredictionsOnly = cms.bool(False),
                Fitter = cms.string('hltESPL3MuKFTrajectoryFitter'),
                HitThreshold = cms.int32(1),
                MuonHitsOption = cms.int32(1),
                MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
                PropDirForCosmics = cms.bool(False),
                Propagator = cms.string('hltESPSmartPropagatorAny'),
                RefitDirection = cms.string('insideOut'),
                RefitFlag = cms.bool(True),
                RefitRPCHits = cms.bool(True),
                SkipStation = cms.int32(-1),
                TrackerRecHitBuilder = cms.string('WithTrackAngle'),
                TrackerSkipSection = cms.int32(-1),
                TrackerSkipSystem = cms.int32(-1)
            ),
            GlobalMuonTrackMatcher = cms.PSet(
                Chi2Cut_1 = cms.double(50.0),
                Chi2Cut_2 = cms.double(50.0),
                Chi2Cut_3 = cms.double(200.0),
                DeltaDCut_1 = cms.double(40.0),
                DeltaDCut_2 = cms.double(10.0),
                DeltaDCut_3 = cms.double(15.0),
                DeltaRCut_1 = cms.double(0.1),
                DeltaRCut_2 = cms.double(0.2),
                DeltaRCut_3 = cms.double(1.0),
                Eta_threshold = cms.double(1.2),
                LocChi2Cut = cms.double(0.001),
                MinP = cms.double(2.5),
                MinPt = cms.double(1.0),
                Propagator = cms.string('hltESPSmartPropagator'),
                Pt_threshold1 = cms.double(0.0),
                Pt_threshold2 = cms.double(999999999.0),
                Quality_1 = cms.double(20.0),
                Quality_2 = cms.double(15.0),
                Quality_3 = cms.double(7.0)
            ),
            MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
            MuonTrackingRegionBuilder = cms.PSet(
                DeltaEta = cms.double(0.2),
                DeltaPhi = cms.double(0.15),
                DeltaR = cms.double(0.025),
                DeltaZ = cms.double(24.2),
                EtaR_UpperLimit_Par1 = cms.double(0.25),
                EtaR_UpperLimit_Par2 = cms.double(0.15),
                Eta_fixed = cms.bool(True),
                Eta_min = cms.double(0.1),
                MeasurementTrackerName = cms.InputTag("hltESPMeasurementTracker"),
                OnDemand = cms.int32(-1),
                PhiR_UpperLimit_Par1 = cms.double(0.6),
                PhiR_UpperLimit_Par2 = cms.double(0.2),
                Phi_fixed = cms.bool(True),
                Phi_min = cms.double(0.1),
                Pt_fixed = cms.bool(False),
                Pt_min = cms.double(3.0),
                Rescale_Dz = cms.double(4.0),
                Rescale_eta = cms.double(3.0),
                Rescale_phi = cms.double(3.0),
                UseVertex = cms.bool(False),
                Z_fixed = cms.bool(False),
                beamSpot = cms.InputTag("offlineBeamSpot"),
                input = cms.InputTag("hltL2MuonsFromL1TkMuon","UpdatedAtVtx"),
                maxRegions = cms.int32(2),
                precise = cms.bool(True),
                vertexCollection = cms.InputTag("pixelVertices")
            ),
            PCut = cms.double(2.5),
            PtCut = cms.double(1.0),
            RefitRPCHits = cms.bool(True),
            ScaleTECxFactor = cms.double(-1.0),
            ScaleTECyFactor = cms.double(-1.0),
            TrackTransformer = cms.PSet(
                DoPredictionsOnly = cms.bool(False),
                Fitter = cms.string('hltESPL3MuKFTrajectoryFitter'),
                MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
                Propagator = cms.string('hltESPSmartPropagatorAny'),
                RefitDirection = cms.string('insideOut'),
                RefitRPCHits = cms.bool(True),
                Smoother = cms.string('hltESPKFTrajectorySmootherForMuonTrackLoader'),
                TrackerRecHitBuilder = cms.string('WithTrackAngle')
            ),
            TrackerPropagator = cms.string('SteppingHelixPropagatorAny'),
            TrackerRecHitBuilder = cms.string('WithTrackAngle'),
            tkTrajBeamSpot = cms.InputTag("offlineBeamSpot"),
            tkTrajLabel = cms.InputTag("hltPhase2L3OIMuonTrackSelectionHighPurity"),
            tkTrajMaxChi2 = cms.double(9999.0),
            tkTrajMaxDXYBeamSpot = cms.double(9999.0),
            tkTrajUseVertex = cms.bool(False),
            tkTrajVertex = cms.InputTag("Notused")
        ),
        MuonCollectionLabel = cms.InputTag("hltL2MuonsFromL1TkMuon","UpdatedAtVtx"),
        ServiceParameters = cms.PSet(
            Propagators = cms.untracked.vstring(
                'hltESPSmartPropagatorAny', 
                'SteppingHelixPropagatorAny', 
                'hltESPSmartPropagator', 
                'hltESPSteppingHelixPropagatorOpposite'
            ),
            RPCLayers = cms.bool(True),
            UseMuonNavigation = cms.untracked.bool(True)
        ),
        TrackLoaderParameters = cms.PSet(
            DoSmoothing = cms.bool(True),
            MuonSeededTracksInstance = cms.untracked.string('L2Seeded'),
            MuonUpdatorAtVertexParameters = cms.PSet(
                BeamSpotPositionErrors = cms.vdouble(0.1, 0.1, 5.3),
                MaxChi2 = cms.double(1000000.0),
                Propagator = cms.string('hltESPSteppingHelixPropagatorOpposite')
            ),
            PutTkTrackIntoEvent = cms.untracked.bool(False),
            SmoothTkTrack = cms.untracked.bool(False),
            Smoother = cms.string('hltESPKFTrajectorySmootherForMuonTrackLoader'),
            TTRHBuilder = cms.string('WithTrackAngle'),
            VertexConstraint = cms.bool(False),
            beamSpot = cms.InputTag("offlineBeamSpot")
        )
    )

    process.hltPhase2L3OIL3MuonsLinksCombination = cms.EDProducer("L3TrackLinksCombiner",
        labels = cms.VInputTag("hltL3MuonsPhase2L3OI")
    )

    process.hltPhase2L3OIL3Muons = cms.EDProducer("L3TrackCombiner",
        labels = cms.VInputTag("hltL3MuonsPhase2L3OI")
    )

    process.hltPhase2L3OIL3MuonCandidates = cms.EDProducer("L3MuonCandidateProducer",
        InputLinksObjects = cms.InputTag("hltPhase2L3OIL3MuonsLinksCombination"),
        InputObjects = cms.InputTag("hltPhase2L3OIL3Muons"),
        MuonPtOption = cms.string('Tracker')
    )

    process.HLTPhase2L3OImuonTkCandidateSequence = cms.Sequence(
        process.hltPhase2L3OISeedsFromL1TkMuons+
        process.hltPhase2L3OITrackCandidates+
        process.hltPhase2L3OIMuCtfWithMaterialTracks+
        process.hltPhase2L3OIMuonTrackCutClassifier+
        process.hltPhase2L3OIMuonTrackSelectionHighPurity+
        process.hltL3MuonsPhase2L3OI+
        process.hltPhase2L3OIL3MuonsLinksCombination+
        process.hltPhase2L3OIL3Muons+
        process.hltPhase2L3OIL3MuonCandidates
    )

    return process

def customizeIO(process, processName = "MYHLT"):

    # -- Iter0 quadruplet -- #
    process.hltPhase2L3MuonPixelTracksFilter = cms.EDProducer("PixelTrackFilterByKinematicsProducer",
        chi2 = cms.double(1000.0),
        nSigmaInvPtTolerance = cms.double(0.0),
        nSigmaTipMaxTolerance = cms.double(0.0),
        ptMin = cms.double(0.9),  ## before was 0.1
        tipMax = cms.double(1.0)
    )

    process.hltPhase2L3MuonPixelTracksFitter = cms.EDProducer("PixelFitterByHelixProjectionsProducer",
        scaleErrorsForBPix1 = cms.bool(False),
        scaleFactor = cms.double(0.65)
    )

    process.hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions = cms.EDProducer("CandidateSeededTrackingRegionsEDProducer",
        RegionPSet = cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            deltaEta = cms.double(0.035),        # modified at the end
            deltaPhi = cms.double(0.02),         # modified at the end
            input = cms.InputTag("L1TkMuons", "", processName),  # new L1 input
            maxNRegions = cms.int32(10000),      # was 2
            maxNVertices = cms.int32(1),
            measurementTrackerName = cms.InputTag(""),
            mode = cms.string('BeamSpotSigma'),
            nSigmaZBeamSpot = cms.double(4.0),
            nSigmaZVertex = cms.double(3.0),
            originRadius = cms.double(0.2),
            precise = cms.bool(True),
            ptMin = cms.double(2.0),            # modified at the end
            searchOpt = cms.bool(False),
            vertexCollection = cms.InputTag("notUsed"),
            whereToUseMeasurementTracker = cms.string('Never'),
            zErrorBeamSpot = cms.double(24.2),
            zErrorVetex = cms.double(0.2)
        )
    )

    process.hltPhase2L3FromL1TkMuonPixelLayerQuadruplets = cms.EDProducer("SeedingLayersEDProducer",
        BPix = cms.PSet(
            HitProducer = cms.string('siPixelRecHits'),
            TTRHBuilder = cms.string('WithTrackAngle'),
            #hitErrorRPhi = cms.double(0.0027),
            #hitErrorRZ = cms.double(0.006),
            #useErrorsFromParam = cms.bool(True)
        ),
        FPix = cms.PSet(
            HitProducer = cms.string('siPixelRecHits'),
            TTRHBuilder = cms.string('WithTrackAngle'),
            #hitErrorRPhi = cms.double(0.0051),
            #hitErrorRZ = cms.double(0.0036),
            #useErrorsFromParam = cms.bool(True)
        ),
        MTEC = cms.PSet(

        ),
        MTIB = cms.PSet(

        ),
        MTID = cms.PSet(

        ),
        MTOB = cms.PSet(

        ),
        TEC = cms.PSet(

        ),
        TIB = cms.PSet(

        ),
        TID = cms.PSet(

        ),
        TOB = cms.PSet(

        ),
        layerList = cms.vstring(
            'BPix1+BPix2+BPix3+BPix4', 
            'BPix1+BPix2+BPix3+FPix1_pos', 
            'BPix1+BPix2+BPix3+FPix1_neg', 
            'BPix1+BPix2+FPix1_pos+FPix2_pos', 
            'BPix1+BPix2+FPix1_neg+FPix2_neg', 
            'BPix1+FPix1_pos+FPix2_pos+FPix3_pos', 
            'BPix1+FPix1_neg+FPix2_neg+FPix3_neg'
        )
    )

    process.hltPhase2L3FromL1TkMuonPixelTracksHitDoublets = cms.EDProducer("HitPairEDProducer",
        clusterCheck = cms.InputTag(""),
        layerPairs = cms.vuint32(0, 1, 2),
        maxElement = cms.uint32(0),
        produceIntermediateHitDoublets = cms.bool(True),
        produceSeedingHitSets = cms.bool(False),
        seedingLayers = cms.InputTag("hltPhase2L3FromL1TkMuonPixelLayerQuadruplets"),
        trackingRegions = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions"),
        trackingRegionsSeedingLayers = cms.InputTag("")
    )

    process.hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets = cms.EDProducer("CAHitQuadrupletEDProducer",
        CAHardPtCut = cms.double(0.0),
        CAPhiCut = cms.double(0.2),
        CAThetaCut = cms.double(0.005),
        SeedComparitorPSet = cms.PSet(
            ComponentName = cms.string('LowPtClusterShapeSeedComparitor'),
            clusterShapeCacheSrc = cms.InputTag("siPixelClusterShapeCache"),
            clusterShapeHitFilter = cms.string('ClusterShapeHitFilter')
        ),
        doublets = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracksHitDoublets"),
        extraHitRPhitolerance = cms.double(0.032),
        fitFastCircle = cms.bool(True),
        fitFastCircleChi2Cut = cms.bool(True),
        maxChi2 = cms.PSet(
            enabled = cms.bool(True),
            pt1 = cms.double(0.7),
            pt2 = cms.double(2.0),
            value1 = cms.double(200.0),
            value2 = cms.double(50.0)
        ),
        useBendingCorrection = cms.bool(True)
    )

    process.hltPhase2L3FromL1TkMuonPixelTracks = cms.EDProducer("PixelTrackProducer",
        Cleaner = cms.string('hltPixelTracksCleanerBySharedHits'),
        Filter = cms.InputTag("hltPhase2L3MuonPixelTracksFilter"),
        Fitter = cms.InputTag("hltPhase2L3MuonPixelTracksFitter"),
        SeedingHitSets = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets"),
        passLabel = cms.string('')
    )

    process.hltPhase2L3FromL1TkMuonPixelVertices = cms.EDProducer("PixelVertexProducer",
        Finder = cms.string('DivisiveVertexFinder'),
        Method2 = cms.bool(True),
        NTrkMin = cms.int32(2),
        PVcomparer = cms.PSet(
            refToPSet_ = cms.string('hltPhase2PSetPvClusterComparerForIT')
        ),
        PtMin = cms.double(1.0),
        TrackCollection = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracks"),
        UseError = cms.bool(True),
        Verbosity = cms.int32(0),
        WtAverage = cms.bool(True),
        ZOffset = cms.double(5.0),
        ZSeparation = cms.double(0.05),
        beamSpot = cms.InputTag("offlineBeamSpot")
    )

    process.hltPhase2L3FromL1TkMuonTrimmedPixelVertices = cms.EDProducer("PixelVertexCollectionTrimmer",
        PVcomparer = cms.PSet(
            refToPSet_ = cms.string('hltPhase2PSetPvClusterComparerForIT')
        ),
        fractionSumPt2 = cms.double(0.3),
        maxVtx = cms.uint32(100),
        minSumPt2 = cms.double(0.0),
        src = cms.InputTag("hltPhase2L3FromL1TkMuonPixelVertices")
    )

    process.hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracks = cms.EDProducer("SeedGeneratorFromProtoTracksEDProducer",
        InputCollection = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracks"),
        InputVertexCollection = cms.InputTag("hltPhase2L3FromL1TkMuonTrimmedPixelVertices"),
        SeedCreatorPSet = cms.PSet(
            refToPSet_ = cms.string('hltPhase2SeedFromProtoTracks')
        ),
        TTRHBuilder = cms.string('WithTrackAngle'),
        originHalfLength = cms.double(0.3),
        originRadius = cms.double(0.1),
        useEventsWithNoVertex = cms.bool(True),
        usePV = cms.bool(False),
        useProtoTrackKinematics = cms.bool(False)
    )

    process.hltIter0Phase2L3FromL1TkMuonCkfTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
        MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
        NavigationSchool = cms.string('SimpleNavigationSchool'),
        RedundantSeedCleaner = cms.string('none'),
        SimpleMagneticField = cms.string('ParabolicMf'),
        TrajectoryBuilder = cms.string('GroupedCkfTrajectoryBuilder'),
        TrajectoryBuilderPSet = cms.PSet(
            refToPSet_ = cms.string('HLTIter0Phase2L3FromL1TkMuonPSetGroupedCkfTrajectoryBuilderIT')
        ),
        TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
        TransientInitialStateEstimatorParameters = cms.PSet(
            numberMeasurementsForFit = cms.int32(4),
            propagatorAlongTISE = cms.string('PropagatorWithMaterialParabolicMf'),
            propagatorOppositeTISE = cms.string('PropagatorWithMaterialParabolicMfOpposite')
        ),
        cleanTrajectoryAfterInOut = cms.bool(False),
        doSeedingRegionRebuilding = cms.bool(True),
        maxNSeeds = cms.uint32(100000),
        maxSeedsBeforeCleaning = cms.uint32(1000),
        src = cms.InputTag("hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracks"),
        useHitsSplitting = cms.bool(True)
    )

    process.hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks = cms.EDProducer("TrackProducer",
        AlgorithmName = cms.string('hltIter0'),
        Fitter = cms.string('FlexibleKFFittingSmoother'),
        GeometricInnerState = cms.bool(True),
        MeasurementTracker = cms.string(''),
        MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
        NavigationSchool = cms.string(''),
        Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
        SimpleMagneticField = cms.string(''),
        TTRHBuilder = cms.string('WithTrackAngle'),
        TrajectoryInEvent = cms.bool(False),
        alias = cms.untracked.string('ctfWithMaterialTracks'),
        beamSpot = cms.InputTag("offlineBeamSpot"),
        clusterRemovalInfo = cms.InputTag(""),
        src = cms.InputTag("hltIter0Phase2L3FromL1TkMuonCkfTrackCandidates"),
        useHitsSplitting = cms.bool(False),
        useSimpleMF = cms.bool(False)
    )

    process.hltIter0Phase2L3FromL1TkMuonTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
        beamspot = cms.InputTag("offlineBeamSpot"),
        ignoreVertices = cms.bool(False),
        mva = cms.PSet(
            dr_par = cms.PSet(
                d0err = cms.vdouble(0.003, 0.003, 3.40282346639e+38),
                d0err_par = cms.vdouble(0.001, 0.001, 3.40282346639e+38),
                dr_exp = cms.vint32(4, 4, 2147483647),
                dr_par1 = cms.vdouble(0.4, 0.4, 3.40282346639e+38),
                dr_par2 = cms.vdouble(0.3, 0.3, 3.40282346639e+38)
            ),
            dz_par = cms.PSet(
                dz_exp = cms.vint32(4, 4, 2147483647),
                dz_par1 = cms.vdouble(0.4, 0.4, 3.40282346639e+38),
                dz_par2 = cms.vdouble(0.35, 0.35, 3.40282346639e+38)
            ),
            maxChi2 = cms.vdouble(3.40282346639e+38, 3.40282346639e+38, 3.40282346639e+38),
            maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
            maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
            maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
            maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 100.0),
            maxLostLayers = cms.vint32(1, 1, 1),
            min3DLayers = cms.vint32(0, 3, 4),
            minLayers = cms.vint32(3, 3, 4),
            minNVtxTrk = cms.int32(3),
            minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
            minPixelHits = cms.vint32(0, 3, 4)
        ),
        qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
        src = cms.InputTag("hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks"),
        vertices = cms.InputTag("hltPhase2L3FromL1TkMuonTrimmedPixelVertices")
    )

    process.hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
        copyExtras = cms.untracked.bool(True),
        copyTrajectories = cms.untracked.bool(False),
        minQuality = cms.string('highPurity'),
        originalMVAVals = cms.InputTag("hltIter0Phase2L3FromL1TkMuonTrackCutClassifier","MVAValues"),
        originalQualVals = cms.InputTag("hltIter0Phase2L3FromL1TkMuonTrackCutClassifier","QualityMasks"),
        originalSource = cms.InputTag("hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks")
    )

    # -- Iter2 triplet -- #
    process.hltIter2Phase2L3FromL1TkMuonClustersRefRemoval = cms.EDProducer("TrackClusterRemoverPhase2",
        TrackQuality = cms.string('highPurity'),
        maxChi2 = cms.double(16.0),
        minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
        oldClusterRemovalInfo = cms.InputTag(""),
        overrideTrkQuals = cms.InputTag(""),
        phase2pixelClusters = cms.InputTag("siPixelClusters"),
        phase2OTClusters = cms.InputTag("siPhase2Clusters"),
        trackClassifier = cms.InputTag("","QualityMasks"),
        trajectories = cms.InputTag("hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity")
    )

    process.hltIter2Phase2L3FromL1TkMuonMaskedMeasurementTrackerEvent = cms.EDProducer("MaskedMeasurementTrackerEventProducer",
        OnDemand = cms.bool(False),
        phase2clustersToSkip = cms.InputTag("hltIter2Phase2L3FromL1TkMuonClustersRefRemoval"),
        src = cms.InputTag("MeasurementTrackerEvent")
    )

    process.hltIter2Phase2L3FromL1TkMuonPixelLayerTriplets = cms.EDProducer("SeedingLayersEDProducer",
        BPix = cms.PSet(
            HitProducer = cms.string('siPixelRecHits'),
            TTRHBuilder = cms.string('WithTrackAngle'),
            #hitErrorRPhi = cms.double(0.0027),
            #hitErrorRZ = cms.double(0.006),
            skipClusters = cms.InputTag("hltIter2Phase2L3FromL1TkMuonClustersRefRemoval"),
            #useErrorsFromParam = cms.bool(True)
        ),
        FPix = cms.PSet(
            HitProducer = cms.string('siPixelRecHits'),
            TTRHBuilder = cms.string('WithTrackAngle'),
            #hitErrorRPhi = cms.double(0.0051),
            #hitErrorRZ = cms.double(0.0036),
            skipClusters = cms.InputTag("hltIter2Phase2L3FromL1TkMuonClustersRefRemoval"),
            #useErrorsFromParam = cms.bool(True)
        ),
        MTEC = cms.PSet(

        ),
        MTIB = cms.PSet(

        ),
        MTID = cms.PSet(

        ),
        MTOB = cms.PSet(

        ),
        TEC = cms.PSet(

        ),
        TIB = cms.PSet(

        ),
        TID = cms.PSet(

        ),
        TOB = cms.PSet(

        ),
        layerList = cms.vstring(
            'BPix1+BPix2+BPix3', 
            'BPix2+BPix3+BPix4', 
            'BPix1+BPix3+BPix4', 
            'BPix1+BPix2+BPix4', 
            'BPix2+BPix3+FPix1_pos', 
            'BPix2+BPix3+FPix1_neg', 
            'BPix1+BPix2+FPix1_pos', 
            'BPix1+BPix2+FPix1_neg', 
            'BPix2+FPix1_pos+FPix2_pos', 
            'BPix2+FPix1_neg+FPix2_neg', 
            'BPix1+FPix1_pos+FPix2_pos', 
            'BPix1+FPix1_neg+FPix2_neg', 
            'FPix1_pos+FPix2_pos+FPix3_pos', 
            'FPix1_neg+FPix2_neg+FPix3_neg', 
            'BPix1+BPix3+FPix1_pos', 
            'BPix1+BPix2+FPix2_pos', 
            'BPix1+BPix3+FPix1_neg', 
            'BPix1+BPix2+FPix2_neg', 
            'BPix1+FPix2_neg+FPix3_neg', 
            'BPix1+FPix1_neg+FPix3_neg', 
            'BPix1+FPix2_pos+FPix3_pos', 
            'BPix1+FPix1_pos+FPix3_pos'
        )
    )

    process.hltIter2Phase2L3FromL1TkMuonPixelClusterCheck = cms.EDProducer("ClusterCheckerEDProducer",
        ClusterCollectionLabel = cms.InputTag("MeasurementTrackerEvent"),
        MaxNumberOfCosmicClusters = cms.uint32(50000),
        MaxNumberOfPixelClusters = cms.uint32(10000),
        PixelClusterCollectionLabel = cms.InputTag("siPixelClusters"),
        cut = cms.string(''),
        doClusterCheck = cms.bool(False),
        silentClusterCheck = cms.untracked.bool(False)
    )

    process.hltIter2Phase2L3FromL1TkMuonPixelHitDoublets = cms.EDProducer("HitPairEDProducer",
        clusterCheck = cms.InputTag("hltIter2Phase2L3FromL1TkMuonPixelClusterCheck"),
        layerPairs = cms.vuint32(0, 1),
        maxElement = cms.uint32(0),
        produceIntermediateHitDoublets = cms.bool(True),
        produceSeedingHitSets = cms.bool(False),
        seedingLayers = cms.InputTag("hltIter2Phase2L3FromL1TkMuonPixelLayerTriplets"),
        trackingRegions = cms.InputTag("hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions"),
        trackingRegionsSeedingLayers = cms.InputTag("")
    )

    process.hltIter2Phase2L3FromL1TkMuonPixelHitTriplets = cms.EDProducer("CAHitTripletEDProducer",
        CAHardPtCut = cms.double(0.3),
        CAPhiCut = cms.double(0.1),
        CAThetaCut = cms.double(0.015),
        SeedComparitorPSet = cms.PSet(
            ComponentName = cms.string('none')
        ),
        doublets = cms.InputTag("hltIter2Phase2L3FromL1TkMuonPixelHitDoublets"),
        extraHitRPhitolerance = cms.double(0.032),
        maxChi2 = cms.PSet(
            enabled = cms.bool(True),
            pt1 = cms.double(0.8),
            pt2 = cms.double(8.0),
            value1 = cms.double(100.0),
            value2 = cms.double(6.0)
        ),
        useBendingCorrection = cms.bool(True)
    )

    process.hltIter2Phase2L3FromL1TkMuonPixelSeeds = cms.EDProducer("SeedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer",
        MinOneOverPtError = cms.double(1.0),
        OriginTransverseErrorMultiplier = cms.double(1.0),
        SeedComparitorPSet = cms.PSet(
            ComponentName = cms.string('none')
        ),
        SeedMomentumForBOFF = cms.double(5.0),
        TTRHBuilder = cms.string('WithTrackAngle'),
        forceKinematicWithRegionDirection = cms.bool(False),
        magneticField = cms.string('ParabolicMf'),
        propagator = cms.string('PropagatorWithMaterialParabolicMf'),
        seedingHitSets = cms.InputTag("hltIter2Phase2L3FromL1TkMuonPixelHitTriplets")
    )

    process.hltIter2Phase2L3FromL1TkMuonCkfTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
        MeasurementTrackerEvent = cms.InputTag("hltIter2Phase2L3FromL1TkMuonMaskedMeasurementTrackerEvent"),
        NavigationSchool = cms.string('SimpleNavigationSchool'),
        RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
        SimpleMagneticField = cms.string('ParabolicMf'),
        TrajectoryBuilder = cms.string(''),
        TrajectoryBuilderPSet = cms.PSet(
            refToPSet_ = cms.string('HLTIter2Phase2L3FromL1TkMuonPSetGroupedCkfTrajectoryBuilderIT')
        ),
        TrajectoryCleaner = cms.string('hltESPTrajectoryCleanerBySharedHits'),
        TransientInitialStateEstimatorParameters = cms.PSet(
            numberMeasurementsForFit = cms.int32(4),
            propagatorAlongTISE = cms.string('PropagatorWithMaterialParabolicMf'),
            propagatorOppositeTISE = cms.string('PropagatorWithMaterialParabolicMfOpposite')
        ),
        cleanTrajectoryAfterInOut = cms.bool(False),
        doSeedingRegionRebuilding = cms.bool(False),
        maxNSeeds = cms.uint32(100000),
        maxSeedsBeforeCleaning = cms.uint32(1000),
        src = cms.InputTag("hltIter2Phase2L3FromL1TkMuonPixelSeeds"),
        useHitsSplitting = cms.bool(False)
    )

    process.hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracks = cms.EDProducer("TrackProducer",
        AlgorithmName = cms.string('hltIter2'),
        Fitter = cms.string('FlexibleKFFittingSmoother'),
        GeometricInnerState = cms.bool(True),
        MeasurementTracker = cms.string(''),
        MeasurementTrackerEvent = cms.InputTag("hltIter2Phase2L3FromL1TkMuonMaskedMeasurementTrackerEvent"),
        NavigationSchool = cms.string(''),
        Propagator = cms.string('hltESPRungeKuttaTrackerPropagator'),
        SimpleMagneticField = cms.string(''),
        TTRHBuilder = cms.string('WithTrackAngle'),
        TrajectoryInEvent = cms.bool(False),
        alias = cms.untracked.string('ctfWithMaterialTracks'),
        beamSpot = cms.InputTag("offlineBeamSpot"),
        clusterRemovalInfo = cms.InputTag(""),
        src = cms.InputTag("hltIter2Phase2L3FromL1TkMuonCkfTrackCandidates"),
        useHitsSplitting = cms.bool(False),
        useSimpleMF = cms.bool(False)
    )

    process.hltIter2Phase2L3FromL1TkMuonTrackCutClassifier = cms.EDProducer("TrackCutClassifier",
        beamspot = cms.InputTag("offlineBeamSpot"),
        ignoreVertices = cms.bool(False),
        mva = cms.PSet(
            dr_par = cms.PSet(
                d0err = cms.vdouble(0.003, 0.003, 3.40282346639e+38),
                d0err_par = cms.vdouble(0.001, 0.001, 3.40282346639e+38),
                dr_exp = cms.vint32(4, 4, 2147483647),
                dr_par1 = cms.vdouble(3.40282346639e+38, 0.4, 3.40282346639e+38),
                dr_par2 = cms.vdouble(3.40282346639e+38, 0.3, 3.40282346639e+38)
            ),
            dz_par = cms.PSet(
                dz_exp = cms.vint32(4, 4, 2147483647),
                dz_par1 = cms.vdouble(3.40282346639e+38, 0.4, 3.40282346639e+38),
                dz_par2 = cms.vdouble(3.40282346639e+38, 0.35, 3.40282346639e+38)
            ),
            maxChi2 = cms.vdouble(9999.0, 25.0, 3.40282346639e+38),
            maxChi2n = cms.vdouble(1.2, 1.0, 0.7),
            maxDr = cms.vdouble(0.5, 0.03, 3.40282346639e+38),
            maxDz = cms.vdouble(0.5, 0.2, 3.40282346639e+38),
            maxDzWrtBS = cms.vdouble(3.40282346639e+38, 24.0, 100.0),
            maxLostLayers = cms.vint32(1, 1, 1),
            min3DLayers = cms.vint32(0, 0, 0),
            minLayers = cms.vint32(3, 3, 3),
            minNVtxTrk = cms.int32(3),
            minNdof = cms.vdouble(1e-05, 1e-05, 1e-05),
            minPixelHits = cms.vint32(0, 0, 0)
        ),
        qualityCuts = cms.vdouble(-0.7, 0.1, 0.7),
        src = cms.InputTag("hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracks"),
        vertices = cms.InputTag("hltPhase2L3FromL1TkMuonTrimmedPixelVertices")
    )

    process.hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity = cms.EDProducer("TrackCollectionFilterCloner",
        copyExtras = cms.untracked.bool(True),
        copyTrajectories = cms.untracked.bool(False),
        minQuality = cms.string('highPurity'),
        originalMVAVals = cms.InputTag("hltIter2Phase2L3FromL1TkMuonTrackCutClassifier","MVAValues"),
        originalQualVals = cms.InputTag("hltIter2Phase2L3FromL1TkMuonTrackCutClassifier","QualityMasks"),
        originalSource = cms.InputTag("hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracks")
    )

    return process

def customizeL3Muon(process, processName = "MYHLT"):

    process.hltIter2Phase2L3FromL1TkMuonMerged = cms.EDProducer("TrackListMerger",
        Epsilon = cms.double(-0.001),
        FoundHitBonus = cms.double(5.0),
        LostHitPenalty = cms.double(20.0),
        MaxNormalizedChisq = cms.double(1000.0),
        MinFound = cms.int32(3),
        MinPT = cms.double(0.05),
        ShareFrac = cms.double(0.19),
        TrackProducers = cms.VInputTag("hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity", "hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity"),
        allowFirstHitShare = cms.bool(True),
        copyExtras = cms.untracked.bool(True),
        copyMVA = cms.bool(False),
        hasSelector = cms.vint32(0, 0),
        indivShareFrac = cms.vdouble(1.0, 1.0),
        newQuality = cms.string('confirmed'),
        selectedTrackQuals = cms.VInputTag("hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity", "hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity"),
        setsToMerge = cms.VPSet(cms.PSet(
            pQual = cms.bool(False),
            tLists = cms.vint32(0, 1)
        )),
        trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
        writeOnlyTrkQuals = cms.bool(False)
    )

    process.hltPhase2L3MuonMerged = cms.EDProducer("TrackListMerger",
        Epsilon = cms.double(-0.001),
        FoundHitBonus = cms.double(5.0),
        LostHitPenalty = cms.double(20.0),
        MaxNormalizedChisq = cms.double(1000.0),
        MinFound = cms.int32(3),
        MinPT = cms.double(0.05),
        ShareFrac = cms.double(0.19),
        TrackProducers = cms.VInputTag("hltPhase2L3OIMuonTrackSelectionHighPurity", "hltIter2Phase2L3FromL1TkMuonMerged"),
        allowFirstHitShare = cms.bool(True),
        copyExtras = cms.untracked.bool(True),
        copyMVA = cms.bool(False),
        hasSelector = cms.vint32(0, 0),
        indivShareFrac = cms.vdouble(1.0, 1.0),
        newQuality = cms.string('confirmed'),
        selectedTrackQuals = cms.VInputTag("hltPhase2L3OIMuonTrackSelectionHighPurity", "hltIter2Phase2L3FromL1TkMuonMerged"),
        setsToMerge = cms.VPSet(cms.PSet(
            pQual = cms.bool(False),
            tLists = cms.vint32(0, 1)
        )),
        trackAlgoPriorityOrder = cms.string('hltESPTrackAlgoPriorityOrder'),
        writeOnlyTrkQuals = cms.bool(False)
    )

    process.hltPhase2L3GlbMuon = cms.EDProducer("L3MuonProducer",
        L3TrajBuilderParameters = cms.PSet(
            GlbRefitterParameters = cms.PSet(
                CSCRecSegmentLabel = cms.InputTag("hltCscSegments"),
                Chi2CutCSC = cms.double(150.0),
                Chi2CutDT = cms.double(10.0),
                Chi2CutRPC = cms.double(1.0),
                DTRecSegmentLabel = cms.InputTag("hltDt4DSegments"),
                DYTthrs = cms.vint32(30, 15),
                DoPredictionsOnly = cms.bool(False),
                Fitter = cms.string('hltESPL3MuKFTrajectoryFitter'),
                HitThreshold = cms.int32(1),
                MuonHitsOption = cms.int32(1),
                MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
                PropDirForCosmics = cms.bool(False),
                Propagator = cms.string('hltESPSmartPropagatorAny'),
                RefitDirection = cms.string('insideOut'),
                RefitFlag = cms.bool(True),
                RefitRPCHits = cms.bool(True),
                SkipStation = cms.int32(-1),
                TrackerRecHitBuilder = cms.string('WithTrackAngle'),
                TrackerSkipSection = cms.int32(-1),
                TrackerSkipSystem = cms.int32(-1)
            ),
            GlobalMuonTrackMatcher = cms.PSet(
                Chi2Cut_1 = cms.double(50.0),
                Chi2Cut_2 = cms.double(50.0),
                Chi2Cut_3 = cms.double(200.0),
                DeltaDCut_1 = cms.double(40.0),
                DeltaDCut_2 = cms.double(10.0),
                DeltaDCut_3 = cms.double(15.0),
                DeltaRCut_1 = cms.double(0.1),
                DeltaRCut_2 = cms.double(0.2),
                DeltaRCut_3 = cms.double(1.0),
                Eta_threshold = cms.double(1.2),
                LocChi2Cut = cms.double(0.001),
                MinP = cms.double(2.5),
                MinPt = cms.double(1.0),
                Propagator = cms.string('hltESPSmartPropagator'),
                Pt_threshold1 = cms.double(0.0),
                Pt_threshold2 = cms.double(999999999.0),
                Quality_1 = cms.double(20.0),
                Quality_2 = cms.double(15.0),
                Quality_3 = cms.double(7.0)
            ),
            MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
            MuonTrackingRegionBuilder = cms.PSet(
                DeltaEta = cms.double(0.2),
                DeltaPhi = cms.double(0.15),
                DeltaR = cms.double(0.025),
                DeltaZ = cms.double(24.2),
                EtaR_UpperLimit_Par1 = cms.double(0.25),
                EtaR_UpperLimit_Par2 = cms.double(0.15),
                Eta_fixed = cms.bool(True),
                Eta_min = cms.double(0.1),
                MeasurementTrackerName = cms.InputTag("hltESPMeasurementTracker"),
                OnDemand = cms.int32(-1),
                PhiR_UpperLimit_Par1 = cms.double(0.6),
                PhiR_UpperLimit_Par2 = cms.double(0.2),
                Phi_fixed = cms.bool(True),
                Phi_min = cms.double(0.1),
                Pt_fixed = cms.bool(False),
                Pt_min = cms.double(3.0),
                Rescale_Dz = cms.double(4.0),
                Rescale_eta = cms.double(3.0),
                Rescale_phi = cms.double(3.0),
                UseVertex = cms.bool(False),
                Z_fixed = cms.bool(False),
                beamSpot = cms.InputTag("offlineBeamSpot"),
                input = cms.InputTag("hltL2MuonsFromL1TkMuon","UpdatedAtVtx"),
                maxRegions = cms.int32(2),
                precise = cms.bool(True),
                vertexCollection = cms.InputTag("pixelVertices")
            ),
            PCut = cms.double(2.5),
            PtCut = cms.double(1.0),
            RefitRPCHits = cms.bool(True),
            ScaleTECxFactor = cms.double(-1.0),
            ScaleTECyFactor = cms.double(-1.0),
            TrackTransformer = cms.PSet(
                DoPredictionsOnly = cms.bool(False),
                Fitter = cms.string('hltESPL3MuKFTrajectoryFitter'),
                MuonRecHitBuilder = cms.string('hltESPMuonTransientTrackingRecHitBuilder'),
                Propagator = cms.string('hltESPSmartPropagatorAny'),
                RefitDirection = cms.string('insideOut'),
                RefitRPCHits = cms.bool(True),
                Smoother = cms.string('hltESPKFTrajectorySmootherForMuonTrackLoader'),
                TrackerRecHitBuilder = cms.string('WithTrackAngle')
            ),
            TrackerPropagator = cms.string('SteppingHelixPropagatorAny'),
            TrackerRecHitBuilder = cms.string('WithTrackAngle'),
            tkTrajBeamSpot = cms.InputTag("offlineBeamSpot"),
            tkTrajLabel = cms.InputTag("hltPhase2L3MuonMerged"),
            tkTrajMaxChi2 = cms.double(9999.0),
            tkTrajMaxDXYBeamSpot = cms.double(9999.0),
            tkTrajUseVertex = cms.bool(False),
            tkTrajVertex = cms.InputTag("Notused")
        ),
        MuonCollectionLabel = cms.InputTag("hltL2MuonsFromL1TkMuon","UpdatedAtVtx"),
        ServiceParameters = cms.PSet(
            Propagators = cms.untracked.vstring(
                'hltESPSmartPropagatorAny', 
                'SteppingHelixPropagatorAny', 
                'hltESPSmartPropagator', 
                'hltESPSteppingHelixPropagatorOpposite'
            ),
            RPCLayers = cms.bool(True),
            UseMuonNavigation = cms.untracked.bool(True)
        ),
        TrackLoaderParameters = cms.PSet(
            DoSmoothing = cms.bool(True),
            MuonSeededTracksInstance = cms.untracked.string('L2Seeded'),
            MuonUpdatorAtVertexParameters = cms.PSet(
                BeamSpotPositionErrors = cms.vdouble(0.1, 0.1, 5.3),
                MaxChi2 = cms.double(1000000.0),
                Propagator = cms.string('hltESPSteppingHelixPropagatorOpposite')
            ),
            PutTkTrackIntoEvent = cms.untracked.bool(False),
            SmoothTkTrack = cms.untracked.bool(False),
            Smoother = cms.string('hltESPKFTrajectorySmootherForMuonTrackLoader'),
            TTRHBuilder = cms.string('WithTrackAngle'),
            VertexConstraint = cms.bool(False),
            beamSpot = cms.InputTag("offlineBeamSpot")
        )
    )

    process.hltPhase2L3MuonsNoID = cms.EDProducer("MuonIdProducer",
        CaloExtractorPSet = cms.PSet(
            CenterConeOnCalIntersection = cms.bool(False),
            ComponentName = cms.string('CaloExtractorByAssociator'),
            DR_Max = cms.double(1.0),
            DR_Veto_E = cms.double(0.07),
            DR_Veto_H = cms.double(0.1),
            DR_Veto_HO = cms.double(0.1),
            DepositInstanceLabels = cms.vstring(
                'ecal', 
                'hcal', 
                'ho'
            ),
            DepositLabel = cms.untracked.string('Cal'),
            NoiseTow_EB = cms.double(0.04),
            NoiseTow_EE = cms.double(0.15),
            Noise_EB = cms.double(0.025),
            Noise_EE = cms.double(0.1),
            Noise_HB = cms.double(0.2),
            Noise_HE = cms.double(0.2),
            Noise_HO = cms.double(0.2),
            PrintTimeReport = cms.untracked.bool(False),
            PropagatorName = cms.string('hltESPFastSteppingHelixPropagatorAny'),
            ServiceParameters = cms.PSet(
                Propagators = cms.untracked.vstring('hltESPFastSteppingHelixPropagatorAny'),
                RPCLayers = cms.bool(False),
                UseMuonNavigation = cms.untracked.bool(False)
            ),
            Threshold_E = cms.double(0.2),
            Threshold_H = cms.double(0.5),
            Threshold_HO = cms.double(0.5),
            TrackAssociatorParameters = cms.PSet(
                CSCSegmentCollectionLabel = cms.InputTag("hltCscSegments"),
                CaloTowerCollectionLabel = cms.InputTag("Notused"),
                DTRecSegment4DCollectionLabel = cms.InputTag("hltDt4DSegments"),
                EBRecHitCollectionLabel = cms.InputTag("Notused"),
                EERecHitCollectionLabel = cms.InputTag("Notused"),
                HBHERecHitCollectionLabel = cms.InputTag("Notused"),
                HORecHitCollectionLabel = cms.InputTag("Notused"),
                accountForTrajectoryChangeCalo = cms.bool(False),
                dREcal = cms.double(1.0),
                dREcalPreselection = cms.double(1.0),
                dRHcal = cms.double(1.0),
                dRHcalPreselection = cms.double(1.0),
                dRMuon = cms.double(9999.0),
                dRMuonPreselection = cms.double(0.2),
                dRPreshowerPreselection = cms.double(0.2),
                muonMaxDistanceSigmaX = cms.double(0.0),
                muonMaxDistanceSigmaY = cms.double(0.0),
                muonMaxDistanceX = cms.double(5.0),
                muonMaxDistanceY = cms.double(5.0),
                propagateAllDirections = cms.bool(True),
                trajectoryUncertaintyTolerance = cms.double(-1.0),
                truthMatch = cms.bool(False),
                useCalo = cms.bool(True),
                useEcal = cms.bool(False),
                useHO = cms.bool(False),
                useHcal = cms.bool(False),
                useMuon = cms.bool(False),
                usePreshower = cms.bool(False)
            ),
            UseRecHitsFlag = cms.bool(False)
        ),
        JetExtractorPSet = cms.PSet(
            ComponentName = cms.string('JetExtractor'),
            DR_Max = cms.double(1.0),
            DR_Veto = cms.double(0.1),
            ExcludeMuonVeto = cms.bool(True),
            JetCollectionLabel = cms.InputTag("Notused"),
            PrintTimeReport = cms.untracked.bool(False),
            PropagatorName = cms.string('hltESPFastSteppingHelixPropagatorAny'),
            ServiceParameters = cms.PSet(
                Propagators = cms.untracked.vstring('hltESPFastSteppingHelixPropagatorAny'),
                RPCLayers = cms.bool(False),
                UseMuonNavigation = cms.untracked.bool(False)
            ),
            Threshold = cms.double(5.0),
            TrackAssociatorParameters = cms.PSet(
                CSCSegmentCollectionLabel = cms.InputTag("hltCscSegments"),
                CaloTowerCollectionLabel = cms.InputTag("Notused"),
                DTRecSegment4DCollectionLabel = cms.InputTag("hltDt4DSegments"),
                EBRecHitCollectionLabel = cms.InputTag("Notused"),
                EERecHitCollectionLabel = cms.InputTag("Notused"),
                HBHERecHitCollectionLabel = cms.InputTag("Notused"),
                HORecHitCollectionLabel = cms.InputTag("Notused"),
                accountForTrajectoryChangeCalo = cms.bool(False),
                dREcal = cms.double(0.5),
                dREcalPreselection = cms.double(0.5),
                dRHcal = cms.double(0.5),
                dRHcalPreselection = cms.double(0.5),
                dRMuon = cms.double(9999.0),
                dRMuonPreselection = cms.double(0.2),
                dRPreshowerPreselection = cms.double(0.2),
                muonMaxDistanceSigmaX = cms.double(0.0),
                muonMaxDistanceSigmaY = cms.double(0.0),
                muonMaxDistanceX = cms.double(5.0),
                muonMaxDistanceY = cms.double(5.0),
                propagateAllDirections = cms.bool(True),
                trajectoryUncertaintyTolerance = cms.double(-1.0),
                truthMatch = cms.bool(False),
                useCalo = cms.bool(True),
                useEcal = cms.bool(False),
                useHO = cms.bool(False),
                useHcal = cms.bool(False),
                useMuon = cms.bool(False),
                usePreshower = cms.bool(False)
            )
        ),
        MuonCaloCompatibility = cms.PSet(
            MuonTemplateFileName = cms.FileInPath('RecoMuon/MuonIdentification/data/MuID_templates_muons_lowPt_3_1_norm.root'),
            PionTemplateFileName = cms.FileInPath('RecoMuon/MuonIdentification/data/MuID_templates_pions_lowPt_3_1_norm.root'),
            allSiPMHO = cms.bool(False),
            delta_eta = cms.double(0.02),
            delta_phi = cms.double(0.02)
        ),
        TimingFillerParameters = cms.PSet(
            CSCTimingParameters = cms.PSet(
                CSCStripError = cms.double(7.0),
                CSCStripTimeOffset = cms.double(0.0),
                CSCTimeOffset = cms.double(0.0),
                CSCWireError = cms.double(8.6),
                CSCWireTimeOffset = cms.double(0.0),
                CSCsegments = cms.InputTag("hltCscSegments"),
                MatchParameters = cms.PSet(
                    CSCsegments = cms.InputTag("hltCscSegments"),
                    DTradius = cms.double(0.01),
                    DTsegments = cms.InputTag("hltDt4DSegments"),
                    TightMatchCSC = cms.bool(True),
                    TightMatchDT = cms.bool(False)
                ),
                PruneCut = cms.double(100.0),
                ServiceParameters = cms.PSet(
                    Propagators = cms.untracked.vstring('hltESPFastSteppingHelixPropagatorAny'),
                    RPCLayers = cms.bool(True)
                ),
                UseStripTime = cms.bool(True),
                UseWireTime = cms.bool(True),
                debug = cms.bool(False)
            ),
            DTTimingParameters = cms.PSet(
                DTTimeOffset = cms.double(2.7),
                DTsegments = cms.InputTag("hltDt4DSegments"),
                DoWireCorr = cms.bool(False),
                DropTheta = cms.bool(True),
                HitError = cms.double(6.0),
                HitsMin = cms.int32(5),
                MatchParameters = cms.PSet(
                    CSCsegments = cms.InputTag("hltCscSegments"),
                    DTradius = cms.double(0.01),
                    DTsegments = cms.InputTag("hltDt4DSegments"),
                    TightMatchCSC = cms.bool(True),
                    TightMatchDT = cms.bool(False)
                ),
                PruneCut = cms.double(10000.0),
                RequireBothProjections = cms.bool(False),
                ServiceParameters = cms.PSet(
                    Propagators = cms.untracked.vstring('hltESPFastSteppingHelixPropagatorAny'),
                    RPCLayers = cms.bool(True)
                ),
                UseSegmentT0 = cms.bool(False),
                debug = cms.bool(False)
            ),
            EcalEnergyCut = cms.double(0.4),
            ErrorCSC = cms.double(7.4),
            ErrorDT = cms.double(6.0),
            ErrorEB = cms.double(2.085),
            ErrorEE = cms.double(6.95),
            UseCSC = cms.bool(True),
            UseDT = cms.bool(True),
            UseECAL = cms.bool(True)
        ),
        TrackAssociatorParameters = cms.PSet(
            CSCSegmentCollectionLabel = cms.InputTag("hltCscSegments"),
            CaloTowerCollectionLabel = cms.InputTag("Notused"),
            DTRecSegment4DCollectionLabel = cms.InputTag("hltDt4DSegments"),
            EBRecHitCollectionLabel = cms.InputTag("Notused"),
            EERecHitCollectionLabel = cms.InputTag("Notused"),
            GEMSegmentCollectionLabel = cms.InputTag("hltGemSegments"),
            HBHERecHitCollectionLabel = cms.InputTag("Notused"),
            HORecHitCollectionLabel = cms.InputTag("Notused"),
            ME0SegmentCollectionLabel = cms.InputTag("hltMe0Segments"),
            accountForTrajectoryChangeCalo = cms.bool(False),
            dREcal = cms.double(9999.0),
            dREcalPreselection = cms.double(0.05),
            dRHcal = cms.double(9999.0),
            dRHcalPreselection = cms.double(0.2),
            dRMuon = cms.double(9999.0),
            dRMuonPreselection = cms.double(0.2),
            dRPreshowerPreselection = cms.double(0.2),
            muonMaxDistanceSigmaX = cms.double(0.0),
            muonMaxDistanceSigmaY = cms.double(0.0),
            muonMaxDistanceX = cms.double(5.0),
            muonMaxDistanceY = cms.double(5.0),
            propagateAllDirections = cms.bool(True),
            trajectoryUncertaintyTolerance = cms.double(-1.0),
            truthMatch = cms.bool(False),
            useCalo = cms.bool(False),
            useEcal = cms.bool(False),
            useGEM = cms.bool(True),
            useHO = cms.bool(False),
            useHcal = cms.bool(False),
            useME0 = cms.bool(True),
            useMuon = cms.bool(True),
            usePreshower = cms.bool(False)
        ),
        TrackExtractorPSet = cms.PSet(
            BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
            BeamlineOption = cms.string('BeamSpotFromEvent'),
            Chi2Ndof_Max = cms.double(1e+64),
            Chi2Prob_Min = cms.double(-1.0),
            ComponentName = cms.string('TrackExtractor'),
            DR_Max = cms.double(1.0),
            DR_Veto = cms.double(0.01),
            Diff_r = cms.double(0.1),
            Diff_z = cms.double(0.2),
            NHits_Min = cms.uint32(0),
            Pt_Min = cms.double(-1.0),
            inputTrackCollection = cms.InputTag("hltPhase2L3MuonMerged")
        ),
        TrackerKinkFinderParameters = cms.PSet(
            diagonalOnly = cms.bool(False),
            usePosition = cms.bool(False)
        ),
        addExtraSoftMuons = cms.bool(False),
        arbitrateTrackerMuons = cms.bool(True),
        arbitrationCleanerOptions = cms.PSet(
            ClusterDPhi = cms.double(0.6),
            ClusterDTheta = cms.double(0.02),
            Clustering = cms.bool(True),
            ME1a = cms.bool(True),
            Overlap = cms.bool(True),
            OverlapDPhi = cms.double(0.0786),
            OverlapDTheta = cms.double(0.02)
        ),
        debugWithTruthMatching = cms.bool(False),
        ecalDepositName = cms.string('ecal'),
        fillCaloCompatibility = cms.bool(False),
        fillEnergy = cms.bool(False),
        fillGlobalTrackQuality = cms.bool(False),
        fillGlobalTrackRefits = cms.bool(False),
        fillIsolation = cms.bool(False),
        fillMatching = cms.bool(True),
        fillTrackerKink = cms.bool(False),
        globalTrackQualityInputTag = cms.InputTag("glbTrackQual"),
        hcalDepositName = cms.string('hcal'),
        hoDepositName = cms.string('ho'),
        # inputCollectionLabels = cms.VInputTag("hltPhase2L3MuonMerged"),
        # inputCollectionTypes = cms.vstring(
        #     'inner tracks'
        # ),
        inputCollectionLabels = cms.VInputTag("hltPhase2L3MuonMerged", "hltPhase2L3GlbMuon", "hltL2MuonsFromL1TkMuon:UpdatedAtVtx"),
        inputCollectionTypes = cms.vstring(
            'inner tracks', 
            'links', 
            'outer tracks'
        ),
        jetDepositName = cms.string('jets'),
        maxAbsDx = cms.double(3.0),
        maxAbsDy = cms.double(9999.0),
        maxAbsEta = cms.double(3.0),
        maxAbsPullX = cms.double(4.0),
        maxAbsPullY = cms.double(9999.0),
        minCaloCompatibility = cms.double(0.6),
        minNumberOfMatches = cms.int32(1),
        minP = cms.double(0.0),
        minPCaloMuon = cms.double(1000000000.0),
        minPt = cms.double(2.0),
        ptThresholdToFillCandidateP4WithGlobalFit = cms.double(200.0),
        runArbitrationCleaner = cms.bool(False),
        sigmaThresholdToFillCandidateP4WithGlobalFit = cms.double(2.0),
        trackDepositName = cms.string('tracker'),
        writeIsoDeposits = cms.bool(False)
    )

    process.hltPhase2L3Muons = cms.EDProducer("MuonIDFilterProducerForHLT",
        allowedTypeMask = cms.uint32(0),
        applyTriggerIdLoose = cms.bool(True),
        inputMuonCollection = cms.InputTag("hltPhase2L3MuonsNoID"),
        maxNormalizedChi2 = cms.double(9999.0),
        minNMuonHits = cms.int32(0),
        minNMuonStations = cms.int32(0),
        minNTrkLayers = cms.int32(0),
        minPixHits = cms.int32(0),
        minPixLayer = cms.int32(0),
        minPt = cms.double(0.0),
        minTrkHits = cms.int32(0),
        requiredTypeMask = cms.uint32(0),
        typeMuon = cms.uint32(0)
    )

    process.hltL3MuonsPhase2L3Links = cms.EDProducer("MuonLinksProducer",
        inputCollection = cms.InputTag("hltPhase2L3Muons")
    )

    process.hltPhase2L3MuonTracks = cms.EDProducer("HLTMuonTrackSelector",
        copyExtras = cms.untracked.bool(True),
        copyMVA = cms.bool(False),
        copyTrajectories = cms.untracked.bool(False),
        muon = cms.InputTag("hltPhase2L3Muons"),
        originalMVAVals = cms.InputTag("none"),
        track = cms.InputTag("hltPhase2L3MuonMerged")
    )

    process.hltPhase2L3MuonCandidates = cms.EDProducer( "L3MuonCandidateProducerFromMuons",
        InputObjects = cms.InputTag("hltPhase2L3Muons")
    )

    process.HLTPhase2L3IOFromL1TkMuonTkCandidateSequence = cms.Sequence(
        process.hltPhase2L3MuonPixelTracksFilter+
        process.hltPhase2L3MuonPixelTracksFitter+
        process.hltPhase2L3FromL1TkMuonPixelTracksTrackingRegions+
        process.hltPhase2L3FromL1TkMuonPixelLayerQuadruplets+
        process.hltPhase2L3FromL1TkMuonPixelTracksHitDoublets+
        process.hltPhase2L3FromL1TkMuonPixelTracksHitQuadruplets+
        process.hltPhase2L3FromL1TkMuonPixelTracks+
        process.hltPhase2L3FromL1TkMuonPixelVertices+
        process.hltPhase2L3FromL1TkMuonTrimmedPixelVertices+
        process.hltIter0Phase2L3FromL1TkMuonPixelSeedsFromPixelTracks+
        process.hltIter0Phase2L3FromL1TkMuonCkfTrackCandidates+
        process.hltIter0Phase2L3FromL1TkMuonCtfWithMaterialTracks+
        process.hltIter0Phase2L3FromL1TkMuonTrackCutClassifier+
        process.hltIter0Phase2L3FromL1TkMuonTrackSelectionHighPurity+
        process.hltIter2Phase2L3FromL1TkMuonClustersRefRemoval+
        process.hltIter2Phase2L3FromL1TkMuonMaskedMeasurementTrackerEvent+
        process.hltIter2Phase2L3FromL1TkMuonPixelLayerTriplets+
        process.hltIter2Phase2L3FromL1TkMuonPixelClusterCheck+
        process.hltIter2Phase2L3FromL1TkMuonPixelHitDoublets+
        process.hltIter2Phase2L3FromL1TkMuonPixelHitTriplets+
        process.hltIter2Phase2L3FromL1TkMuonPixelSeeds+
        process.hltIter2Phase2L3FromL1TkMuonCkfTrackCandidates+
        process.hltIter2Phase2L3FromL1TkMuonCtfWithMaterialTracks+
        process.hltIter2Phase2L3FromL1TkMuonTrackCutClassifier+
        process.hltIter2Phase2L3FromL1TkMuonTrackSelectionHighPurity+
        process.hltIter2Phase2L3FromL1TkMuonMerged
    )

    process.HLTPhase2L3MuonRecoNocandSequence = cms.Sequence(
        process.HLTPhase2L3OImuonTkCandidateSequence+
        process.HLTPhase2L3IOFromL1TkMuonTkCandidateSequence+
        process.hltPhase2L3MuonMerged+
        process.hltPhase2L3GlbMuon+
        process.hltPhase2L3MuonsNoID+
        process.hltPhase2L3Muons+
        process.hltL3MuonsPhase2L3Links+
        process.hltPhase2L3MuonTracks
    )

    process.HLTPhase2L3MuonRecoSequence = cms.Sequence(
        process.HLTPhase2L3MuonRecoNocandSequence+
        process.hltPhase2L3MuonCandidates
    )

    return process

def customizePhase2MuonHLTReconstruction(process, processName = "MYHLT"):

    # -- ESProducers -- #
    process = loadPhase2MuonHLTESProducers(process)

    # -- PSet -- #
    process = loadPhase2MuonHLTPSets(process)

    # -- Begin sequence -- #
    process.hltTriggerType = cms.EDFilter("HLTTriggerTypeFilter",
        SelectedTriggerType = cms.int32(1)
    )

    process.hltScalersRawToDigi = cms.EDProducer("ScalersRawToDigi",
        scalersInputTag = cms.InputTag("rawDataCollector")
    )

    # process.hltOnlineBeamSpot = cms.EDProducer("BeamSpotOnlineProducer",
    #     changeToCMSCoordinates = cms.bool(False),
    #     gtEvmLabel = cms.InputTag(""),
    #     maxRadius = cms.double(2.0),
    #     maxZ = cms.double(40.0),
    #     setSigmaZ = cms.double(0.0),
    #     src = cms.InputTag("hltScalersRawToDigi")
    # )

    process.HLTBeamSpot = cms.Sequence(
        process.hltScalersRawToDigi+
        process.offlineBeamSpot
    )

    process.HLTBeginSequence = cms.Sequence(
        process.hltTriggerType+
        process.HLTBeamSpot
    )

    # -- Muon Local Reconstruction -- #
    process = customizeMuonLocalReco(process)

    # -- Tracker Local Reconstruction -- #
    process = customizeTrackerLocalReco(process)

    # -- L2 Muon -- #
    process = customizeL2MuonReco(process)

    # -- OI -- #
    process = customizeOI(process)

    # -- IO -- #
    process = customizeIO(process)

    # -- L3 Muon -- #
    process = customizeL3Muon(process)

    # -- End sequence -- #
    process.hltBoolEnd = cms.EDFilter("HLTBool",
        result = cms.bool(True)
    )

    process.HLTEndSequence = cms.Sequence(
        process.hltBoolEnd
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

    return process

def customizePhase2MuonHLTCaloLocalReco(process, processName = "MYHLT"):

    process.bunchSpacingProducer = cms.EDProducer("BunchSpacingProducer")

    process.hltEcalDigis = cms.EDProducer("EcalRawToDigi",
        DoRegional = cms.bool(False),
        FEDs = cms.vint32(
            601, 602, 603, 604, 605, 
            606, 607, 608, 609, 610, 
            611, 612, 613, 614, 615, 
            616, 617, 618, 619, 620, 
            621, 622, 623, 624, 625, 
            626, 627, 628, 629, 630, 
            631, 632, 633, 634, 635, 
            636, 637, 638, 639, 640, 
            641, 642, 643, 644, 645, 
            646, 647, 648, 649, 650, 
            651, 652, 653, 654
        ),
        FedLabel = cms.InputTag("listfeds"),
        InputLabel = cms.InputTag("rawDataCollector"),
        eventPut = cms.bool(True),
        feIdCheck = cms.bool(True),
        feUnpacking = cms.bool(True),
        forceToKeepFRData = cms.bool(False),
        headerUnpacking = cms.bool(True),
        memUnpacking = cms.bool(True),
        mightGet = cms.optional.untracked.vstring,
        numbTriggerTSamples = cms.int32(1),
        numbXtalTSamples = cms.int32(10),
        orderedDCCIdList = cms.vint32(
            1, 2, 3, 4, 5, 
            6, 7, 8, 9, 10, 
            11, 12, 13, 14, 15, 
            16, 17, 18, 19, 20, 
            21, 22, 23, 24, 25, 
            26, 27, 28, 29, 30, 
            31, 32, 33, 34, 35, 
            36, 37, 38, 39, 40, 
            41, 42, 43, 44, 45, 
            46, 47, 48, 49, 50, 
            51, 52, 53, 54
        ),
        orderedFedList = cms.vint32(
            601, 602, 603, 604, 605, 
            606, 607, 608, 609, 610, 
            611, 612, 613, 614, 615, 
            616, 617, 618, 619, 620, 
            621, 622, 623, 624, 625, 
            626, 627, 628, 629, 630, 
            631, 632, 633, 634, 635, 
            636, 637, 638, 639, 640, 
            641, 642, 643, 644, 645, 
            646, 647, 648, 649, 650, 
            651, 652, 653, 654
        ),
        silentMode = cms.untracked.bool(True),
        srpUnpacking = cms.bool(True),
        syncCheck = cms.bool(True),
        tccUnpacking = cms.bool(True)
    )

    process.hltEcalPreshowerDigis = cms.EDProducer("ESRawToDigi",
        ESdigiCollection = cms.string(''),
        InstanceES = cms.string(''),
        LookupTable = cms.FileInPath('EventFilter/ESDigiToRaw/data/ES_lookup_table.dat'),
        debugMode = cms.untracked.bool(False),
        mightGet = cms.optional.untracked.vstring,
        sourceTag = cms.InputTag("rawDataCollector")
    )

    process.hltEcalUncalibRecHit = cms.EDProducer("EcalUncalibRecHitProducer",
        EBdigiCollection = cms.InputTag("hltEcalDigis","ebDigis"),
        EBhitCollection = cms.string('EcalUncalibRecHitsEB'),
        EEdigiCollection = cms.InputTag("hltEcalDigis","eeDigis"),
        EEhitCollection = cms.string('EcalUncalibRecHitsEE'),
        algo = cms.string('EcalUncalibRecHitWorkerMultiFit'),
        algoPSet = cms.PSet(
            EBamplitudeFitParameters = cms.vdouble(1.138, 1.652),
            EBtimeConstantTerm = cms.double(0.6),
            EBtimeFitLimits_Lower = cms.double(0.2),
            EBtimeFitLimits_Upper = cms.double(1.4),
            EBtimeFitParameters = cms.vdouble(
                -2.015452, 3.130702, -12.3473, 41.88921, -82.83944, 
                91.01147, -50.35761, 11.05621
            ),
            EBtimeNconst = cms.double(28.5),
            EEamplitudeFitParameters = cms.vdouble(1.89, 1.4),
            EEtimeConstantTerm = cms.double(1.0),
            EEtimeFitLimits_Lower = cms.double(0.2),
            EEtimeFitLimits_Upper = cms.double(1.4),
            EEtimeFitParameters = cms.vdouble(
                -2.390548, 3.553628, -17.62341, 67.67538, -133.213, 
                140.7432, -75.41106, 16.20277
            ),
            EEtimeNconst = cms.double(31.8),
            EcalPulseShapeParameters = cms.PSet(
                EBCorrNoiseMatrixG01 = cms.vdouble(
                    1.0, 0.73354, 0.64442, 0.58851, 0.55425, 
                    0.53082, 0.51916, 0.51097, 0.50732, 0.50409
                ),
                EBCorrNoiseMatrixG06 = cms.vdouble(
                    1.0, 0.70946, 0.58021, 0.49846, 0.45006, 
                    0.41366, 0.39699, 0.38478, 0.37847, 0.37055
                ),
                EBCorrNoiseMatrixG12 = cms.vdouble(
                    1.0, 0.71073, 0.55721, 0.46089, 0.40449, 
                    0.35931, 0.33924, 0.32439, 0.31581, 0.30481
                ),
                EBPulseShapeCovariance = cms.vdouble(
                    3.001e-06, 1.233e-05, 0.0, -4.416e-06, -4.571e-06, 
                    -3.614e-06, -2.636e-06, -1.286e-06, -8.41e-07, -5.296e-07, 
                    0.0, 0.0, 1.233e-05, 6.154e-05, 0.0, 
                    -2.2e-05, -2.309e-05, -1.838e-05, -1.373e-05, -7.334e-06, 
                    -5.088e-06, -3.745e-06, -2.428e-06, 0.0, 0.0, 
                    0.0, 0.0, 0.0, 0.0, 0.0, 
                    0.0, 0.0, 0.0, 0.0, 0.0, 
                    0.0, -4.416e-06, -2.2e-05, 0.0, 8.319e-06, 
                    8.545e-06, 6.792e-06, 5.059e-06, 2.678e-06, 1.816e-06, 
                    1.223e-06, 8.245e-07, 5.589e-07, -4.571e-06, -2.309e-05, 
                    0.0, 8.545e-06, 9.182e-06, 7.219e-06, 5.388e-06, 
                    2.853e-06, 1.944e-06, 1.324e-06, 9.083e-07, 6.335e-07, 
                    -3.614e-06, -1.838e-05, 0.0, 6.792e-06, 7.219e-06, 
                    6.016e-06, 4.437e-06, 2.385e-06, 1.636e-06, 1.118e-06, 
                    7.754e-07, 5.556e-07, -2.636e-06, -1.373e-05, 0.0, 
                    5.059e-06, 5.388e-06, 4.437e-06, 3.602e-06, 1.917e-06, 
                    1.322e-06, 9.079e-07, 6.529e-07, 4.752e-07, -1.286e-06, 
                    -7.334e-06, 0.0, 2.678e-06, 2.853e-06, 2.385e-06, 
                    1.917e-06, 1.375e-06, 9.1e-07, 6.455e-07, 4.693e-07, 
                    3.657e-07, -8.41e-07, -5.088e-06, 0.0, 1.816e-06, 
                    1.944e-06, 1.636e-06, 1.322e-06, 9.1e-07, 9.115e-07, 
                    6.062e-07, 4.436e-07, 3.422e-07, -5.296e-07, -3.745e-06, 
                    0.0, 1.223e-06, 1.324e-06, 1.118e-06, 9.079e-07, 
                    6.455e-07, 6.062e-07, 7.217e-07, 4.862e-07, 3.768e-07, 
                    0.0, -2.428e-06, 0.0, 8.245e-07, 9.083e-07, 
                    7.754e-07, 6.529e-07, 4.693e-07, 4.436e-07, 4.862e-07, 
                    6.509e-07, 4.418e-07, 0.0, 0.0, 0.0, 
                    5.589e-07, 6.335e-07, 5.556e-07, 4.752e-07, 3.657e-07, 
                    3.422e-07, 3.768e-07, 4.418e-07, 6.142e-07
                ),
                EBPulseShapeTemplate = cms.vdouble(
                    0.0113979, 0.758151, 1.0, 0.887744, 0.673548, 
                    0.474332, 0.319561, 0.215144, 0.147464, 0.101087, 
                    0.0693181, 0.0475044
                ),
                EBdigiCollection = cms.string(''),
                EECorrNoiseMatrixG01 = cms.vdouble(
                    1.0, 0.72698, 0.62048, 0.55691, 0.51848, 
                    0.49147, 0.47813, 0.47007, 0.46621, 0.46265
                ),
                EECorrNoiseMatrixG06 = cms.vdouble(
                    1.0, 0.71217, 0.47464, 0.34056, 0.26282, 
                    0.20287, 0.17734, 0.16256, 0.15618, 0.14443
                ),
                EECorrNoiseMatrixG12 = cms.vdouble(
                    1.0, 0.71373, 0.44825, 0.30152, 0.21609, 
                    0.14786, 0.11772, 0.10165, 0.09465, 0.08098
                ),
                EEPulseShapeCovariance = cms.vdouble(
                    3.941e-05, 3.333e-05, 0.0, -1.449e-05, -1.661e-05, 
                    -1.424e-05, -1.183e-05, -6.842e-06, -4.915e-06, -3.411e-06, 
                    0.0, 0.0, 3.333e-05, 2.862e-05, 0.0, 
                    -1.244e-05, -1.431e-05, -1.233e-05, -1.032e-05, -5.883e-06, 
                    -4.154e-06, -2.902e-06, -2.128e-06, 0.0, 0.0, 
                    0.0, 0.0, 0.0, 0.0, 0.0, 
                    0.0, 0.0, 0.0, 0.0, 0.0, 
                    0.0, -1.449e-05, -1.244e-05, 0.0, 5.84e-06, 
                    6.649e-06, 5.72e-06, 4.812e-06, 2.708e-06, 1.869e-06, 
                    1.33e-06, 9.186e-07, 6.446e-07, -1.661e-05, -1.431e-05, 
                    0.0, 6.649e-06, 7.966e-06, 6.898e-06, 5.794e-06, 
                    3.157e-06, 2.184e-06, 1.567e-06, 1.084e-06, 7.575e-07, 
                    -1.424e-05, -1.233e-05, 0.0, 5.72e-06, 6.898e-06, 
                    6.341e-06, 5.347e-06, 2.859e-06, 1.991e-06, 1.431e-06, 
                    9.839e-07, 6.886e-07, -1.183e-05, -1.032e-05, 0.0, 
                    4.812e-06, 5.794e-06, 5.347e-06, 4.854e-06, 2.628e-06, 
                    1.809e-06, 1.289e-06, 9.02e-07, 6.146e-07, -6.842e-06, 
                    -5.883e-06, 0.0, 2.708e-06, 3.157e-06, 2.859e-06, 
                    2.628e-06, 1.863e-06, 1.296e-06, 8.882e-07, 6.108e-07, 
                    4.283e-07, -4.915e-06, -4.154e-06, 0.0, 1.869e-06, 
                    2.184e-06, 1.991e-06, 1.809e-06, 1.296e-06, 1.217e-06, 
                    8.669e-07, 5.751e-07, 3.882e-07, -3.411e-06, -2.902e-06, 
                    0.0, 1.33e-06, 1.567e-06, 1.431e-06, 1.289e-06, 
                    8.882e-07, 8.669e-07, 9.522e-07, 6.717e-07, 4.293e-07, 
                    0.0, -2.128e-06, 0.0, 9.186e-07, 1.084e-06, 
                    9.839e-07, 9.02e-07, 6.108e-07, 5.751e-07, 6.717e-07, 
                    7.911e-07, 5.493e-07, 0.0, 0.0, 0.0, 
                    6.446e-07, 7.575e-07, 6.886e-07, 6.146e-07, 4.283e-07, 
                    3.882e-07, 4.293e-07, 5.493e-07, 7.027e-07
                ),
                EEPulseShapeTemplate = cms.vdouble(
                    0.116442, 0.756246, 1.0, 0.897182, 0.686831, 
                    0.491506, 0.344111, 0.245731, 0.174115, 0.123361, 
                    0.0874288, 0.061957
                ),
                EEdigiCollection = cms.string(''),
                ESdigiCollection = cms.string(''),
                EcalPreMixStage1 = cms.bool(False),
                EcalPreMixStage2 = cms.bool(False),
                UseLCcorrection = cms.untracked.bool(True)
            ),
            activeBXs = cms.vint32(
                -5, -4, -3, -2, -1, 
                0, 1, 2, 3, 4
            ),
            addPedestalUncertaintyEB = cms.double(0.0),
            addPedestalUncertaintyEE = cms.double(0.0),
            ampErrorCalculation = cms.bool(True),
            amplitudeThresholdEB = cms.double(10),
            amplitudeThresholdEE = cms.double(10),
            chi2ThreshEB_ = cms.double(65.0),
            chi2ThreshEE_ = cms.double(50.0),
            doPrefitEB = cms.bool(False),
            doPrefitEE = cms.bool(False),
            dynamicPedestalsEB = cms.bool(False),
            dynamicPedestalsEE = cms.bool(False),
            ebPulseShape = cms.vdouble(
                5.2e-05, -5.26e-05, 6.66e-05, 0.1168, 0.7575, 
                1.0, 0.8876, 0.6732, 0.4741, 0.3194
            ),
            ebSpikeThreshold = cms.double(1.042),
            eePulseShape = cms.vdouble(
                5.2e-05, -5.26e-05, 6.66e-05, 0.1168, 0.7575, 
                1.0, 0.8876, 0.6732, 0.4741, 0.3194
            ),
            gainSwitchUseMaxSampleEB = cms.bool(True),
            gainSwitchUseMaxSampleEE = cms.bool(False),
            kPoorRecoFlagEB = cms.bool(True),
            kPoorRecoFlagEE = cms.bool(False),
            mitigateBadSamplesEB = cms.bool(False),
            mitigateBadSamplesEE = cms.bool(False),
            outOfTimeThresholdGain12mEB = cms.double(5),
            outOfTimeThresholdGain12mEE = cms.double(1000),
            outOfTimeThresholdGain12pEB = cms.double(5),
            outOfTimeThresholdGain12pEE = cms.double(1000),
            outOfTimeThresholdGain61mEB = cms.double(5),
            outOfTimeThresholdGain61mEE = cms.double(1000),
            outOfTimeThresholdGain61pEB = cms.double(5),
            outOfTimeThresholdGain61pEE = cms.double(1000),
            prefitMaxChiSqEB = cms.double(25.0),
            prefitMaxChiSqEE = cms.double(10.0),
            selectiveBadSampleCriteriaEB = cms.bool(False),
            selectiveBadSampleCriteriaEE = cms.bool(False),
            simplifiedNoiseModelForGainSwitch = cms.bool(True),
            timealgo = cms.string('None'),
            useLumiInfoRunHeader = cms.bool(False)
        )
    )

    process.hltEcalDetIdToBeRecovered = cms.EDProducer("EcalDetIdToBeRecoveredProducer",
        ebDetIdToBeRecovered = cms.string('ebDetId'),
        ebFEToBeRecovered = cms.string('ebFE'),
        ebIntegrityChIdErrors = cms.InputTag("hltEcalDigis","EcalIntegrityChIdErrors"),
        ebIntegrityGainErrors = cms.InputTag("hltEcalDigis","EcalIntegrityGainErrors"),
        ebIntegrityGainSwitchErrors = cms.InputTag("hltEcalDigis","EcalIntegrityGainSwitchErrors"),
        ebSrFlagCollection = cms.InputTag("hltEcalDigis"),
        eeDetIdToBeRecovered = cms.string('eeDetId'),
        eeFEToBeRecovered = cms.string('eeFE'),
        eeIntegrityChIdErrors = cms.InputTag("hltEcalDigis","EcalIntegrityChIdErrors"),
        eeIntegrityGainErrors = cms.InputTag("hltEcalDigis","EcalIntegrityGainErrors"),
        eeIntegrityGainSwitchErrors = cms.InputTag("hltEcalDigis","EcalIntegrityGainSwitchErrors"),
        eeSrFlagCollection = cms.InputTag("hltEcalDigis"),
        integrityBlockSizeErrors = cms.InputTag("hltEcalDigis","EcalIntegrityBlockSizeErrors"),
        integrityTTIdErrors = cms.InputTag("hltEcalDigis","EcalIntegrityTTIdErrors")
    )

    process.hltEcalRecHit = cms.EDProducer("EcalRecHitProducer",
        ChannelStatusToBeExcluded = cms.vstring(
            'kDAC', 
            'kNoisy', 
            'kNNoisy', 
            'kFixedG6', 
            'kFixedG1', 
            'kFixedG0', 
            'kNonRespondingIsolated', 
            'kDeadVFE', 
            'kDeadFE', 
            'kNoDataNoTP'
        ),
        EBLaserMAX = cms.double(3.0),
        EBLaserMIN = cms.double(0.5),
        EBrechitCollection = cms.string('EcalRecHitsEB'),
        EBuncalibRecHitCollection = cms.InputTag("hltEcalUncalibRecHit","EcalUncalibRecHitsEB"),
        EELaserMAX = cms.double(8.0),
        EELaserMIN = cms.double(0.5),
        EErechitCollection = cms.string('EcalRecHitsEE'),
        EEuncalibRecHitCollection = cms.InputTag("hltEcalUncalibRecHit","EcalUncalibRecHitsEE"),
        algo = cms.string('EcalRecHitWorkerSimple'),
        algoRecover = cms.string('EcalRecHitWorkerRecover'),
        bdtWeightFileCracks = cms.FileInPath('RecoLocalCalo/EcalDeadChannelRecoveryAlgos/data/BDTWeights/bdtgAllRH_8GT700MeV_onlyCracks_ZskimData2017_v1.xml'),
        bdtWeightFileNoCracks = cms.FileInPath('RecoLocalCalo/EcalDeadChannelRecoveryAlgos/data/BDTWeights/bdtgAllRH_8GT700MeV_noCracks_ZskimData2017_v1.xml'),
        cleaningConfig = cms.PSet(
            cThreshold_barrel = cms.double(4),
            cThreshold_double = cms.double(10),
            cThreshold_endcap = cms.double(15),
            e4e1Threshold_barrel = cms.double(0.08),
            e4e1Threshold_endcap = cms.double(0.3),
            e4e1_a_barrel = cms.double(0.02),
            e4e1_a_endcap = cms.double(0.02),
            e4e1_b_barrel = cms.double(0.02),
            e4e1_b_endcap = cms.double(-0.0125),
            e6e2thresh = cms.double(0.04),
            ignoreOutOfTimeThresh = cms.double(1000000000.0),
            tightenCrack_e1_double = cms.double(2),
            tightenCrack_e1_single = cms.double(1),
            tightenCrack_e4e1_single = cms.double(2.5),
            tightenCrack_e6e2_double = cms.double(3)
        ),
        dbStatusToBeExcludedEB = cms.vint32(14, 78, 142),
        dbStatusToBeExcludedEE = cms.vint32(14, 78, 142),
        ebDetIdToBeRecovered = cms.InputTag("hltEcalDetIdToBeRecovered","ebDetId"),
        ebFEToBeRecovered = cms.InputTag("hltEcalDetIdToBeRecovered","ebFE"),
        eeDetIdToBeRecovered = cms.InputTag("hltEcalDetIdToBeRecovered","eeDetId"),
        eeFEToBeRecovered = cms.InputTag("hltEcalDetIdToBeRecovered","eeFE"),
        flagsMapDBReco = cms.PSet(
            kDead = cms.vstring('kNoDataNoTP'),
            kGood = cms.vstring(
                'kOk', 
                'kDAC', 
                'kNoLaser', 
                'kNoisy'
            ),
            kNeighboursRecovered = cms.vstring(
                'kFixedG0', 
                'kNonRespondingIsolated', 
                'kDeadVFE'
            ),
            kNoisy = cms.vstring(
                'kNNoisy', 
                'kFixedG6', 
                'kFixedG1'
            ),
            kTowerRecovered = cms.vstring('kDeadFE')
        ),
        killDeadChannels = cms.bool(True),
        laserCorrection = cms.bool(True),
        logWarningEtThreshold_EB_FE = cms.double(50),
        logWarningEtThreshold_EE_FE = cms.double(50),
        recoverEBFE = cms.bool(True),
        recoverEBIsolatedChannels = cms.bool(False),
        recoverEBVFE = cms.bool(False),
        recoverEEFE = cms.bool(True),
        recoverEEIsolatedChannels = cms.bool(False),
        recoverEEVFE = cms.bool(False),
        singleChannelRecoveryMethod = cms.string('BDTG'),
        singleChannelRecoveryThreshold = cms.double(0.7),
        skipTimeCalib = cms.bool(False),
        sum8ChannelRecoveryThreshold = cms.double(0.0),
        triggerPrimitiveDigiCollection = cms.InputTag("hltEcalDigis","EcalTriggerPrimitives")
    )

    process.hltEcalPreshowerRecHit = cms.EDProducer("ESRecHitProducer",
        ESRecoAlgo = cms.int32(0),
        ESdigiCollection = cms.InputTag("hltEcalPreshowerDigis"),
        ESrechitCollection = cms.string('EcalRecHitsES'),
        algo = cms.string('ESRecHitWorker')
    )

    process.HLTDoFullUnpackingEgammaEcalSequence = cms.Sequence(
        process.bunchSpacingProducer + 
        process.hltEcalDigis + 
        process.hltEcalPreshowerDigis + 
        process.hltEcalUncalibRecHit + 
        process.hltEcalDetIdToBeRecovered + 
        process.hltEcalRecHit + 
        process.hltEcalPreshowerRecHit
    )

    process.hltHcalDigis = cms.EDProducer("HcalRawToDigi",
        ComplainEmptyData = cms.untracked.bool(False),
        ElectronicsMap = cms.string(''),
        ExpectedOrbitMessageTime = cms.untracked.int32(-1),
        FEDs = cms.untracked.vint32(),
        FilterDataQuality = cms.bool(True),
        HcalFirstFED = cms.untracked.int32(700),
        InputLabel = cms.InputTag("rawDataCollector"),
        UnpackCalib = cms.untracked.bool(True),
        UnpackTTP = cms.untracked.bool(True),
        UnpackUMNio = cms.untracked.bool(True),
        UnpackZDC = cms.untracked.bool(True),
        UnpackerMode = cms.untracked.int32(0),
        firstSample = cms.int32(0),
        lastSample = cms.int32(9),
        mightGet = cms.optional.untracked.vstring,
        saveQIE10DataNSamples = cms.untracked.vint32(),
        saveQIE10DataTags = cms.untracked.vstring(),
        saveQIE11DataNSamples = cms.untracked.vint32(),
        saveQIE11DataTags = cms.untracked.vstring(),
        silent = cms.untracked.bool(True)
    )

    process.hltHbhereco = cms.EDProducer("HBHEPhase1Reconstructor",
        algoConfigClass = cms.string(''),
        algorithm = cms.PSet(
            Class = cms.string('SimpleHBHEPhase1Algo'),
            activeBXs = cms.vint32(
                -3, -2, -1, 0, 1, 
                2, 3, 4
            ),
            applyLegacyHBMCorrection = cms.bool(False),
            applyPedConstraint = cms.bool(True),
            applyPulseJitter = cms.bool(False),
            applyTimeConstraint = cms.bool(True),
            applyTimeSlew = cms.bool(True),
            applyTimeSlewM3 = cms.bool(True),
            calculateArrivalTime = cms.bool(True),
            chiSqSwitch = cms.double(15.0),
            correctForPhaseContainment = cms.bool(True),
            correctionPhaseNS = cms.double(6.0),
            deltaChiSqThresh = cms.double(0.001),
            dynamicPed = cms.bool(False),
            firstSampleShift = cms.int32(0),
            fitTimes = cms.int32(1),
            meanPed = cms.double(0.0),
            meanTime = cms.double(0.0),
            nMaxItersMin = cms.int32(500),
            nMaxItersNNLS = cms.int32(500),
            nnlsThresh = cms.double(1e-11),
            pulseJitter = cms.double(1.0),
            respCorrM3 = cms.double(1.0),
            samplesToAdd = cms.int32(2),
            tdcTimeShift = cms.double(0.0),
            timeMax = cms.double(12.5),
            timeMin = cms.double(-12.5),
            timeSigmaHPD = cms.double(5.0),
            timeSigmaSiPM = cms.double(2.5),
            timeSlewParsType = cms.int32(3),
            ts4Max = cms.vdouble(100.0, 20000.0, 30000.0),
            ts4Min = cms.double(0.0),
            ts4Thresh = cms.double(0.0),
            ts4chi2 = cms.vdouble(15.0, 15.0),
            useM2 = cms.bool(False),
            useM3 = cms.bool(True),
            useMahi = cms.bool(True)
        ),
        digiLabelQIE11 = cms.InputTag("hltHcalDigis"),
        digiLabelQIE8 = cms.InputTag("hltHcalDigis"),
        dropZSmarkedPassed = cms.bool(True),
        flagParametersQIE11 = cms.PSet(

        ),
        flagParametersQIE8 = cms.PSet(
            hitEnergyMinimum = cms.double(1.0),
            hitMultiplicityThreshold = cms.int32(17),
            nominalPedestal = cms.double(3.0),
            pulseShapeParameterSets = cms.VPSet(
                cms.PSet(
                    pulseShapeParameters = cms.vdouble(
                        0.0, 100.0, -50.0, 0.0, -15.0, 
                        0.15
                    )
                ), 
                cms.PSet(
                    pulseShapeParameters = cms.vdouble(
                        100.0, 2000.0, -50.0, 0.0, -5.0, 
                        0.05
                    )
                ), 
                cms.PSet(
                    pulseShapeParameters = cms.vdouble(
                        2000.0, 1000000.0, -50.0, 0.0, 95.0, 
                        0.0
                    )
                ), 
                cms.PSet(
                    pulseShapeParameters = cms.vdouble(
                        -1000000.0, 1000000.0, 45.0, 0.1, 1000000.0, 
                        0.0
                    )
                )
            )
        ),
        makeRecHits = cms.bool(True),
        processQIE11 = cms.bool(True),
        processQIE8 = cms.bool(True),
        pulseShapeParametersQIE11 = cms.PSet(

        ),
        pulseShapeParametersQIE8 = cms.PSet(
            LeftSlopeCut = cms.vdouble(5.0, 2.55, 2.55),
            LeftSlopeThreshold = cms.vdouble(250.0, 500.0, 100000.0),
            LinearCut = cms.vdouble(-3.0, -0.054, -0.054),
            LinearThreshold = cms.vdouble(20.0, 100.0, 100000.0),
            MinimumChargeThreshold = cms.double(20.0),
            MinimumTS4TS5Threshold = cms.double(100.0),
            R45MinusOneRange = cms.double(0.2),
            R45PlusOneRange = cms.double(0.2),
            RMS8MaxCut = cms.vdouble(-13.5, -11.5, -11.5),
            RMS8MaxThreshold = cms.vdouble(20.0, 100.0, 100000.0),
            RightSlopeCut = cms.vdouble(5.0, 4.15, 4.15),
            RightSlopeSmallCut = cms.vdouble(1.08, 1.16, 1.16),
            RightSlopeSmallThreshold = cms.vdouble(150.0, 200.0, 100000.0),
            RightSlopeThreshold = cms.vdouble(250.0, 400.0, 100000.0),
            TS3TS4ChargeThreshold = cms.double(70.0),
            TS3TS4UpperChargeThreshold = cms.double(20.0),
            TS4TS5ChargeThreshold = cms.double(70.0),
            TS4TS5LowerCut = cms.vdouble(
                -1.0, -0.7, -0.5, -0.4, -0.3, 
                0.1
            ),
            TS4TS5LowerThreshold = cms.vdouble(
                100.0, 120.0, 160.0, 200.0, 300.0, 
                500.0
            ),
            TS4TS5UpperCut = cms.vdouble(1.0, 0.8, 0.75, 0.72),
            TS4TS5UpperThreshold = cms.vdouble(70.0, 90.0, 100.0, 400.0),
            TS5TS6ChargeThreshold = cms.double(70.0),
            TS5TS6UpperChargeThreshold = cms.double(20.0),
            TriangleIgnoreSlow = cms.bool(False),
            TrianglePeakTS = cms.uint32(10000),
            UseDualFit = cms.bool(True)
        ),
        recoParamsFromDB = cms.bool(True),
        saveDroppedInfos = cms.bool(False),
        saveEffectivePedestal = cms.bool(True),
        saveInfos = cms.bool(False),
        setLegacyFlagsQIE11 = cms.bool(False),
        setLegacyFlagsQIE8 = cms.bool(True),
        setNegativeFlagsQIE11 = cms.bool(False),
        setNegativeFlagsQIE8 = cms.bool(True),
        setNoiseFlagsQIE11 = cms.bool(False),
        setNoiseFlagsQIE8 = cms.bool(True),
        setPulseShapeFlagsQIE11 = cms.bool(False),
        setPulseShapeFlagsQIE8 = cms.bool(True),
        sipmQNTStoSum = cms.int32(3),
        sipmQTSShift = cms.int32(0),
        tsFromDB = cms.bool(False),
        use8ts = cms.bool(True)
    )

    process.hltHfprereco = cms.EDProducer("HFPreReconstructor",
        digiLabel = cms.InputTag("hltHcalDigis"),
        dropZSmarkedPassed = cms.bool(True),
        forceSOI = cms.int32(-1),
        soiShift = cms.int32(0),
        sumAllTimeSlices = cms.bool(False),
        tsFromDB = cms.bool(False)
    )

    process.hltHfreco = cms.EDProducer("HFPhase1Reconstructor",
        HFStripFilter = cms.PSet(
            gap = cms.int32(2),
            lstrips = cms.int32(2),
            maxStripTime = cms.double(10.0),
            maxThreshold = cms.double(100.0),
            seedHitIetaMax = cms.int32(35),
            stripThreshold = cms.double(40.0),
            timeMax = cms.double(6.0),
            verboseLevel = cms.untracked.int32(10),
            wedgeCut = cms.double(0.05)
        ),
        PETstat = cms.PSet(
            HcalAcceptSeverityLevel = cms.int32(9),
            longETParams = cms.vdouble(
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0
            ),
            longEnergyParams = cms.vdouble(
                43.5, 45.7, 48.32, 51.36, 54.82, 
                58.7, 63.0, 67.72, 72.86, 78.42, 
                84.4, 90.8, 97.62
            ),
            long_R = cms.vdouble(0.98),
            long_R_29 = cms.vdouble(0.8),
            shortETParams = cms.vdouble(
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0
            ),
            shortEnergyParams = cms.vdouble(
                35.1773, 35.37, 35.7933, 36.4472, 37.3317, 
                38.4468, 39.7925, 41.3688, 43.1757, 45.2132, 
                47.4813, 49.98, 52.7093
            ),
            short_R = cms.vdouble(0.8),
            short_R_29 = cms.vdouble(0.8)
        ),
        S8S1stat = cms.PSet(
            HcalAcceptSeverityLevel = cms.int32(9),
            isS8S1 = cms.bool(True),
            longETParams = cms.vdouble(
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0
            ),
            longEnergyParams = cms.vdouble(
                40.0, 100.0, 100.0, 100.0, 100.0, 
                100.0, 100.0, 100.0, 100.0, 100.0, 
                100.0, 100.0, 100.0
            ),
            long_optimumSlope = cms.vdouble(
                0.3, 0.1, 0.1, 0.1, 0.1, 
                0.1, 0.1, 0.1, 0.1, 0.1, 
                0.1, 0.1, 0.1
            ),
            shortETParams = cms.vdouble(
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0
            ),
            shortEnergyParams = cms.vdouble(
                40.0, 100.0, 100.0, 100.0, 100.0, 
                100.0, 100.0, 100.0, 100.0, 100.0, 
                100.0, 100.0, 100.0
            ),
            short_optimumSlope = cms.vdouble(
                0.3, 0.1, 0.1, 0.1, 0.1, 
                0.1, 0.1, 0.1, 0.1, 0.1, 
                0.1, 0.1, 0.1
            )
        ),
        S9S1stat = cms.PSet(
            HcalAcceptSeverityLevel = cms.int32(9),
            isS8S1 = cms.bool(False),
            longETParams = cms.vdouble(
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0
            ),
            longEnergyParams = cms.vdouble(
                43.5, 45.7, 48.32, 51.36, 54.82, 
                58.7, 63.0, 67.72, 72.86, 78.42, 
                84.4, 90.8, 97.62
            ),
            long_optimumSlope = cms.vdouble(
                -99999.0, 0.0164905, 0.0238698, 0.0321383, 0.041296, 
                0.0513428, 0.0622789, 0.0741041, 0.0868186, 0.100422, 
                0.135313, 0.136289, 0.0589927
            ),
            shortETParams = cms.vdouble(
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0, 0.0, 0.0, 
                0.0, 0.0, 0.0
            ),
            shortEnergyParams = cms.vdouble(
                35.1773, 35.37, 35.7933, 36.4472, 37.3317, 
                38.4468, 39.7925, 41.3688, 43.1757, 45.2132, 
                47.4813, 49.98, 52.7093
            ),
            short_optimumSlope = cms.vdouble(
                -99999.0, 0.0164905, 0.0238698, 0.0321383, 0.041296, 
                0.0513428, 0.0622789, 0.0741041, 0.0868186, 0.100422, 
                0.135313, 0.136289, 0.0589927
            )
        ),
        algoConfigClass = cms.string('HFPhase1PMTParams'),
        algorithm = cms.PSet(
            Class = cms.string('HFFlexibleTimeCheck'),
            energyWeights = cms.vdouble(
                1.0, 1.0, 1.0, 0.0, 1.0, 
                0.0, 2.0, 0.0, 2.0, 0.0, 
                2.0, 0.0, 1.0, 0.0, 0.0, 
                1.0, 0.0, 1.0, 0.0, 2.0, 
                0.0, 2.0, 0.0, 2.0, 0.0, 
                1.0
            ),
            rejectAllFailures = cms.bool(True),
            soiPhase = cms.uint32(1),
            tfallIfNoTDC = cms.double(-101.0),
            timeShift = cms.double(0.0),
            tlimits = cms.vdouble(-1000.0, 1000.0, -1000.0, 1000.0),
            triseIfNoTDC = cms.double(-100.0)
        ),
        checkChannelQualityForDepth3and4 = cms.bool(False),
        inputLabel = cms.InputTag("hltHfprereco"),
        runHFStripFilter = cms.bool(False),
        setNoiseFlags = cms.bool(True),
        useChannelQualityFromDB = cms.bool(False)
    )

    process.hltHoreco = cms.EDProducer("HcalHitReconstructor",
        HFInWindowStat = cms.PSet(

        ),
        PETstat = cms.PSet(

        ),
        S8S1stat = cms.PSet(

        ),
        S9S1stat = cms.PSet(

        ),
        Subdetector = cms.string('HO'),
        correctForPhaseContainment = cms.bool(True),
        correctForTimeslew = cms.bool(True),
        correctTiming = cms.bool(False),
        correctionPhaseNS = cms.double(13.0),
        dataOOTCorrectionCategory = cms.string('Data'),
        dataOOTCorrectionName = cms.string(''),
        digiLabel = cms.InputTag("hltHcalDigis"),
        digiTimeFromDB = cms.bool(True),
        digistat = cms.PSet(

        ),
        dropZSmarkedPassed = cms.bool(True),
        firstAuxTS = cms.int32(4),
        firstSample = cms.int32(4),
        hfTimingTrustParameters = cms.PSet(

        ),
        mcOOTCorrectionCategory = cms.string('MC'),
        mcOOTCorrectionName = cms.string(''),
        recoParamsFromDB = cms.bool(True),
        samplesToAdd = cms.int32(4),
        saturationParameters = cms.PSet(
            maxADCvalue = cms.int32(127)
        ),
        setHSCPFlags = cms.bool(False),
        setNegativeFlags = cms.bool(False),
        setNoiseFlags = cms.bool(False),
        setPulseShapeFlags = cms.bool(False),
        setSaturationFlags = cms.bool(False),
        setTimingTrustFlags = cms.bool(False),
        tsFromDB = cms.bool(True),
        useLeakCorrection = cms.bool(False)
    )

    process.HLTDoLocalHcalSequence = cms.Sequence(
        process.hltHcalDigis + 
        process.hltHbhereco + 
        process.hltHfprereco + 
        process.hltHfreco + 
        process.hltHoreco
    )

    return process

def customizePhase2MuonHLTFixedGridRho(process, processName = "MYHLT"):

    process.hltTowerMakerForAll = cms.EDProducer("CaloTowersCreator",
        AllowMissingInputs = cms.bool(False),
        EBGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
        EBSumThreshold = cms.double(0.2),
        EBThreshold = cms.double(0.07),
        EBWeight = cms.double(1.0),
        EBWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
        EEGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
        EESumThreshold = cms.double(0.45),
        EEThreshold = cms.double(0.3),
        EEWeight = cms.double(1.0),
        EEWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
        EcalRecHitSeveritiesToBeExcluded = cms.vstring(
            'kTime', 
            'kWeird', 
            'kBad'
        ),
        EcalSeveritiesToBeUsedInBadTowers = cms.vstring(),
        EcutTower = cms.double(-1000.0),
        HBGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
        HBThreshold = cms.double(0.3),
        HBThreshold1 = cms.double(0.1),
        HBThreshold2 = cms.double(0.2),
        HBWeight = cms.double(1.0),
        HBWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
        HEDGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
        HEDThreshold = cms.double(0.2),
        HEDThreshold1 = cms.double(0.1),
        HEDWeight = cms.double(1.0),
        HEDWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
        HESGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
        HESThreshold = cms.double(0.2),
        HESThreshold1 = cms.double(0.1),
        HESWeight = cms.double(1.0),
        HESWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
        HF1Grid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
        HF1Threshold = cms.double(0.5),
        HF1Weight = cms.double(1.0),
        HF1Weights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
        HF2Grid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
        HF2Threshold = cms.double(0.85),
        HF2Weight = cms.double(1.0),
        HF2Weights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
        HOGrid = cms.vdouble(-1.0, 1.0, 10.0, 100.0, 1000.0),
        HOThreshold0 = cms.double(1.1),
        HOThresholdMinus1 = cms.double(3.5),
        HOThresholdMinus2 = cms.double(3.5),
        HOThresholdPlus1 = cms.double(3.5),
        HOThresholdPlus2 = cms.double(3.5),
        HOWeight = cms.double(1.0),
        HOWeights = cms.vdouble(1.0, 1.0, 1.0, 1.0, 1.0),
        HcalAcceptSeverityLevel = cms.uint32(9),
        HcalAcceptSeverityLevelForRejectedHit = cms.uint32(9999),
        HcalPhase = cms.int32(1),
        HcalThreshold = cms.double(-1000.0),
        MomConstrMethod = cms.int32(1),
        MomEBDepth = cms.double(0.3),
        MomEEDepth = cms.double(0.0),
        MomHBDepth = cms.double(0.2),
        MomHEDepth = cms.double(0.4),
        UseEcalRecoveredHits = cms.bool(False),
        UseEtEBTreshold = cms.bool(False),
        UseEtEETreshold = cms.bool(False),
        UseHO = cms.bool(False),
        UseHcalRecoveredHits = cms.bool(True),
        UseRejectedHitsOnly = cms.bool(False),
        UseRejectedRecoveredEcalHits = cms.bool(False),
        UseRejectedRecoveredHcalHits = cms.bool(True),
        UseSymEBTreshold = cms.bool(True),
        UseSymEETreshold = cms.bool(True),
        ecalInputs = cms.VInputTag("hltEcalRecHit:EcalRecHitsEB", "hltEcalRecHit:EcalRecHitsEE"),
        hbheInput = cms.InputTag("hltHbhereco"),
        hfInput = cms.InputTag("hltHfreco"),
        hoInput = cms.InputTag("hltHoreco"),
        missingHcalRescaleFactorForEcal = cms.double(0)
    )

    process.hltFixedGridRhoFastjetAllCaloForMuons = cms.EDProducer("FixedGridRhoProducerFastjet",
        gridSpacing = cms.double(0.55),
        maxRapidity = cms.double(2.5),
        pfCandidatesTag = cms.InputTag("hltTowerMakerForAll")
    )

    process.HLTFastJetForMuons = cms.Sequence(
        process.hltTowerMakerForAll + 
        process.hltFixedGridRhoFastjetAllCaloForMuons
    )

    return process

def customizePhase2MuonHLTEcalIsolation(process, processName = "MYHLT"):

    process.hltParticleFlowRecHitECALUnseeded = cms.EDProducer("PFRecHitProducer",
        navigator = cms.PSet(
            barrel = cms.PSet(

            ),
            endcap = cms.PSet(

            ),
            name = cms.string('PFRecHitECALNavigator')
        ),
        producers = cms.VPSet(
            cms.PSet(
                name = cms.string('PFEBRecHitCreator'),
                qualityTests = cms.VPSet(
                    cms.PSet(
                        applySelectionsToAllCrystals = cms.bool(True),
                        name = cms.string('PFRecHitQTestDBThreshold')
                    ), 
                    cms.PSet(
                        cleaningThreshold = cms.double(2.0),
                        name = cms.string('PFRecHitQTestECAL'),
                        skipTTRecoveredHits = cms.bool(True),
                        timingCleaning = cms.bool(True),
                        topologicalCleaning = cms.bool(True)
                    )
                ),
                srFlags = cms.InputTag(""),
                src = cms.InputTag("hltEcalRecHit","EcalRecHitsEB")
            ), 
            cms.PSet(
                name = cms.string('PFEERecHitCreator'),
                qualityTests = cms.VPSet(
                    cms.PSet(
                        applySelectionsToAllCrystals = cms.bool(True),
                        name = cms.string('PFRecHitQTestDBThreshold')
                    ), 
                    cms.PSet(
                        cleaningThreshold = cms.double(2.0),
                        name = cms.string('PFRecHitQTestECAL'),
                        skipTTRecoveredHits = cms.bool(True),
                        timingCleaning = cms.bool(True),
                        topologicalCleaning = cms.bool(True)
                    )
                ),
                srFlags = cms.InputTag(""),
                src = cms.InputTag("hltEcalRecHit","EcalRecHitsEE")
            )
        )
    )

    process.hltParticleFlowRecHitPSUnseeded = cms.EDProducer("PFRecHitProducer",
        navigator = cms.PSet(
            name = cms.string('PFRecHitPreshowerNavigator')
        ),
        producers = cms.VPSet(cms.PSet(
            name = cms.string('PFPSRecHitCreator'),
            qualityTests = cms.VPSet(
                cms.PSet(
                    name = cms.string('PFRecHitQTestThreshold'),
                    threshold = cms.double(0.0)
                ), 
                cms.PSet(
                    cleaningThreshold = cms.double(0.0),
                    name = cms.string('PFRecHitQTestES'),
                    topologicalCleaning = cms.bool(True)
                )
            ),
            src = cms.InputTag("hltEcalPreshowerRecHit","EcalRecHitsES")
        ))
    )

    process.hltParticleFlowClusterPSUnseeded = cms.EDProducer("PFClusterProducer",
        energyCorrector = cms.PSet(

        ),
        initialClusteringStep = cms.PSet(
            algoName = cms.string('Basic2DGenericTopoClusterizer'),
            thresholdsByDetector = cms.VPSet(
                cms.PSet(
                    detector = cms.string('PS1'),
                    gatheringThreshold = cms.double(6e-05),
                    gatheringThresholdPt = cms.double(0.0)
                ), 
                cms.PSet(
                    detector = cms.string('PS2'),
                    gatheringThreshold = cms.double(6e-05),
                    gatheringThresholdPt = cms.double(0.0)
                )
            ),
            useCornerCells = cms.bool(False)
        ),
        pfClusterBuilder = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowClusterizer'),
            excludeOtherSeeds = cms.bool(True),
            maxIterations = cms.uint32(50),
            minFracTot = cms.double(1e-20),
            minFractionToKeep = cms.double(1e-07),
            positionCalc = cms.PSet(
                algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
                logWeightDenominator = cms.double(6e-05),
                minAllowedNormalization = cms.double(1e-09),
                minFractionInCalc = cms.double(1e-09),
                posCalcNCrystals = cms.int32(-1)
            ),
            recHitEnergyNorms = cms.VPSet(
                cms.PSet(
                    detector = cms.string('PS1'),
                    recHitEnergyNorm = cms.double(6e-05)
                ), 
                cms.PSet(
                    detector = cms.string('PS2'),
                    recHitEnergyNorm = cms.double(6e-05)
                )
            ),
            showerSigma = cms.double(0.3),
            stoppingTolerance = cms.double(1e-08)
        ),
        positionReCalc = cms.PSet(

        ),
        recHitCleaners = cms.VPSet(),
        recHitsSource = cms.InputTag("hltParticleFlowRecHitPSUnseeded"),
        seedCleaners = cms.VPSet(),
        seedFinder = cms.PSet(
            algoName = cms.string('LocalMaximumSeedFinder'),
            nNeighbours = cms.int32(4),
            thresholdsByDetector = cms.VPSet(
                cms.PSet(
                    detector = cms.string('PS1'),
                    seedingThreshold = cms.double(0.00012),
                    seedingThresholdPt = cms.double(0.0)
                ), 
                cms.PSet(
                    detector = cms.string('PS2'),
                    seedingThreshold = cms.double(0.00012),
                    seedingThresholdPt = cms.double(0.0)
                )
            )
        )
    )

    process.hltParticleFlowClusterECALUncorrectedUnseeded = cms.EDProducer("PFClusterProducer",
        energyCorrector = cms.PSet(

        ),
        initialClusteringStep = cms.PSet(
            algoName = cms.string('Basic2DGenericTopoClusterizer'),
            thresholdsByDetector = cms.VPSet(
                cms.PSet(
                    detector = cms.string('ECAL_BARREL'),
                    gatheringThreshold = cms.double(0.08),
                    gatheringThresholdPt = cms.double(0.0)
                ), 
                cms.PSet(
                    detector = cms.string('ECAL_ENDCAP'),
                    gatheringThreshold = cms.double(0.3),
                    gatheringThresholdPt = cms.double(0.0)
                )
            ),
            useCornerCells = cms.bool(True)
        ),
        pfClusterBuilder = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowClusterizer'),
            allCellsPositionCalc = cms.PSet(
                algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
                logWeightDenominator = cms.double(0.08),
                minAllowedNormalization = cms.double(1e-09),
                minFractionInCalc = cms.double(1e-09),
                posCalcNCrystals = cms.int32(-1),
                timeResolutionCalcBarrel = cms.PSet(
                    constantTerm = cms.double(0.428192),
                    constantTermLowE = cms.double(0.0),
                    corrTermLowE = cms.double(0.0510871),
                    noiseTerm = cms.double(1.10889),
                    noiseTermLowE = cms.double(1.31883),
                    threshHighE = cms.double(5.0),
                    threshLowE = cms.double(0.5)
                ),
                timeResolutionCalcEndcap = cms.PSet(
                    constantTerm = cms.double(0.0),
                    constantTermLowE = cms.double(0.0),
                    corrTermLowE = cms.double(0.0),
                    noiseTerm = cms.double(5.72489999999),
                    noiseTermLowE = cms.double(6.92683000001),
                    threshHighE = cms.double(10.0),
                    threshLowE = cms.double(1.0)
                )
            ),
            excludeOtherSeeds = cms.bool(True),
            maxIterations = cms.uint32(50),
            minFracTot = cms.double(1e-20),
            minFractionToKeep = cms.double(1e-07),
            positionCalc = cms.PSet(
                algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
                logWeightDenominator = cms.double(0.08),
                minAllowedNormalization = cms.double(1e-09),
                minFractionInCalc = cms.double(1e-09),
                posCalcNCrystals = cms.int32(9),
                timeResolutionCalcBarrel = cms.PSet(
                    constantTerm = cms.double(0.428192),
                    constantTermLowE = cms.double(0.0),
                    corrTermLowE = cms.double(0.0510871),
                    noiseTerm = cms.double(1.10889),
                    noiseTermLowE = cms.double(1.31883),
                    threshHighE = cms.double(5.0),
                    threshLowE = cms.double(0.5)
                ),
                timeResolutionCalcEndcap = cms.PSet(
                    constantTerm = cms.double(0.0),
                    constantTermLowE = cms.double(0.0),
                    corrTermLowE = cms.double(0.0),
                    noiseTerm = cms.double(5.72489999999),
                    noiseTermLowE = cms.double(6.92683000001),
                    threshHighE = cms.double(10.0),
                    threshLowE = cms.double(1.0)
                )
            ),
            positionCalcForConvergence = cms.PSet(
                T0_EB = cms.double(7.4),
                T0_EE = cms.double(3.1),
                T0_ES = cms.double(1.2),
                W0 = cms.double(4.2),
                X0 = cms.double(0.89),
                algoName = cms.string('ECAL2DPositionCalcWithDepthCorr'),
                minAllowedNormalization = cms.double(0.0),
                minFractionInCalc = cms.double(0.0)
            ),
            recHitEnergyNorms = cms.VPSet(
                cms.PSet(
                    detector = cms.string('ECAL_BARREL'),
                    recHitEnergyNorm = cms.double(0.08)
                ), 
                cms.PSet(
                    detector = cms.string('ECAL_ENDCAP'),
                    recHitEnergyNorm = cms.double(0.3)
                )
            ),
            showerSigma = cms.double(1.5),
            stoppingTolerance = cms.double(1e-08)
        ),
        positionReCalc = cms.PSet(
            T0_EB = cms.double(7.4),
            T0_EE = cms.double(3.1),
            T0_ES = cms.double(1.2),
            W0 = cms.double(4.2),
            X0 = cms.double(0.89),
            algoName = cms.string('ECAL2DPositionCalcWithDepthCorr'),
            minAllowedNormalization = cms.double(0.0),
            minFractionInCalc = cms.double(0.0)
        ),
        recHitCleaners = cms.VPSet(),
        recHitsSource = cms.InputTag("hltParticleFlowRecHitECALUnseeded"),
        seedCleaners = cms.VPSet(cms.PSet(
            RecHitFlagsToBeExcluded = cms.vstring(),
            algoName = cms.string('FlagsCleanerECAL')
        )),
        seedFinder = cms.PSet(
            algoName = cms.string('LocalMaximumSeedFinder'),
            nNeighbours = cms.int32(8),
            thresholdsByDetector = cms.VPSet(
                cms.PSet(
                    detector = cms.string('ECAL_ENDCAP'),
                    seedingThreshold = cms.double(0.6),
                    seedingThresholdPt = cms.double(0.15)
                ), 
                cms.PSet(
                    detector = cms.string('ECAL_BARREL'),
                    seedingThreshold = cms.double(0.23),
                    seedingThresholdPt = cms.double(0.0)
                )
            )
        )
    )

    process.hltParticleFlowClusterECALUnseeded = cms.EDProducer("CorrectedECALPFClusterProducer",
        energyCorrector = cms.PSet(
            applyCrackCorrections = cms.bool(False),
            applyMVACorrections = cms.bool(True),
            autoDetectBunchSpacing = cms.bool(True),
            bunchSpacing = cms.int32(25),
            ebSrFlagLabel = cms.InputTag("hltEcalDigis"),
            eeSrFlagLabel = cms.InputTag("hltEcalDigis"),
            maxPtForMVAEvaluation = cms.double(300.0),
            recHitsEBLabel = cms.InputTag("hltEcalRecHit","EcalRecHitsEB"),
            recHitsEELabel = cms.InputTag("hltEcalRecHit","EcalRecHitsEE"),
            setEnergyUncertainty = cms.bool(False),
            srfAwareCorrection = cms.bool(True)
        ),
        inputECAL = cms.InputTag("hltParticleFlowClusterECALUncorrectedUnseeded"),
        inputPS = cms.InputTag("hltParticleFlowClusterPSUnseeded"),
        mightGet = cms.optional.untracked.vstring,
        minimumPSEnergy = cms.double(0)
    )

    process.HLTPFClusteringForMuonsUnseeded = cms.Sequence(
        process.hltParticleFlowRecHitECALUnseeded + 
        process.hltParticleFlowRecHitPSUnseeded + 
        process.hltParticleFlowClusterPSUnseeded + 
        process.hltParticleFlowClusterECALUncorrectedUnseeded + 
        process.hltParticleFlowClusterECALUnseeded
        # hltParticleFlowSuperClusterECALUnseeded
    )

    process.hltPhase2L3MuonsEcalIsoBase = cms.EDProducer("MuonHLTEcalPFClusterIsolationProducer",
        absEtaLowEdges = cms.vdouble(0.0, 1.479),
        doRhoCorrection = cms.bool(False),
        drMax = cms.double( 1.0 ),  ### 0.3
        drVetoBarrel = cms.double( 0.05 ),
        drVetoEndcap = cms.double( 0.05 ),
        effectiveAreas = cms.vdouble( 0.35, 0.193 ),
        energyBarrel = cms.double(0.0),
        energyEndcap = cms.double(0.0),
        etaStripBarrel = cms.double(0.0),
        etaStripEndcap = cms.double(0.0),
        pfClusterProducer = cms.InputTag("hltParticleFlowClusterECALUnseeded"),
        recoCandidateProducer = cms.InputTag("hltPhase2L3MuonCandidates"),  # HERE
        rhoMax = cms.double(99999999.0),
        rhoProducer = cms.InputTag("hltFixedGridRhoFastjetAllCaloForMuons"),
        rhoScale = cms.double(1.0)
    )

    return process

def customizePhase2MuonHLTHcalIsolation(process, processName = "MYHLT"):

    process.hltParticleFlowRecHitHBHEForMuons = cms.EDProducer("PFRecHitProducer",
        navigator = cms.PSet(
            name = cms.string('PFRecHitHCALNavigator'),
            sigmaCut = cms.double(4.0),
            timeResolutionCalc = cms.PSet(
                constantTerm = cms.double(1.92),
                constantTermLowE = cms.double(6.0),
                corrTermLowE = cms.double(0.0),
                noiseTerm = cms.double(8.64),
                noiseTermLowE = cms.double(0.0),
                threshHighE = cms.double(8.0),
                threshLowE = cms.double(2.0)
            )
        ),
        producers = cms.VPSet(cms.PSet(
            name = cms.string('PFHBHERecHitCreator'),
            qualityTests = cms.VPSet(
                cms.PSet(
                    cuts = cms.VPSet(
                        cms.PSet(
                            depth = cms.vint32(1, 2, 3, 4),
                            detectorEnum = cms.int32(1),
                            threshold = cms.vdouble(0.8, 0.8, 0.8, 0.8)
                        ), 
                        cms.PSet(
                            depth = cms.vint32(
                                1, 2, 3, 4, 5, 
                                6, 7
                            ),
                            detectorEnum = cms.int32(2),
                            threshold = cms.vdouble(
                                0.8, 0.8, 0.8, 0.8, 0.8, 
                                0.8, 0.8
                            )
                        )
                    ),
                    name = cms.string('PFRecHitQTestThreshold'),
                    threshold = cms.double(0.8)
                ), 
                cms.PSet(
                    cleaningThresholds = cms.vdouble(0.0),
                    flags = cms.vstring('Standard'),
                    maxSeverities = cms.vint32(11),
                    name = cms.string('PFRecHitQTestHCALChannel')
                )
            ),
            src = cms.InputTag("hltHbhereco")
        ))
    )

    process.hltParticleFlowClusterHBHEForMuons = cms.EDProducer("PFClusterProducer",
        energyCorrector = cms.PSet(

        ),
        initialClusteringStep = cms.PSet(
            algoName = cms.string('Basic2DGenericTopoClusterizer'),
            thresholdsByDetector = cms.VPSet(
                cms.PSet(
                    depths = cms.vint32(1, 2, 3, 4),
                    detector = cms.string('HCAL_BARREL1'),
                    gatheringThreshold = cms.vdouble(0.8, 0.8, 0.8, 0.8),
                    gatheringThresholdPt = cms.vdouble(0.0, 0.0, 0.0, 0.0)
                ), 
                cms.PSet(
                    depths = cms.vint32(
                        1, 2, 3, 4, 5, 
                        6, 7
                    ),
                    detector = cms.string('HCAL_ENDCAP'),
                    gatheringThreshold = cms.vdouble(
                        0.8, 0.8, 0.8, 0.8, 0.8, 
                        0.8, 0.8
                    ),
                    gatheringThresholdPt = cms.vdouble(
                        0.0, 0.0, 0.0, 0.0, 0.0, 
                        0.0, 0.0
                    )
                )
            ),
            useCornerCells = cms.bool(True)
        ),
        pfClusterBuilder = cms.PSet(
            algoName = cms.string('Basic2DGenericPFlowClusterizer'),
            allCellsPositionCalc = cms.PSet(
                algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
                logWeightDenominator = cms.double(0.8),
                minAllowedNormalization = cms.double(1e-09),
                minFractionInCalc = cms.double(1e-09),
                posCalcNCrystals = cms.int32(-1)
            ),
            clusterTimeResFromSeed = cms.bool(False),
            excludeOtherSeeds = cms.bool(True),
            maxIterations = cms.uint32(50),
            maxNSigmaTime = cms.double(10.0),
            minChi2Prob = cms.double(0.0),
            minFracTot = cms.double(1e-20),
            minFractionToKeep = cms.double(1e-07),
            positionCalc = cms.PSet(
                algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
                logWeightDenominator = cms.double(0.8),
                minAllowedNormalization = cms.double(1e-09),
                minFractionInCalc = cms.double(1e-09),
                posCalcNCrystals = cms.int32(5)
            ),
            recHitEnergyNorms = cms.VPSet(
                cms.PSet(
                    depths = cms.vint32(1, 2, 3, 4),
                    detector = cms.string('HCAL_BARREL1'),
                    recHitEnergyNorm = cms.vdouble(0.8, 0.8, 0.8, 0.8)
                ), 
                cms.PSet(
                    depths = cms.vint32(
                        1, 2, 3, 4, 5, 
                        6, 7
                    ),
                    detector = cms.string('HCAL_ENDCAP'),
                    recHitEnergyNorm = cms.vdouble(
                        0.8, 0.8, 0.8, 0.8, 0.8, 
                        0.8, 0.8
                    )
                )
            ),
            showerSigma = cms.double(10.0),
            stoppingTolerance = cms.double(1e-08),
            timeResolutionCalcBarrel = cms.PSet(
                constantTerm = cms.double(2.82),
                constantTermLowE = cms.double(4.24),
                corrTermLowE = cms.double(0.0),
                noiseTerm = cms.double(21.86),
                noiseTermLowE = cms.double(8.0),
                threshHighE = cms.double(15.0),
                threshLowE = cms.double(6.0)
            ),
            timeResolutionCalcEndcap = cms.PSet(
                constantTerm = cms.double(2.82),
                constantTermLowE = cms.double(4.24),
                corrTermLowE = cms.double(0.0),
                noiseTerm = cms.double(21.86),
                noiseTermLowE = cms.double(8.0),
                threshHighE = cms.double(15.0),
                threshLowE = cms.double(6.0)
            ),
            timeSigmaEB = cms.double(10.0),
            timeSigmaEE = cms.double(10.0)
        ),
        positionReCalc = cms.PSet(

        ),
        recHitCleaners = cms.VPSet(),
        recHitsSource = cms.InputTag("hltParticleFlowRecHitHBHEForMuons"),
        seedFinder = cms.PSet(
            algoName = cms.string('LocalMaximumSeedFinder'),
            nNeighbours = cms.int32(4),
            thresholdsByDetector = cms.VPSet(
                cms.PSet(
                    depths = cms.vint32(1, 2, 3, 4),
                    detector = cms.string('HCAL_BARREL1'),
                    seedingThreshold = cms.vdouble(1.0, 1.0, 1.0, 1.0),
                    seedingThresholdPt = cms.vdouble(0.0, 0.0, 0.0, 0.0)
                ), 
                cms.PSet(
                    depths = cms.vint32(
                        1, 2, 3, 4, 5, 
                        6, 7
                    ),
                    detector = cms.string('HCAL_ENDCAP'),
                    seedingThreshold = cms.vdouble(
                        1.1, 1.1, 1.1, 1.1, 1.1, 
                        1.1, 1.1
                    ),
                    seedingThresholdPt = cms.vdouble(
                        0.0, 0.0, 0.0, 0.0, 0.0, 
                        0.0, 0.0
                    )
                )
            )
        )
    )

    process.hltParticleFlowClusterHCALForMuons = cms.EDProducer("PFMultiDepthClusterProducer",
        clustersSource = cms.InputTag("hltParticleFlowClusterHBHEForMuons"),
        energyCorrector = cms.PSet(

        ),
        pfClusterBuilder = cms.PSet(
            algoName = cms.string('PFMultiDepthClusterizer'),
            allCellsPositionCalc = cms.PSet(
                algoName = cms.string('Basic2DGenericPFlowPositionCalc'),
                logWeightDenominator = cms.double(0.8),
                minAllowedNormalization = cms.double(1e-09),
                minFractionInCalc = cms.double(1e-09),
                posCalcNCrystals = cms.int32(-1)
            ),
            minFractionToKeep = cms.double(1e-07),
            nSigmaEta = cms.double(2.0),
            nSigmaPhi = cms.double(2.0)
        ),
        positionReCalc = cms.PSet(

        )
    )

    process.HLTPFHcalClusteringForMuons = cms.Sequence(
        # hltRegionalTowerForMuons + 
        process.hltParticleFlowRecHitHBHEForMuons + 
        process.hltParticleFlowClusterHBHEForMuons + 
        process.hltParticleFlowClusterHCALForMuons
    )

    process.hltPhase2L3MuonsHcalIsoBase = cms.EDProducer("MuonHLTHcalPFClusterIsolationProducer",
        absEtaLowEdges = cms.vdouble(0.0, 1.479),
        doRhoCorrection = cms.bool(False),
        drMax = cms.double( 1.0 ),  ### 0.3
        drVetoBarrel = cms.double( 0.1 ),
        drVetoEndcap = cms.double( 0.1 ),
        effectiveAreas = cms.vdouble( 0.227, 0.372 ),
        energyBarrel = cms.double(0.0),
        energyEndcap = cms.double(0.0),
        etaStripBarrel = cms.double(0.0),
        etaStripEndcap = cms.double(0.0),
        pfClusterProducerHCAL = cms.InputTag("hltParticleFlowClusterHCALForMuons"),
        pfClusterProducerHFEM = cms.InputTag(""),
        pfClusterProducerHFHAD = cms.InputTag(""),
        recoCandidateProducer = cms.InputTag("hltPhase2L3MuonCandidates"),  # HERE
        rhoMax = cms.double(99999999.0),
        rhoProducer = cms.InputTag("hltFixedGridRhoFastjetAllCaloForMuons"),
        rhoScale = cms.double(1.0),
        useEt = cms.bool(True),
        useHF = cms.bool(False)
    )

    return process

def customizePhase2MuonHLTHgcalPFIsolation(process, processName = "MYHLT"):

    print "not using anymore"
    import sys
    sys.exit(1)

    # process.hgcalDigis = cms.EDProducer("HGCalRawToDigiFake",
    #     bhDigis = cms.InputTag("simHGCalUnsuppressedDigis","HEback"),
    #     eeDigis = cms.InputTag("simHGCalUnsuppressedDigis","EE"),
    #     fhDigis = cms.InputTag("simHGCalUnsuppressedDigis","HEfront"),
    #     mightGet = cms.optional.untracked.vstring
    # )

    # process.HGCalUncalibRecHit = cms.EDProducer("HGCalUncalibRecHitProducer",
    #     HGCEEConfig = cms.PSet(
    #         adcNbits = cms.uint32(10),
    #         adcSaturation = cms.double(100),
    #         fCPerMIP = cms.vdouble(2.06, 3.43, 5.15),
    #         isSiFE = cms.bool(True),
    #         tdcNbits = cms.uint32(12),
    #         tdcOnset = cms.double(60),
    #         tdcSaturation = cms.double(10000),
    #         toaLSB_ns = cms.double(0.0244)
    #     ),
    #     HGCEEdigiCollection = cms.InputTag("hgcalDigis","EE"),
    #     HGCEEhitCollection = cms.string('HGCEEUncalibRecHits'),
    #     HGCHEBConfig = cms.PSet(
    #         adcNbits = cms.uint32(10),
    #         adcSaturation = cms.double(68.75),
    #         fCPerMIP = cms.vdouble(1.0, 1.0, 1.0),
    #         isSiFE = cms.bool(True),
    #         tdcNbits = cms.uint32(12),
    #         tdcOnset = cms.double(55),
    #         tdcSaturation = cms.double(1000),
    #         toaLSB_ns = cms.double(0.0244)
    #     ),
    #     HGCHEBdigiCollection = cms.InputTag("hgcalDigis","HEback"),
    #     HGCHEBhitCollection = cms.string('HGCHEBUncalibRecHits'),
    #     HGCHEFConfig = cms.PSet(
    #         adcNbits = cms.uint32(10),
    #         adcSaturation = cms.double(100),
    #         fCPerMIP = cms.vdouble(2.06, 3.43, 5.15),
    #         isSiFE = cms.bool(True),
    #         tdcNbits = cms.uint32(12),
    #         tdcOnset = cms.double(60),
    #         tdcSaturation = cms.double(10000),
    #         toaLSB_ns = cms.double(0.0244)
    #     ),
    #     HGCHEFdigiCollection = cms.InputTag("hgcalDigis","HEfront"),
    #     HGCHEFhitCollection = cms.string('HGCHEFUncalibRecHits'),
    #     HGCHFNoseConfig = cms.PSet(
    #         adcNbits = cms.uint32(10),
    #         adcSaturation = cms.double(100),
    #         fCPerMIP = cms.vdouble(1.25, 2.57, 3.88),
    #         isSiFE = cms.bool(False),
    #         tdcNbits = cms.uint32(12),
    #         tdcOnset = cms.double(60),
    #         tdcSaturation = cms.double(10000),
    #         toaLSB_ns = cms.double(0.0244)
    #     ),
    #     HGCHFNosedigiCollection = cms.InputTag("hfnoseDigis","HFNose"),
    #     HGCHFNosehitCollection = cms.string('HGCHFNoseUncalibRecHits'),
    #     algo = cms.string('HGCalUncalibRecHitWorkerWeights')
    # )

    # process.HGCalRecHit = cms.EDProducer("HGCalRecHitProducer",
    #     HGCEE_cce = cms.PSet(
    #         refToPSet_ = cms.string('HGCAL_chargeCollectionEfficiencies')
    #     ),
    #     HGCEE_fCPerMIP = cms.vdouble(2.06, 3.43, 5.15),
    #     HGCEE_isSiFE = cms.bool(True),
    #     HGCEE_keV2DIGI = cms.double(0.044259),
    #     HGCEE_noise_fC = cms.PSet(
    #         refToPSet_ = cms.string('HGCAL_noise_fC')
    #     ),
    #     HGCEErechitCollection = cms.string('HGCEERecHits'),
    #     HGCEEuncalibRecHitCollection = cms.InputTag("HGCalUncalibRecHit","HGCEEUncalibRecHits"),
    #     HGCHEB_isSiFE = cms.bool(True),
    #     HGCHEB_keV2DIGI = cms.double(0.00148148148148),
    #     HGCHEB_noise_MIP = cms.PSet(
    #         refToPSet_ = cms.string('HGCAL_noise_heback')
    #     ),
    #     HGCHEBrechitCollection = cms.string('HGCHEBRecHits'),
    #     HGCHEBuncalibRecHitCollection = cms.InputTag("HGCalUncalibRecHit","HGCHEBUncalibRecHits"),
    #     HGCHEF_cce = cms.PSet(
    #         refToPSet_ = cms.string('HGCAL_chargeCollectionEfficiencies')
    #     ),
    #     HGCHEF_fCPerMIP = cms.vdouble(2.06, 3.43, 5.15),
    #     HGCHEF_isSiFE = cms.bool(True),
    #     HGCHEF_keV2DIGI = cms.double(0.044259),
    #     HGCHEF_noise_fC = cms.PSet(
    #         refToPSet_ = cms.string('HGCAL_noise_fC')
    #     ),
    #     HGCHEFrechitCollection = cms.string('HGCHEFRecHits'),
    #     HGCHEFuncalibRecHitCollection = cms.InputTag("HGCalUncalibRecHit","HGCHEFUncalibRecHits"),
    #     HGCHFNose_cce = cms.PSet(
    #         refToPSet_ = cms.string('HGCAL_chargeCollectionEfficiencies')
    #     ),
    #     HGCHFNose_fCPerMIP = cms.vdouble(1.25, 2.57, 3.88),
    #     HGCHFNose_isSiFE = cms.bool(False),
    #     HGCHFNose_keV2DIGI = cms.double(0.044259),
    #     HGCHFNose_noise_fC = cms.PSet(
    #         refToPSet_ = cms.string('HGCAL_noise_fC')
    #     ),
    #     HGCHFNoserechitCollection = cms.string('HGCHFNoseRecHits'),
    #     HGCHFNoseuncalibRecHitCollection = cms.InputTag("HGCalUncalibRecHit","HGCHFNoseUncalibRecHits"),
    #     algo = cms.string('HGCalRecHitWorkerSimple'),
    #     constSiPar = cms.double(0.02),
    #     deltasi_index_regemfac = cms.int32(3),
    #     layerNoseWeights = cms.vdouble(
    #         0.0, 39.500245, 39.756638, 39.756638, 39.756638, 
    #         39.756638, 66.020266, 92.283895, 92.283895
    #     ),
    #     layerWeights = cms.vdouble(
    #         0.0, 8.894541, 10.937907, 10.937907, 10.937907, 
    #         10.937907, 10.937907, 10.937907, 10.937907, 10.937907, 
    #         10.932882, 10.932882, 10.937907, 10.937907, 10.938169, 
    #         10.938169, 10.938169, 10.938169, 10.938169, 10.938169, 
    #         10.938169, 10.938169, 10.938169, 10.938169, 10.938169, 
    #         10.938169, 10.938169, 10.938169, 32.332097, 51.574301, 
    #         51.444192, 51.444192, 51.444192, 51.444192, 51.444192, 
    #         51.444192, 51.444192, 51.444192, 51.444192, 51.444192, 
    #         69.513118, 87.582044, 87.582044, 87.582044, 87.582044, 
    #         87.582044, 87.214571, 86.888309, 86.92952, 86.92952, 
    #         86.92952
    #     ),
    #     maxValSiPar = cms.double(10000.0),
    #     minValSiPar = cms.double(10.0),
    #     noiseSiPar = cms.double(5.5),
    #     rangeMask = cms.uint32(4294442496),
    #     rangeMatch = cms.uint32(1161838592),
    #     sciThicknessCorrection = cms.double(0.9),
    #     thicknessCorrection = cms.vdouble(
    #         0.77, 0.77, 0.77, 0.84, 0.84, 
    #         0.84
    #     ),
    #     thicknessNoseCorrection = cms.vdouble(1.132, 1.092, 1.084)
    # )

    # process.hgcalLayerClusters = cms.EDProducer("HGCalLayerClusterProducer",
    #     HFNoseInput = cms.InputTag("HGCalRecHit","HGCHFNoseRecHits"),
    #     HGCBHInput = cms.InputTag("HGCalRecHit","HGCHEBRecHits"),
    #     HGCEEInput = cms.InputTag("HGCalRecHit","HGCEERecHits"),
    #     HGCFHInput = cms.InputTag("HGCalRecHit","HGCHEFRecHits"),
    #     detector = cms.string('all'),
    #     doSharing = cms.bool(False),
    #     mightGet = cms.optional.untracked.vstring,
    #     nHitsTime = cms.uint32(3),
    #     plugin = cms.PSet(
    #         dEdXweights = cms.vdouble(
    #             0.0, 8.894541, 10.937907, 10.937907, 10.937907, 
    #             10.937907, 10.937907, 10.937907, 10.937907, 10.937907, 
    #             10.932882, 10.932882, 10.937907, 10.937907, 10.938169, 
    #             10.938169, 10.938169, 10.938169, 10.938169, 10.938169, 
    #             10.938169, 10.938169, 10.938169, 10.938169, 10.938169, 
    #             10.938169, 10.938169, 10.938169, 32.332097, 51.574301, 
    #             51.444192, 51.444192, 51.444192, 51.444192, 51.444192, 
    #             51.444192, 51.444192, 51.444192, 51.444192, 51.444192, 
    #             69.513118, 87.582044, 87.582044, 87.582044, 87.582044, 
    #             87.582044, 87.214571, 86.888309, 86.92952, 86.92952, 
    #             86.92952
    #         ),
    #         deltac = cms.vdouble(1.3, 1.3, 5, 0.0315),
    #         dependSensor = cms.bool(True),
    #         ecut = cms.double(3),
    #         fcPerEle = cms.double(0.00016020506),
    #         fcPerMip = cms.vdouble(2.06, 3.43, 5.15),
    #         kappa = cms.double(9),
    #         noiseMip = cms.PSet(
    #             refToPSet_ = cms.string('HGCAL_noise_heback')
    #         ),
    #         noises = cms.PSet(
    #             refToPSet_ = cms.string('HGCAL_noises')
    #         ),
    #         thicknessCorrection = cms.vdouble(
    #             0.77, 0.77, 0.77, 0.84, 0.84, 
    #             0.84
    #         ),
    #         thresholdW0 = cms.vdouble(2.9, 2.9, 2.9),
    #         type = cms.string('CLUE'),
    #         use2x2 = cms.bool(True),
    #         verbosity = cms.untracked.uint32(3)
    #     ),
    #     timeClname = cms.string('timeLayerCluster'),
    #     timeOffset = cms.double(5)
    # )

    # process.filteredLayerClustersTrk = cms.EDProducer("FilteredLayerClustersProducer",
    #     HGCLayerClusters = cms.InputTag("hgcalLayerClusters"),
    #     LayerClustersInputMask = cms.InputTag("hgcalLayerClusters","InitialLayerClustersMask"),
    #     algo_number = cms.int32(8),
    #     clusterFilter = cms.string('ClusterFilterByAlgo'),
    #     iteration_label = cms.string('Trk'),
    #     max_cluster_size = cms.int32(9999),
    #     mightGet = cms.optional.untracked.vstring,
    #     min_cluster_size = cms.int32(0)
    # )

    # process.ticlLayerTileProducer = cms.EDProducer("TICLLayerTileProducer",
    #     layer_clusters = cms.InputTag("hgcalLayerClusters"),
    #     mightGet = cms.optional.untracked.vstring
    # )

    # process.ticlSeedingTrk = cms.EDProducer("TICLSeedingRegionProducer",
    #     algoId = cms.int32(2),
    #     algo_verbosity = cms.int32(0),
    #     cutTk = cms.string('1.48 < abs(eta) < 3.0 && pt > 1. && quality("highPurity") && hitPattern().numberOfLostHits("MISSING_OUTER_HITS") < 10'),
    #     mightGet = cms.optional.untracked.vstring,
    #     propagator = cms.string('PropagatorWithMaterial'),
    #     tracks = cms.InputTag("")
    # )

    # process.ticlTrackstersTrk = cms.EDProducer("TrackstersProducer",
    #     algo_verbosity = cms.int32(2),
    #     eid_graph_path = cms.string('RecoHGCal/TICL/data/tf_models/energy_id_v0.pb'),
    #     eid_input_name = cms.string('input'),
    #     eid_min_cluster_energy = cms.double(1),
    #     eid_n_clusters = cms.int32(10),
    #     eid_n_layers = cms.int32(50),
    #     eid_output_name_energy = cms.string('output/regressed_energy'),
    #     eid_output_name_id = cms.string('output/id_probabilities'),
    #     etaLimitIncreaseWindow = cms.double(2.1),
    #     filter_on_categories = cms.vint32(2, 4),
    #     filtered_mask = cms.InputTag("filteredLayerClustersTrk","Trk"),
    #     itername = cms.string('TRK'),
    #     layer_clusters = cms.InputTag("hgcalLayerClusters"),
    #     layer_clusters_tiles = cms.InputTag("ticlLayerTileProducer"),
    #     max_delta_time = cms.double(-1.0),
    #     max_out_in_hops = cms.int32(10),
    #     mightGet = cms.optional.untracked.vstring,
    #     min_clusters_per_ntuplet = cms.int32(10),
    #     min_cos_pointing = cms.double(0.798),
    #     min_cos_theta = cms.double(0.866),
    #     missing_layers = cms.int32(3),
    #     oneTracksterPerTrackSeed = cms.bool(True),
    #     original_mask = cms.InputTag("hgcalLayerClusters","InitialLayerClustersMask"),
    #     out_in_dfs = cms.bool(True),
    #     pid_threshold = cms.double(0),
    #     promoteEmptyRegionToTrackster = cms.bool(True),
    #     seeding_regions = cms.InputTag("ticlSeedingTrk"),
    #     time_layerclusters = cms.InputTag("hgcalLayerClusters","timeLayerCluster")
    # )

    # process.filteredLayerClustersEM = cms.EDProducer("FilteredLayerClustersProducer",
    #     HGCLayerClusters = cms.InputTag("hgcalLayerClusters"),
    #     LayerClustersInputMask = cms.InputTag("ticlTrackstersTrk"),
    #     algo_number = cms.int32(8),
    #     clusterFilter = cms.string('ClusterFilterByAlgoAndSize'),
    #     iteration_label = cms.string('EM'),
    #     max_cluster_size = cms.int32(9999),
    #     mightGet = cms.optional.untracked.vstring,
    #     min_cluster_size = cms.int32(2)
    # )

    # process.ticlSeedingGlobal = cms.EDProducer("TICLSeedingRegionProducer",
    #     algoId = cms.int32(2),
    #     algo_verbosity = cms.int32(0),
    #     cutTk = cms.string('1.48 < abs(eta) < 3.0 && pt > 1. && quality("highPurity") && hitPattern().numberOfLostHits("MISSING_OUTER_HITS") < 10'),
    #     mightGet = cms.optional.untracked.vstring,
    #     propagator = cms.string('PropagatorWithMaterial'),
    #     tracks = cms.InputTag("")
    # )

    # process.ticlTrackstersEM = cms.EDProducer("TrackstersProducer",
    #     algo_verbosity = cms.int32(0),
    #     eid_graph_path = cms.string('RecoHGCal/TICL/data/tf_models/energy_id_v0.pb'),
    #     eid_input_name = cms.string('input'),
    #     eid_min_cluster_energy = cms.double(1),
    #     eid_n_clusters = cms.int32(10),
    #     eid_n_layers = cms.int32(50),
    #     eid_output_name_energy = cms.string('output/regressed_energy'),
    #     eid_output_name_id = cms.string('output/id_probabilities'),
    #     etaLimitIncreaseWindow = cms.double(2.1),
    #     filter_on_categories = cms.vint32(0, 1),
    #     filtered_mask = cms.InputTag("filteredLayerClustersEM","EM"),
    #     itername = cms.string('EM'),
    #     layer_clusters = cms.InputTag("hgcalLayerClusters"),
    #     layer_clusters_tiles = cms.InputTag("ticlLayerTileProducer"),
    #     max_delta_time = cms.double(3),
    #     max_out_in_hops = cms.int32(4),
    #     mightGet = cms.optional.untracked.vstring,
    #     min_clusters_per_ntuplet = cms.int32(10),
    #     min_cos_pointing = cms.double(0.9),
    #     min_cos_theta = cms.double(0.978),
    #     missing_layers = cms.int32(1),
    #     oneTracksterPerTrackSeed = cms.bool(False),
    #     original_mask = cms.InputTag("ticlTrackstersTrk"),
    #     out_in_dfs = cms.bool(True),
    #     pid_threshold = cms.double(0.8),
    #     promoteEmptyRegionToTrackster = cms.bool(False),
    #     seeding_regions = cms.InputTag("ticlSeedingGlobal"),
    #     time_layerclusters = cms.InputTag("hgcalLayerClusters","timeLayerCluster")
    # )

    # process.ticlMultiClustersFromTrackstersEM = cms.EDProducer("MultiClustersFromTrackstersProducer",
    #     LayerClusters = cms.InputTag("hgcalLayerClusters"),
    #     Tracksters = cms.InputTag("ticlTrackstersEM"),
    #     mightGet = cms.optional.untracked.vstring,
    #     verbosity = cms.untracked.uint32(3)
    # )

    # process.particleFlowRecHitHGC = cms.EDProducer("PFRecHitProducer",
    #     navigator = cms.PSet(
    #         hgcee = cms.PSet(
    #             name = cms.string('PFRecHitHGCEENavigator'),
    #             topologySource = cms.string('HGCalEESensitive')
    #         ),
    #         hgcheb = cms.PSet(
    #             name = cms.string('PFRecHitHGCHENavigator'),
    #             topologySource = cms.string('HGCalHEScintillatorSensitive')
    #         ),
    #         hgchef = cms.PSet(
    #             name = cms.string('PFRecHitHGCHENavigator'),
    #             topologySource = cms.string('HGCalHESiliconSensitive')
    #         ),
    #         name = cms.string('PFRecHitHGCNavigator')
    #     ),
    #     producers = cms.VPSet(
    #         cms.PSet(
    #             geometryInstance = cms.string('HGCalEESensitive'),
    #             name = cms.string('PFHGCalEERecHitCreator'),
    #             qualityTests = cms.VPSet(cms.PSet(
    #                 name = cms.string('PFRecHitQTestHGCalThresholdSNR'),
    #                 thresholdSNR = cms.double(5.0)
    #             )),
    #             src = cms.InputTag("HGCalRecHit","HGCEERecHits")
    #         ), 
    #         cms.PSet(
    #             geometryInstance = cms.string('HGCalHESiliconSensitive'),
    #             name = cms.string('PFHGCalHSiRecHitCreator'),
    #             qualityTests = cms.VPSet(cms.PSet(
    #                 name = cms.string('PFRecHitQTestHGCalThresholdSNR'),
    #                 thresholdSNR = cms.double(5.0)
    #             )),
    #             src = cms.InputTag("HGCalRecHit","HGCHEFRecHits")
    #         ), 
    #         cms.PSet(
    #             geometryInstance = cms.string(''),
    #             name = cms.string('PFHGCalHScRecHitCreator'),
    #             qualityTests = cms.VPSet(cms.PSet(
    #                 name = cms.string('PFRecHitQTestHGCalThresholdSNR'),
    #                 thresholdSNR = cms.double(5.0)
    #             )),
    #             src = cms.InputTag("HGCalRecHit","HGCHEBRecHits")
    #         )
    #     )
    # )

    # process.particleFlowClusterHGCalFromTICL = cms.EDProducer("PFClusterProducer",
    #     energyCorrector = cms.PSet(

    #     ),
    #     initialClusteringStep = cms.PSet(
    #         algoName = cms.string('PFClusterFromHGCalMultiCluster'),
    #         clusterSrc = cms.InputTag("ticlMultiClustersFromTrackstersEM"),
    #         thresholdsByDetector = cms.VPSet()
    #     ),
    #     pfClusterBuilder = cms.PSet(

    #     ),
    #     positionReCalc = cms.PSet(
    #         algoName = cms.string('Cluster3DPCACalculator'),
    #         minFractionInCalc = cms.double(1e-09),
    #         updateTiming = cms.bool(False)
    #     ),
    #     recHitCleaners = cms.VPSet(),
    #     recHitsSource = cms.InputTag("particleFlowRecHitHGC"),
    #     seedCleaners = cms.VPSet(),
    #     seedFinder = cms.PSet(
    #         algoName = cms.string('PassThruSeedFinder'),
    #         nNeighbours = cms.int32(8),
    #         thresholdsByDetector = cms.VPSet()
    #     )
    # )

    # process.HLTHgcalTiclPFClusteringForMuons = cms.Sequence(
    #     process.hgcalDigis + 
    #     process.HGCalUncalibRecHit + 
    #     process.HGCalRecHit + 
    #     process.hgcalLayerClusters + 
    #     process.filteredLayerClustersTrk + 
    #     process.ticlLayerTileProducer + 
    #     process.ticlSeedingTrk + 
    #     process.ticlTrackstersTrk + 
    #     process.filteredLayerClustersEM + 
    #     process.ticlSeedingGlobal + 
    #     process.ticlTrackstersEM + 
    #     process.ticlMultiClustersFromTrackstersEM + 
    #     process.particleFlowRecHitHGC + 
    #     process.particleFlowClusterHGCalFromTICL
    #     # offlineBeamSpot + 
    #     # particleFlowSuperClusterHGCalFromTICL
    # )

    # process.hltPhase2L3MuonsHgcalIsoBase = cms.EDProducer("MuonHLTEcalPFClusterIsolationProducer",
    #     absEtaLowEdges = cms.vdouble(0.0, 1.479),
    #     doRhoCorrection = cms.bool(False),
    #     drMax = cms.double( 1.0 ),  ### 0.3
    #     drVetoBarrel = cms.double(0.05),
    #     drVetoEndcap = cms.double(0.05),
    #     effectiveAreas = cms.vdouble(0.29, 0.21),
    #     energyBarrel = cms.double(0.0),
    #     energyEndcap = cms.double(0.0),
    #     etaStripBarrel = cms.double(0.0),
    #     etaStripEndcap = cms.double(0.0),
    #     pfClusterProducer = cms.InputTag("particleFlowClusterHGCalFromTICL"),
    #     recoCandidateProducer = cms.InputTag("hltPhase2L3MuonCandidates"),  # HERE
    #     rhoMax = cms.double(99999999.0),
    #     rhoProducer = cms.InputTag("hltFixedGridRhoFastjetAllCaloForMuons"),
    #     rhoScale = cms.double(1.0)
    # )

    return process

def customizePhase2MuonHLTHgcalLCIsolation(process, processName = "MYHLT"):

    process.HLTHgcalTiclLayerClusteringForMuons = cms.Sequence(
        process.hgcalDigis + 
        process.HGCalUncalibRecHit + 
        process.HGCalRecHit + 
        process.hgcalLayerClusters
    )

    process.hltPhase2L3MuonsHgcalLCIsoBase = cms.EDProducer("MuonHLTHGCalLayerClusterIsolationProducer",
        doRhoCorrection = cms.bool(False),
        drMax = cms.double( 1.0 ),
        drVetoEM = cms.double(0.0),
        drVetoHad = cms.double(0.0),
        effectiveAreas = cms.vdouble(0.0, 0.0),
        minEnergyEM = cms.double(0.0),
        minEnergyHad = cms.double(0.0),
        minEtEM = cms.double(0.0),
        minEtHad = cms.double(0.0),
        layerClusterProducer = cms.InputTag("hgcalLayerClusters"),
        recoCandidateProducer = cms.InputTag("hltPhase2L3MuonCandidates"),  # HERE
        rhoMax = cms.double(99999999.0),
        rhoProducer = cms.InputTag("hltFixedGridRhoFastjetAllCaloForMuons"),
        rhoScale = cms.double(1.0),
        useEt = cms.bool(False)
    )

    return process

def customizePhase2MuonHLTTrkIsolationRegional(process, processName = "MYHLT"):

    process.hltPhase2L3MuonVertex = cms.EDProducer( "VertexFromTrackProducer",
        verbose = cms.untracked.bool( False ),
        useTriggerFilterElectrons = cms.bool( False ),
        beamSpotLabel = cms.InputTag( "offlineBeamSpot" ),
        isRecoCandidate = cms.bool( True ),
        trackLabel = cms.InputTag( "hltPhase2L3MuonCandidates" ),
        useTriggerFilterMuons = cms.bool( False ),
        useBeamSpot = cms.bool( True ),
        vertexLabel = cms.InputTag( "notUsed" ),
        triggerFilterElectronsSrc = cms.InputTag( "notUsed" ),
        triggerFilterMuonsSrc = cms.InputTag( "notUsed" ),
        useVertex = cms.bool( False )
    )

    process.hltPixelLayerQuadruplets = cms.EDProducer("SeedingLayersEDProducer",
        BPix = cms.PSet(
            HitProducer = cms.string('siPixelRecHits'),
            TTRHBuilder = cms.string('WithTrackAngle'),
            #hitErrorRPhi = cms.double(0.0027),
            #hitErrorRZ = cms.double(0.006),
            #useErrorsFromParam = cms.bool(True)
        ),
        FPix = cms.PSet(
            HitProducer = cms.string('siPixelRecHits'),
            TTRHBuilder = cms.string('WithTrackAngle'),
            #hitErrorRPhi = cms.double(0.0051),
            #hitErrorRZ = cms.double(0.0036),
            #useErrorsFromParam = cms.bool(True)
        ),
        MTEC = cms.PSet(

        ),
        MTIB = cms.PSet(

        ),
        MTID = cms.PSet(

        ),
        MTOB = cms.PSet(

        ),
        TEC = cms.PSet(

        ),
        TIB = cms.PSet(

        ),
        TID = cms.PSet(

        ),
        TOB = cms.PSet(

        ),
        layerList = cms.vstring(
            'BPix1+BPix2+BPix3+BPix4', 
            'BPix1+BPix2+BPix3+FPix1_pos', 
            'BPix1+BPix2+BPix3+FPix1_neg', 
            'BPix1+BPix2+FPix1_pos+FPix2_pos', 
            'BPix1+BPix2+FPix1_neg+FPix2_neg', 
            'BPix1+FPix1_pos+FPix2_pos+FPix3_pos', 
            'BPix1+FPix1_neg+FPix2_neg+FPix3_neg', 
            'FPix1_pos+FPix2_pos+FPix3_pos+FPix4_pos', 
            'FPix1_neg+FPix2_neg+FPix3_neg+FPix4_neg', 
            'FPix2_pos+FPix3_pos+FPix4_pos+FPix5_pos', 
            'FPix2_neg+FPix3_neg+FPix4_neg+FPix5_neg', 
            'FPix3_pos+FPix4_pos+FPix5_pos+FPix6_pos', 
            'FPix3_neg+FPix4_neg+FPix5_neg+FPix6_neg', 
            'FPix4_pos+FPix5_pos+FPix6_pos+FPix7_pos', 
            'FPix4_neg+FPix5_neg+FPix6_neg+FPix7_neg', 
            'FPix5_pos+FPix6_pos+FPix7_pos+FPix8_pos', 
            'FPix5_neg+FPix6_neg+FPix7_neg+FPix8_neg'
        )
    )

    process.hltPixelTracksPhase2L3MuonFilter = cms.EDProducer( "PixelTrackFilterByKinematicsProducer",
        chi2 = cms.double( 1000.0 ),
        nSigmaTipMaxTolerance = cms.double( 0.0 ),
        ptMin = cms.double( 0.9 ), ##before it was 0.1
        nSigmaInvPtTolerance = cms.double( 0.0 ),
        tipMax = cms.double( 1.0 )
    )

    process.hltPixelTracksPhase2L3MuonFitter = cms.EDProducer( "PixelFitterByHelixProjectionsProducer",
        scaleErrorsForBPix1 = cms.bool( False ),
        scaleFactor = cms.double( 0.65 )
    )

    process.hltPixelTracksTrackingRegionsPhase2L3Muon = cms.EDProducer( "GlobalTrackingRegionWithVerticesEDProducer",
        RegionPSet = cms.PSet( 
          useFixedError = cms.bool( True ),
          nSigmaZ = cms.double( 4.0 ),
          VertexCollection = cms.InputTag( "hltPhase2L3MuonVertex" ),
          beamSpot = cms.InputTag( "offlineBeamSpot" ),
          useFoundVertices = cms.bool( True ),
          fixedError = cms.double( 0.5 ),
          sigmaZVertex = cms.double( 4.0 ),
          useFakeVertices = cms.bool( True ),
          ptMin = cms.double( 0.9 ),
          originRadius = cms.double( 0.2 ),
          precise = cms.bool( True ),
          useMultipleScattering = cms.bool( False )
        )
    )

    process.hltPixelTracksHitDoubletsPhase2L3Muon = cms.EDProducer( "HitPairEDProducer",
        trackingRegions = cms.InputTag( "hltPixelTracksTrackingRegionsPhase2L3Muon" ),
        layerPairs = cms.vuint32( 0, 1, 2 ),
        clusterCheck = cms.InputTag( "" ),
        produceSeedingHitSets = cms.bool( False ),
        produceIntermediateHitDoublets = cms.bool( True ),
        trackingRegionsSeedingLayers = cms.InputTag( "" ),
        maxElement = cms.uint32( 0 ),
        seedingLayers = cms.InputTag( "hltPixelLayerQuadruplets" )
    )

    process.hltPixelTracksHitQuadrupletsPhase2L3Muon = cms.EDProducer( "CAHitQuadrupletEDProducer",
        CAThetaCut = cms.double( 0.002 ),
        SeedComparitorPSet = cms.PSet( 
          clusterShapeHitFilter = cms.string( "ClusterShapeHitFilter" ),
          ComponentName = cms.string( "LowPtClusterShapeSeedComparitor" ),
          clusterShapeCacheSrc = cms.InputTag( "siPixelClusterShapeCache" )
        ),
        extraHitRPhitolerance = cms.double( 0.032 ),
        doublets = cms.InputTag( "hltPixelTracksHitDoubletsPhase2L3Muon" ),
        fitFastCircle = cms.bool( True ),
        CAHardPtCut = cms.double( 0.0 ),
        maxChi2 = cms.PSet( 
          value2 = cms.double( 50.0 ),
          value1 = cms.double( 200.0 ),
          pt1 = cms.double( 0.7 ),
          enabled = cms.bool( True ),
          pt2 = cms.double( 2.0 )
        ),
        CAPhiCut = cms.double( 0.2 ),
        useBendingCorrection = cms.bool( True ),
        fitFastCircleChi2Cut = cms.bool( True )
    )

    process.hltPixelTracksPhase2L3Muon = cms.EDProducer( "PixelTrackProducer",
        Filter = cms.InputTag( "hltPixelTracksPhase2L3MuonFilter" ),
        Cleaner = cms.string( "hltPixelTracksCleanerBySharedHits" ),
        passLabel = cms.string( "" ),
        Fitter = cms.InputTag( "hltPixelTracksPhase2L3MuonFitter" ),
        SeedingHitSets = cms.InputTag( "hltPixelTracksHitQuadrupletsPhase2L3Muon" )
    )

    process.hltPixelVerticesPhase2L3Muon = cms.EDProducer( "PixelVertexProducer",
        WtAverage = cms.bool( True ),
        Method2 = cms.bool( True ),
        beamSpot = cms.InputTag( "offlineBeamSpot" ),
        PVcomparer = cms.PSet(  refToPSet_ = cms.string( "hltPhase2PSetPvClusterComparerForIT" ) ),
        Verbosity = cms.int32( 0 ),
        UseError = cms.bool( True ),
        TrackCollection = cms.InputTag( "hltPixelTracksPhase2L3Muon" ),
        PtMin = cms.double( 1.0 ),
        NTrkMin = cms.int32( 2 ),
        ZOffset = cms.double( 5.0 ),
        Finder = cms.string( "DivisiveVertexFinder" ),
        ZSeparation = cms.double( 0.05 )
    )

    process.HLTPixelTrackingPhase2L3Muon = cms.Sequence(
        process.hltPhase2L3MuonVertex +
        process.HLTDoLocalPixelSequence + 
        process.hltPixelLayerQuadruplets + 
        process.hltPixelTracksPhase2L3MuonFilter +
        process.hltPixelTracksPhase2L3MuonFitter + 
        process.hltPixelTracksTrackingRegionsPhase2L3Muon + 
        process.hltPixelTracksHitDoubletsPhase2L3Muon + 
        process.hltPixelTracksHitQuadrupletsPhase2L3Muon + 
        process.hltPixelTracksPhase2L3Muon + 
        process.hltPixelVerticesPhase2L3Muon
    )

    process.hltPixelTracksForSeedsPhase2L3MuonFilter = cms.EDProducer( "PixelTrackFilterByKinematicsProducer",
        chi2 = cms.double( 1000.0 ),
        nSigmaTipMaxTolerance = cms.double( 0.0 ),
        ptMin = cms.double( 0.9 ), ##before it was 0.1
        nSigmaInvPtTolerance = cms.double( 0.0 ),
        tipMax = cms.double( 1.0 )
    )

    process.hltPixelTracksForSeedsPhase2L3MuonFitter = cms.EDProducer( "PixelFitterByHelixProjectionsProducer",
        scaleErrorsForBPix1 = cms.bool( False ),
        scaleFactor = cms.double( 0.65 )
    )

    process.hltPixelTracksTrackingRegionsForSeedsPhase2L3Muon = cms.EDProducer( "CandidateSeededTrackingRegionsEDProducer",
        RegionPSet = cms.PSet( 
          vertexCollection = cms.InputTag( "hltPixelVerticesPhase2L3Muon" ),
          zErrorVetex = cms.double( 0.2 ),
          beamSpot = cms.InputTag( "offlineBeamSpot" ),
          zErrorBeamSpot = cms.double( 24.2 ),
          maxNVertices = cms.int32( 1 ),
          maxNRegions = cms.int32( 10 ),
          nSigmaZVertex = cms.double( 3.0 ),
          nSigmaZBeamSpot = cms.double( 4.0 ),
          ptMin = cms.double( 0.9 ),
          mode = cms.string( "VerticesFixed" ),
          input = cms.InputTag( "hltPhase2L3MuonCandidates" ),
          searchOpt = cms.bool( False ),
          whereToUseMeasurementTracker = cms.string( "Never" ),
          originRadius = cms.double( 0.1 ),
          measurementTrackerName = cms.InputTag( "" ),
          precise = cms.bool( True ),
          deltaEta = cms.double( 0.3 ),
          deltaPhi = cms.double( 0.3 )
        )
    )

    process.hltPixelTracksHitDoubletsForSeedsPhase2L3Muon = cms.EDProducer( "HitPairEDProducer",
        trackingRegions = cms.InputTag( "hltPixelTracksTrackingRegionsForSeedsPhase2L3Muon" ),
        layerPairs = cms.vuint32( 0, 1, 2 ),
        clusterCheck = cms.InputTag( "" ),
        produceSeedingHitSets = cms.bool( False ),
        produceIntermediateHitDoublets = cms.bool( True ),
        trackingRegionsSeedingLayers = cms.InputTag( "" ),
        maxElement = cms.uint32( 0 ),
        seedingLayers = cms.InputTag( "hltPixelLayerQuadruplets" )
    )

    process.hltPixelTracksHitQuadrupletsForSeedsPhase2L3Muon = cms.EDProducer( "CAHitQuadrupletEDProducer",
        CAThetaCut = cms.double( 0.002 ),
        SeedComparitorPSet = cms.PSet( 
          clusterShapeHitFilter = cms.string( "ClusterShapeHitFilter" ),
          ComponentName = cms.string( "LowPtClusterShapeSeedComparitor" ),
          clusterShapeCacheSrc = cms.InputTag( "siPixelClusterShapeCache" )
        ),
        extraHitRPhitolerance = cms.double( 0.032 ),
        doublets = cms.InputTag( "hltPixelTracksHitDoubletsForSeedsPhase2L3Muon" ),
        fitFastCircle = cms.bool( True ),
        CAHardPtCut = cms.double( 0.0 ),
        maxChi2 = cms.PSet( 
          value2 = cms.double( 50.0 ),
          value1 = cms.double( 200.0 ),
          pt1 = cms.double( 0.7 ),
          enabled = cms.bool( True ),
          pt2 = cms.double( 2.0 )
        ),
        CAPhiCut = cms.double( 0.2 ),
        useBendingCorrection = cms.bool( True ),
        fitFastCircleChi2Cut = cms.bool( True )
    )

    process.hltPixelTracksForSeedsPhase2L3Muon = cms.EDProducer( "PixelTrackProducer",
        Filter = cms.InputTag( "hltPixelTracksForSeedsPhase2L3MuonFilter" ),
        Cleaner = cms.string( "hltPixelTracksCleanerBySharedHits" ),
        passLabel = cms.string( "" ),
        Fitter = cms.InputTag( "hltPixelTracksForSeedsPhase2L3MuonFitter" ),
        SeedingHitSets = cms.InputTag( "hltPixelTracksHitQuadrupletsForSeedsPhase2L3Muon" )
    )

    process.hltIter0Phase2L3MuonPixelSeedsFromPixelTracks = cms.EDProducer( "SeedGeneratorFromProtoTracksEDProducer",
        useEventsWithNoVertex = cms.bool( True ),
        originHalfLength = cms.double( 0.2 ),
        useProtoTrackKinematics = cms.bool( False ),
        usePV = cms.bool( False ),
        SeedCreatorPSet = cms.PSet(  refToPSet_ = cms.string( "hltPhase2SeedFromProtoTracks" ) ),
        InputVertexCollection = cms.InputTag( "hltPixelVerticesPhase2L3Muon" ),
        TTRHBuilder = cms.string( "WithTrackAngle" ),
        InputCollection = cms.InputTag( "hltPixelTracksForSeedsPhase2L3Muon" ),
        originRadius = cms.double( 0.1 )
    )

    process.hltIter0Phase2L3MuonCkfTrackCandidates = cms.EDProducer( "CkfTrackCandidateMaker",
        src = cms.InputTag( "hltIter0Phase2L3MuonPixelSeedsFromPixelTracks" ),
        maxSeedsBeforeCleaning = cms.uint32( 1000 ),
        SimpleMagneticField = cms.string( "ParabolicMf" ),
        TransientInitialStateEstimatorParameters = cms.PSet( 
          propagatorAlongTISE = cms.string( "PropagatorWithMaterialParabolicMf" ),
          numberMeasurementsForFit = cms.int32( 4 ),
          propagatorOppositeTISE = cms.string( "PropagatorWithMaterialParabolicMfOpposite" )
        ),
        TrajectoryCleaner = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
        MeasurementTrackerEvent = cms.InputTag( "MeasurementTrackerEvent" ),
        cleanTrajectoryAfterInOut = cms.bool( False ),
        useHitsSplitting = cms.bool( False ),
        RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
        doSeedingRegionRebuilding = cms.bool( False ),
        maxNSeeds = cms.uint32( 100000 ),
        TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTIter0GroupedCkfTrajectoryBuilderIT" ) ),
        NavigationSchool = cms.string( "SimpleNavigationSchool" ),
        TrajectoryBuilder = cms.string( "" )
    )

    process.hltIter0Phase2L3MuonCtfWithMaterialTracks = cms.EDProducer( "TrackProducer",
        src = cms.InputTag( "hltIter0Phase2L3MuonCkfTrackCandidates" ),
        SimpleMagneticField = cms.string(''),
        # SimpleMagneticField = cms.string( "ParabolicMf" ),
        clusterRemovalInfo = cms.InputTag( "" ),
        beamSpot = cms.InputTag( "offlineBeamSpot" ),
        MeasurementTrackerEvent = cms.InputTag( "MeasurementTrackerEvent" ),
        Fitter = cms.string('FlexibleKFFittingSmoother'),
        # Fitter = cms.string('KFFittingSmootherForInOut'),
        useHitsSplitting = cms.bool( False ),
        MeasurementTracker = cms.string( "" ),
        AlgorithmName = cms.string( "hltIterX" ),
        alias = cms.untracked.string( "ctfWithMaterialTracks" ),
        NavigationSchool = cms.string( "" ),
        TrajectoryInEvent = cms.bool( False ),
        TTRHBuilder = cms.string( "WithTrackAngle" ),
        GeometricInnerState = cms.bool( True ),
        useSimpleMF = cms.bool(False),
        # useSimpleMF = cms.bool( True ),
        Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" )
    )

    process.hltIter0Phase2L3MuonTrackCutClassifier = cms.EDProducer( "TrackCutClassifier",
        src = cms.InputTag( "hltIter0Phase2L3MuonCtfWithMaterialTracks" ),
        beamspot = cms.InputTag( "offlineBeamSpot" ),
        vertices = cms.InputTag( "hltPixelVerticesPhase2L3Muon" ),
        qualityCuts = cms.vdouble( -0.7, 0.1, 0.7 ),
        mva = cms.PSet( 
          minPixelHits = cms.vint32( 0, 3, 4 ),
          maxDzWrtBS = cms.vdouble( 3.40282346639E38, 24.0, 15.0 ),
          dr_par = cms.PSet( 
            d0err = cms.vdouble( 0.003, 0.003, 0.003 ),
            dr_par2 = cms.vdouble( 0.3, 0.3, 0.3 ),
            dr_par1 = cms.vdouble( 0.4, 0.4, 0.4 ),
            dr_exp = cms.vint32( 4, 4, 4 ),
            d0err_par = cms.vdouble( 0.001, 0.001, 0.001 )
          ),
          maxLostLayers = cms.vint32( 1, 1, 1 ),
          min3DLayers = cms.vint32( 0, 3, 4 ),
          dz_par = cms.PSet( 
            dz_par1 = cms.vdouble( 0.4, 0.4, 0.4 ),
            dz_par2 = cms.vdouble( 0.35, 0.35, 0.35 ),
            dz_exp = cms.vint32( 4, 4, 4 )
          ),
          minNVtxTrk = cms.int32( 3 ),
          maxDz = cms.vdouble( 0.5, 0.2, 3.40282346639E38 ),
          minNdof = cms.vdouble( 1.0E-5, 1.0E-5, 1.0E-5 ),
          maxChi2 = cms.vdouble( 9999.0, 25.0, 16.0 ),
          maxChi2n = cms.vdouble( 1.2, 1.0, 0.7 ),
          maxDr = cms.vdouble( 0.5, 0.03, 3.40282346639E38 ),
          minLayers = cms.vint32( 3, 3, 4 )
        ),
        ignoreVertices = cms.bool( False )
    )

    process.hltIter0Phase2L3MuonTrackSelectionHighPurity = cms.EDProducer( "TrackCollectionFilterCloner",
        minQuality = cms.string( "highPurity" ),
        copyExtras = cms.untracked.bool( True ),
        copyTrajectories = cms.untracked.bool( False ),
        originalSource = cms.InputTag( "hltIter0Phase2L3MuonCtfWithMaterialTracks" ),
        originalQualVals = cms.InputTag( 'hltIter0Phase2L3MuonTrackCutClassifier','QualityMasks' ),
        originalMVAVals = cms.InputTag( 'hltIter0Phase2L3MuonTrackCutClassifier','MVAValues' )
    )

    process.HLTIterativeTrackingPhase2L3MuonIteration0 = cms.Sequence(
        process.hltPixelTracksForSeedsPhase2L3MuonFilter +
        process.hltPixelTracksForSeedsPhase2L3MuonFitter +
        process.hltPixelTracksTrackingRegionsForSeedsPhase2L3Muon +
        process.hltPixelTracksHitDoubletsForSeedsPhase2L3Muon +
        process.hltPixelTracksHitQuadrupletsForSeedsPhase2L3Muon +
        process.hltPixelTracksForSeedsPhase2L3Muon +
        process.hltIter0Phase2L3MuonPixelSeedsFromPixelTracks +
        process.hltIter0Phase2L3MuonCkfTrackCandidates +
        process.hltIter0Phase2L3MuonCtfWithMaterialTracks +
        process.hltIter0Phase2L3MuonTrackCutClassifier +
        process.hltIter0Phase2L3MuonTrackSelectionHighPurity
    )

    process.hltIter2Phase2L3MuonClustersRefRemoval = cms.EDProducer("TrackClusterRemoverPhase2",
        trackClassifier = cms.InputTag( '','QualityMasks' ),
        minNumberOfLayersWithMeasBeforeFiltering = cms.int32( 0 ),
        maxChi2 = cms.double( 16.0 ),
        trajectories = cms.InputTag( "hltIter0Phase2L3MuonTrackSelectionHighPurity" ),
        oldClusterRemovalInfo = cms.InputTag( "" ),
        # trajectories = cms.InputTag( "hltIter1Phase2L3MuonTrackSelectionHighPurity" ),
        # oldClusterRemovalInfo = cms.InputTag( "hltIter1Phase2L3MuonClustersRefRemoval" ),
        phase2OTClusters = cms.InputTag( "siPhase2Clusters" ),
        overrideTrkQuals = cms.InputTag( "" ),
        phase2pixelClusters = cms.InputTag( "siPixelClusters" ),
        TrackQuality = cms.string( "highPurity" )
    )

    process.hltIter2Phase2L3MuonMaskedMeasurementTrackerEvent = cms.EDProducer( "MaskedMeasurementTrackerEventProducer",
        phase2clustersToSkip = cms.InputTag( "hltIter2Phase2L3MuonClustersRefRemoval" ),
        OnDemand = cms.bool( False ),
        src = cms.InputTag( "MeasurementTrackerEvent" )
    )

    process.hltIter2Phase2L3MuonPixelLayerTriplets = cms.EDProducer( "SeedingLayersEDProducer",
        layerList = cms.vstring( 'BPix1+BPix2+BPix3',
          'BPix2+BPix3+BPix4',
          'BPix1+BPix3+BPix4',
          'BPix1+BPix2+BPix4',
          'BPix2+BPix3+FPix1_pos',
          'BPix2+BPix3+FPix1_neg',
          'BPix1+BPix2+FPix1_pos',
          'BPix1+BPix2+FPix1_neg',
          'BPix2+FPix1_pos+FPix2_pos',
          'BPix2+FPix1_neg+FPix2_neg',
          'BPix1+FPix1_pos+FPix2_pos',
          'BPix1+FPix1_neg+FPix2_neg',
          'FPix1_pos+FPix2_pos+FPix3_pos',
          'FPix1_neg+FPix2_neg+FPix3_neg' ),
        MTOB = cms.PSet(  ),
        TEC = cms.PSet(  ),
        MTID = cms.PSet(  ),
        FPix = cms.PSet( 
          #hitErrorRPhi = cms.double( 0.0051 ),
          TTRHBuilder = cms.string( "WithTrackAngle" ),
          skipClusters = cms.InputTag( "hltIter2Phase2L3MuonClustersRefRemoval" ),
          #useErrorsFromParam = cms.bool( True ),
          #hitErrorRZ = cms.double( 0.0036 ),
          HitProducer = cms.string( "siPixelRecHits" )
        ),
        MTEC = cms.PSet(  ),
        MTIB = cms.PSet(  ),
        TID = cms.PSet(  ),
        TOB = cms.PSet(  ),
        BPix = cms.PSet( 
          #hitErrorRPhi = cms.double( 0.0027 ),
          TTRHBuilder = cms.string( "WithTrackAngle" ),
          skipClusters = cms.InputTag( "hltIter2Phase2L3MuonClustersRefRemoval" ),
          #useErrorsFromParam = cms.bool( True ),
          #hitErrorRZ = cms.double( 0.006 ),
          HitProducer = cms.string( "siPixelRecHits" )
        ),
        TIB = cms.PSet(  )
    )

    process.hltIter2Phase2L3MuonPixelTrackingRegions = cms.EDProducer( "CandidateSeededTrackingRegionsEDProducer",
        RegionPSet = cms.PSet( 
          vertexCollection = cms.InputTag( "hltPixelVerticesPhase2L3Muon" ),
          zErrorVetex = cms.double( 0.05 ),
          beamSpot = cms.InputTag( "offlineBeamSpot" ),
          zErrorBeamSpot = cms.double( 24.2 ),
          maxNVertices = cms.int32( 1 ),
          maxNRegions = cms.int32( 10 ),
          nSigmaZVertex = cms.double( 3.0 ),
          nSigmaZBeamSpot = cms.double( 4.0 ),
          ptMin = cms.double( 0.9 ), ## it was 0.8
          mode = cms.string( "VerticesFixed" ),
          input = cms.InputTag( "hltPhase2L3MuonCandidates" ),
          searchOpt = cms.bool( False ),
          whereToUseMeasurementTracker = cms.string( "ForSiStrips" ),
          originRadius = cms.double( 0.025 ),
          measurementTrackerName = cms.InputTag( "hltIter2Phase2L3MuonMaskedMeasurementTrackerEvent" ),
          precise = cms.bool( True ),
          deltaEta = cms.double( 0.3 ),
          deltaPhi = cms.double( 0.3 )
        )
    )

    process.hltIter2Phase2L3MuonPixelClusterCheck = cms.EDProducer( "ClusterCheckerEDProducer",
        cut = cms.string( "" ),
        silentClusterCheck = cms.untracked.bool( False ),
        MaxNumberOfCosmicClusters = cms.uint32( 50000 ),
        PixelClusterCollectionLabel = cms.InputTag( "siPixelClusters" ),
        doClusterCheck = cms.bool( False ),
        MaxNumberOfPixelClusters = cms.uint32( 10000 ),
        ClusterCollectionLabel = cms.InputTag( "MeasurementTrackerEvent" )
    )

    process.hltIter2Phase2L3MuonPixelHitDoublets = cms.EDProducer( "HitPairEDProducer",
        trackingRegions = cms.InputTag( "hltIter2Phase2L3MuonPixelTrackingRegions" ),
        layerPairs = cms.vuint32( 0, 1 ),
        clusterCheck = cms.InputTag( "hltIter2Phase2L3MuonPixelClusterCheck" ),
        produceSeedingHitSets = cms.bool( False ),
        produceIntermediateHitDoublets = cms.bool( True ),
        trackingRegionsSeedingLayers = cms.InputTag( "" ),
        maxElement = cms.uint32( 0 ),
        seedingLayers = cms.InputTag( "hltIter2Phase2L3MuonPixelLayerTriplets" )
    )

    process.hltIter2Phase2L3MuonPixelHitTriplets = cms.EDProducer( "CAHitTripletEDProducer",
        CAHardPtCut = cms.double( 0.3 ),
        SeedComparitorPSet = cms.PSet(  ComponentName = cms.string( "none" ) ),
        extraHitRPhitolerance = cms.double( 0.032 ),
        doublets = cms.InputTag( "hltIter2Phase2L3MuonPixelHitDoublets" ),
        CAThetaCut = cms.double( 0.004 ),
        maxChi2 = cms.PSet( 
          value2 = cms.double( 6.0 ),
          value1 = cms.double( 100.0 ),
          pt1 = cms.double( 0.8 ),
          enabled = cms.bool( True ),
          pt2 = cms.double( 8.0 )
        ),
        CAPhiCut = cms.double( 0.1 ),
        useBendingCorrection = cms.bool( True )
    )

    process.hltIter2Phase2L3MuonPixelSeeds = cms.EDProducer( "SeedCreatorFromRegionConsecutiveHitsEDProducer",
        SeedComparitorPSet = cms.PSet(  ComponentName = cms.string( "none" ) ),
        forceKinematicWithRegionDirection = cms.bool( False ),
        magneticField = cms.string( "ParabolicMf" ),
        SeedMomentumForBOFF = cms.double( 5.0 ),
        OriginTransverseErrorMultiplier = cms.double( 1.0 ),
        TTRHBuilder = cms.string( "WithTrackAngle" ),
        MinOneOverPtError = cms.double( 1.0 ),
        seedingHitSets = cms.InputTag( "hltIter2Phase2L3MuonPixelHitTriplets" ),
        propagator = cms.string( "PropagatorWithMaterialParabolicMf" )
    )

    process.hltIter2Phase2L3MuonCkfTrackCandidates = cms.EDProducer( "CkfTrackCandidateMaker",
        src = cms.InputTag( "hltIter2Phase2L3MuonPixelSeeds" ),
        maxSeedsBeforeCleaning = cms.uint32( 1000 ),
        SimpleMagneticField = cms.string( "ParabolicMf" ),
        TransientInitialStateEstimatorParameters = cms.PSet( 
          propagatorAlongTISE = cms.string( "PropagatorWithMaterialParabolicMf" ),
          numberMeasurementsForFit = cms.int32( 4 ),
          propagatorOppositeTISE = cms.string( "PropagatorWithMaterialParabolicMfOpposite" )
        ),
        TrajectoryCleaner = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
        MeasurementTrackerEvent = cms.InputTag( "hltIter2Phase2L3MuonMaskedMeasurementTrackerEvent" ),
        cleanTrajectoryAfterInOut = cms.bool( False ),
        useHitsSplitting = cms.bool( False ),
        RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
        doSeedingRegionRebuilding = cms.bool( False ),
        maxNSeeds = cms.uint32( 100000 ),
        TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTIter2GroupedCkfTrajectoryBuilderIT" ) ),
        NavigationSchool = cms.string( "SimpleNavigationSchool" ),
        TrajectoryBuilder = cms.string( "" )
    )

    process.hltIter2Phase2L3MuonCtfWithMaterialTracks = cms.EDProducer( "TrackProducer",
        src = cms.InputTag( "hltIter2Phase2L3MuonCkfTrackCandidates" ),
        SimpleMagneticField = cms.string(''),
        # SimpleMagneticField = cms.string( "ParabolicMf" ),
        clusterRemovalInfo = cms.InputTag( "" ),
        beamSpot = cms.InputTag( "offlineBeamSpot" ),
        MeasurementTrackerEvent = cms.InputTag( "hltIter2Phase2L3MuonMaskedMeasurementTrackerEvent" ),
        Fitter = cms.string('FlexibleKFFittingSmoother'),
        # Fitter = cms.string( "KFFittingSmootherForInOut" ),
        useHitsSplitting = cms.bool( False ),
        MeasurementTracker = cms.string( "" ),
        AlgorithmName = cms.string( "hltIterX" ),
        alias = cms.untracked.string( "ctfWithMaterialTracks" ),
        NavigationSchool = cms.string( "" ),
        TrajectoryInEvent = cms.bool( False ),
        TTRHBuilder = cms.string( "WithTrackAngle" ),
        GeometricInnerState = cms.bool( True ),
        useSimpleMF = cms.bool(False),
        # useSimpleMF = cms.bool( True ),
        Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" )
    )

    process.hltIter2Phase2L3MuonTrackCutClassifier = cms.EDProducer( "TrackCutClassifier",
        src = cms.InputTag( "hltIter2Phase2L3MuonCtfWithMaterialTracks" ),
        beamspot = cms.InputTag( "offlineBeamSpot" ),
        vertices = cms.InputTag( "hltPixelVerticesPhase2L3Muon" ),
        qualityCuts = cms.vdouble( -0.7, 0.1, 0.7 ),
        mva = cms.PSet( 
          minPixelHits = cms.vint32( 0, 0, 0 ),
          maxDzWrtBS = cms.vdouble( 3.40282346639E38, 24.0, 15.0 ),
          dr_par = cms.PSet( 
            d0err = cms.vdouble( 0.003, 0.003, 0.003 ),
            dr_par2 = cms.vdouble( 3.40282346639E38, 0.3, 0.3 ),
            dr_par1 = cms.vdouble( 3.40282346639E38, 0.4, 0.4 ),
            dr_exp = cms.vint32( 4, 4, 4 ),
            d0err_par = cms.vdouble( 0.001, 0.001, 0.001 )
          ),
          maxLostLayers = cms.vint32( 1, 1, 1 ),
          min3DLayers = cms.vint32( 0, 0, 0 ),
          dz_par = cms.PSet( 
            dz_par1 = cms.vdouble( 3.40282346639E38, 0.4, 0.4 ),
            dz_par2 = cms.vdouble( 3.40282346639E38, 0.35, 0.35 ),
            dz_exp = cms.vint32( 4, 4, 4 )
          ),
          minNVtxTrk = cms.int32( 3 ),
          maxDz = cms.vdouble( 0.5, 0.2, 3.40282346639E38 ),
          minNdof = cms.vdouble( 1.0E-5, 1.0E-5, 1.0E-5 ),
          maxChi2 = cms.vdouble( 9999.0, 25.0, 16.0 ),
          maxChi2n = cms.vdouble( 1.2, 1.0, 0.7 ),
          maxDr = cms.vdouble( 0.5, 0.03, 3.40282346639E38 ),
          minLayers = cms.vint32( 3, 3, 3 )
        ),
        ignoreVertices = cms.bool( False )
    )

    process.hltIter2Phase2L3MuonTrackSelectionHighPurity = cms.EDProducer( "TrackCollectionFilterCloner",
        minQuality = cms.string( "highPurity" ),
        copyExtras = cms.untracked.bool( True ),
        copyTrajectories = cms.untracked.bool( False ),
        originalSource = cms.InputTag( "hltIter2Phase2L3MuonCtfWithMaterialTracks" ),
        originalQualVals = cms.InputTag( 'hltIter2Phase2L3MuonTrackCutClassifier','QualityMasks' ),
        originalMVAVals = cms.InputTag( 'hltIter2Phase2L3MuonTrackCutClassifier','MVAValues' )
    )

    process.HLTIterativeTrackingPhase2L3MuonIteration2 = cms.Sequence(
        process.hltIter2Phase2L3MuonClustersRefRemoval + 
        process.hltIter2Phase2L3MuonMaskedMeasurementTrackerEvent + 
        process.hltIter2Phase2L3MuonPixelLayerTriplets + 
        process.hltIter2Phase2L3MuonPixelTrackingRegions + 
        process.hltIter2Phase2L3MuonPixelClusterCheck + 
        process.hltIter2Phase2L3MuonPixelHitDoublets + 
        process.hltIter2Phase2L3MuonPixelHitTriplets + 
        process.hltIter2Phase2L3MuonPixelSeeds + 
        process.hltIter2Phase2L3MuonCkfTrackCandidates + 
        process.hltIter2Phase2L3MuonCtfWithMaterialTracks + 
        process.hltIter2Phase2L3MuonTrackCutClassifier + 
        process.hltIter2Phase2L3MuonTrackSelectionHighPurity
    )

    process.hltIter2Phase2L3MuonMerged = cms.EDProducer( "TrackListMerger",
        ShareFrac = cms.double( 0.19 ),
        writeOnlyTrkQuals = cms.bool( False ),
        MinPT = cms.double( 0.05 ),
        allowFirstHitShare = cms.bool( True ),
        copyExtras = cms.untracked.bool( True ),
        Epsilon = cms.double( -0.001 ),
        #selectedTrackQuals = cms.VInputTag( 'hltIter1Phase2L3MuonMerged','hltIter2Phase2L3MuonTrackSelectionHighPurity' ),
        selectedTrackQuals = cms.VInputTag( 'hltIter0Phase2L3MuonTrackSelectionHighPurity','hltIter2Phase2L3MuonTrackSelectionHighPurity' ),
        indivShareFrac = cms.vdouble( 1.0, 1.0 ),
        MaxNormalizedChisq = cms.double( 1000.0 ),
        copyMVA = cms.bool( False ),
        FoundHitBonus = cms.double( 5.0 ),
        LostHitPenalty = cms.double( 20.0 ),
        setsToMerge = cms.VPSet( 
          cms.PSet(  pQual = cms.bool( False ),
            tLists = cms.vint32( 0, 1 )
          )
        ),
        MinFound = cms.int32( 3 ),
        hasSelector = cms.vint32( 0, 0 ),
        #TrackProducers = cms.VInputTag( 'hltIter1Phase2L3MuonMerged','hltIter2Phase2L3MuonTrackSelectionHighPurity' ),
        TrackProducers = cms.VInputTag( 'hltIter0Phase2L3MuonTrackSelectionHighPurity','hltIter2Phase2L3MuonTrackSelectionHighPurity' ),
        trackAlgoPriorityOrder = cms.string( "hltESPTrackAlgoPriorityOrder" ),
        newQuality = cms.string( "confirmed" )
    )

    process.HLTIterativeTrackingPhase2L3MuonIter02 = cms.Sequence(
        process.HLTIterativeTrackingPhase2L3MuonIteration0 +
        process.HLTIterativeTrackingPhase2L3MuonIteration2 +
        process.hltIter2Phase2L3MuonMerged
    )

    process.HLTTrackReconstructionForIsoPhase2L3MuonIter02 = cms.Sequence( 
        process.HLTPixelTrackingPhase2L3Muon +
        process.HLTIterativeTrackingPhase2L3MuonIter02
    )

    process.hltMuonTkRelIsolationCut0p07Map = cms.EDProducer( "L3MuonCombinedRelativeIsolationProducer",
        printDebug = cms.bool( False ),
        CutsPSet = cms.PSet( 
          applyCutsORmaxNTracks = cms.bool( False ),
          maxNTracks = cms.int32( -1 ),
          Thresholds = cms.vdouble( 0.07 ),
          EtaBounds = cms.vdouble( 2.411 ),
          ComponentName = cms.string( "SimpleCuts" ),
          ConeSizes = cms.vdouble( 0.3 )
        ),
        OutputMuIsoDeposits = cms.bool( True ),
        TrackPt_Min = cms.double( -1.0 ),
        CaloDepositsLabel = cms.InputTag( "notUsed" ),
        CaloExtractorPSet = cms.PSet( 
          DR_Veto_H = cms.double( 0.1 ),
          Vertex_Constraint_Z = cms.bool( False ),
          DR_Veto_E = cms.double( 0.07 ),
          Weight_H = cms.double( 1.0 ),
          CaloTowerCollectionLabel = cms.InputTag( "hltTowerMakerForAll" ),
          DR_Max = cms.double( 0.3 ),
          DepositLabel = cms.untracked.string( "EcalPlusHcal" ),
          Vertex_Constraint_XY = cms.bool( False ),
          Threshold_H = cms.double( 0.5 ),
          Threshold_E = cms.double( 0.2 ),
          ComponentName = cms.string( "CaloExtractor" ),
          Weight_E = cms.double( 1.0 )
        ),
        inputMuonCollection = cms.InputTag( "hltPhase2L3MuonCandidates" ),
        TrkExtractorPSet = cms.PSet( 
          Diff_z = cms.double( 0.2 ),
          inputTrackCollection = cms.InputTag( "hltIter2Phase2L3MuonMerged" ),
          Chi2Ndof_Max = cms.double( 1.0E64 ),
          BeamSpotLabel = cms.InputTag( "offlineBeamSpot" ),
          DR_Veto = cms.double( 0.01 ),
          Pt_Min = cms.double( -1.0 ),
          VetoLeadingTrack = cms.bool( True ),
          DR_Max = cms.double( 0.3 ),
          DepositLabel = cms.untracked.string( "PXLS" ),
          PtVeto_Min = cms.double( 2.0 ),
          NHits_Min = cms.uint32( 0 ),
          PropagateTracksToRadius = cms.bool( True ),
          ReferenceRadius = cms.double( 6.0 ),
          Chi2Prob_Min = cms.double( -1.0 ),
          Diff_r = cms.double( 0.1 ),
          BeamlineOption = cms.string( "BeamSpotFromEvent" ),
          ComponentName = cms.string( "PixelTrackExtractor" ),
          DR_VetoPt = cms.double( 0.025 )
        ),
        UseRhoCorrectedCaloDeposits = cms.bool( False ),
        UseCaloIso = cms.bool( False )
    )

    process.hltPhase2L3MuonsTrkIsoBase = cms.EDProducer( "L3MuonCombinedRelativeIsolationProducer",
        printDebug = cms.bool( False ),
        CutsPSet = cms.PSet( 
          applyCutsORmaxNTracks = cms.bool( False ),
          maxNTracks = cms.int32( -1 ),
          Thresholds = cms.vdouble( 1e9 ),
          EtaBounds = cms.vdouble( 2.411 ),
          ComponentName = cms.string( "SimpleCuts" ),
          ConeSizes = cms.vdouble( 1.0 )  ### 0.3
        ),
        OutputMuIsoDeposits = cms.bool( True ),
        TrackPt_Min = cms.double( -1.0 ),
        CaloDepositsLabel = cms.InputTag( "notUsed" ),
        CaloExtractorPSet = cms.PSet( 
          DR_Veto_H = cms.double( 0.1 ),
          Vertex_Constraint_Z = cms.bool( False ),
          DR_Veto_E = cms.double( 0.07 ),
          Weight_H = cms.double( 1.0 ),
          CaloTowerCollectionLabel = cms.InputTag( "hltTowerMakerForAll" ),
          DR_Max = cms.double( 0.3 ),
          DepositLabel = cms.untracked.string( "EcalPlusHcal" ),
          Vertex_Constraint_XY = cms.bool( False ),
          Threshold_H = cms.double( 0.5 ),
          Threshold_E = cms.double( 0.2 ),
          ComponentName = cms.string( "CaloExtractor" ),
          Weight_E = cms.double( 1.0 )
        ),
        inputMuonCollection = cms.InputTag( "hltPhase2L3MuonCandidates" ),
        TrkExtractorPSet = cms.PSet( 
          Diff_z = cms.double( 0.2 ),
          inputTrackCollection = cms.InputTag( "hltIter2Phase2L3MuonMerged" ),
          Chi2Ndof_Max = cms.double( 1.0E64 ),
          BeamSpotLabel = cms.InputTag( "offlineBeamSpot" ),
          DR_Veto = cms.double( 0.01 ),
          Pt_Min = cms.double( -1.0 ),
          VetoLeadingTrack = cms.bool( True ),
          DR_Max = cms.double( 1.0 ),  ### 0.3
          DepositLabel = cms.untracked.string( "PXLS" ),
          PtVeto_Min = cms.double( 2.0 ),
          NHits_Min = cms.uint32( 0 ),
          PropagateTracksToRadius = cms.bool( True ),
          ReferenceRadius = cms.double( 6.0 ),
          Chi2Prob_Min = cms.double( -1.0 ),
          Diff_r = cms.double( 0.1 ),
          BeamlineOption = cms.string( "BeamSpotFromEvent" ),
          ComponentName = cms.string( "PixelTrackExtractor" ),
          DR_VetoPt = cms.double( 0.025 )
        ),
        UseRhoCorrectedCaloDeposits = cms.bool( False ),
        UseCaloIso = cms.bool( False )
    )

    return process

def customizePhase2MuonHLTTrkIsolationFull(process, processName = "MYHLT"):

    from HLTrigger.PhaseII.Muon.Customizers.phase2_tracking import customise_hltPhase2_TRKv06_1
    process = customise_hltPhase2_TRKv06_1(process)

    process.tracking_v6_1_NoVertexReco = cms.Sequence(
        process.trackerClusterCheck
        + process.hltPhase2PixelTracksSequence
        + process.hltPhase2PixelVerticesSequence
        + process.hltPhase2InitialStepSequence
        + process.hltPhase2HighPtTripletStepSequence
        + process.hltPhase2GeneralTracks
    )

    process.hltPhase2L3MuonsTrkIsoFullBase = cms.EDProducer( "L3MuonCombinedRelativeIsolationProducer",
        printDebug = cms.bool( False ),
        CutsPSet = cms.PSet( 
          applyCutsORmaxNTracks = cms.bool( False ),
          maxNTracks = cms.int32( -1 ),
          Thresholds = cms.vdouble( 1e9 ),
          EtaBounds = cms.vdouble( 2.411 ),
          ComponentName = cms.string( "SimpleCuts" ),
          ConeSizes = cms.vdouble( 1.0 )  ### 0.3
        ),
        OutputMuIsoDeposits = cms.bool( True ),
        TrackPt_Min = cms.double( -1.0 ),
        CaloDepositsLabel = cms.InputTag( "notUsed" ),
        CaloExtractorPSet = cms.PSet( 
          DR_Veto_H = cms.double( 0.1 ),
          Vertex_Constraint_Z = cms.bool( False ),
          DR_Veto_E = cms.double( 0.07 ),
          Weight_H = cms.double( 1.0 ),
          CaloTowerCollectionLabel = cms.InputTag( "hltTowerMakerForAll" ),
          DR_Max = cms.double( 0.3 ),
          DepositLabel = cms.untracked.string( "EcalPlusHcal" ),
          Vertex_Constraint_XY = cms.bool( False ),
          Threshold_H = cms.double( 0.5 ),
          Threshold_E = cms.double( 0.2 ),
          ComponentName = cms.string( "CaloExtractor" ),
          Weight_E = cms.double( 1.0 )
        ),
        inputMuonCollection = cms.InputTag( "hltPhase2L3MuonCandidates" ),
        TrkExtractorPSet = cms.PSet( 
          Diff_z = cms.double( 0.2 ),
          inputTrackCollection = cms.InputTag( "hltPhase2GeneralTracks" ),
          Chi2Ndof_Max = cms.double( 1.0E64 ),
          BeamSpotLabel = cms.InputTag( "offlineBeamSpot" ),
          DR_Veto = cms.double( 0.01 ),
          Pt_Min = cms.double( -1.0 ),
          VetoLeadingTrack = cms.bool( True ),
          DR_Max = cms.double( 1.0 ),  ### 0.3
          DepositLabel = cms.untracked.string( "PXLS" ),
          PtVeto_Min = cms.double( 2.0 ),
          NHits_Min = cms.uint32( 0 ),
          PropagateTracksToRadius = cms.bool( True ),
          ReferenceRadius = cms.double( 6.0 ),
          Chi2Prob_Min = cms.double( -1.0 ),
          Diff_r = cms.double( 0.1 ),
          BeamlineOption = cms.string( "BeamSpotFromEvent" ),
          ComponentName = cms.string( "PixelTrackExtractor" ),
          DR_VetoPt = cms.double( 0.025 )
        ),
        UseRhoCorrectedCaloDeposits = cms.bool( False ),
        UseCaloIso = cms.bool( False )
    )

    return process

def customizePhase2MuonHLTTrkIsolationRegionalNew(process, processName = "MYHLT"):

    from HLTrigger.PhaseII.Muon.Customizers.phase2_tracking_L3Muon import customise_hltPhase2_TRKv06_1_L3Muon
    process = customise_hltPhase2_TRKv06_1_L3Muon(process)

    process.tracking_v6_1_L3Muon_NoVertexReco = cms.Sequence(
        process.trackerClusterCheck
        + process.hltPhase2L3MuonPixelTracksSequence
        + process.hltPhase2L3MuonPixelVerticesSequence
        + process.hltPhase2L3MuonInitialStepSequence
        + process.hltPhase2L3MuonHighPtTripletStepSequence
        + process.hltPhase2L3MuonGeneralTracks
    )

    process.hltPhase2L3MuonsTrkIsoRegionalNewBase = cms.EDProducer( "L3MuonCombinedRelativeIsolationProducer",
        printDebug = cms.bool( False ),
        CutsPSet = cms.PSet( 
          applyCutsORmaxNTracks = cms.bool( False ),
          maxNTracks = cms.int32( -1 ),
          Thresholds = cms.vdouble( 1e9 ),
          EtaBounds = cms.vdouble( 2.411 ),
          ComponentName = cms.string( "SimpleCuts" ),
          ConeSizes = cms.vdouble( 1.0 )  ### 0.3
        ),
        OutputMuIsoDeposits = cms.bool( True ),
        TrackPt_Min = cms.double( -1.0 ),
        CaloDepositsLabel = cms.InputTag( "notUsed" ),
        CaloExtractorPSet = cms.PSet( 
          DR_Veto_H = cms.double( 0.1 ),
          Vertex_Constraint_Z = cms.bool( False ),
          DR_Veto_E = cms.double( 0.07 ),
          Weight_H = cms.double( 1.0 ),
          CaloTowerCollectionLabel = cms.InputTag( "hltTowerMakerForAll" ),
          DR_Max = cms.double( 0.3 ),
          DepositLabel = cms.untracked.string( "EcalPlusHcal" ),
          Vertex_Constraint_XY = cms.bool( False ),
          Threshold_H = cms.double( 0.5 ),
          Threshold_E = cms.double( 0.2 ),
          ComponentName = cms.string( "CaloExtractor" ),
          Weight_E = cms.double( 1.0 )
        ),
        inputMuonCollection = cms.InputTag( "hltPhase2L3MuonCandidates" ),
        TrkExtractorPSet = cms.PSet( 
          Diff_z = cms.double( 0.2 ),
          inputTrackCollection = cms.InputTag( "hltPhase2L3MuonGeneralTracks" ),
          Chi2Ndof_Max = cms.double( 1.0E64 ),
          BeamSpotLabel = cms.InputTag( "offlineBeamSpot" ),
          DR_Veto = cms.double( 0.01 ),
          Pt_Min = cms.double( -1.0 ),
          VetoLeadingTrack = cms.bool( True ),
          DR_Max = cms.double( 1.0 ),  ### 0.3
          DepositLabel = cms.untracked.string( "PXLS" ),
          PtVeto_Min = cms.double( 2.0 ),
          NHits_Min = cms.uint32( 0 ),
          PropagateTracksToRadius = cms.bool( True ),
          ReferenceRadius = cms.double( 6.0 ),
          Chi2Prob_Min = cms.double( -1.0 ),
          Diff_r = cms.double( 0.1 ),
          BeamlineOption = cms.string( "BeamSpotFromEvent" ),
          ComponentName = cms.string( "PixelTrackExtractor" ),
          DR_VetoPt = cms.double( 0.025 )
        ),
        UseRhoCorrectedCaloDeposits = cms.bool( False ),
        UseCaloIso = cms.bool( False )
    )

    return process

def customizePhase2MuonHLTTrkIsolationOffline(process, processName = "MYHLT"):

    process.hltPhase2L3MuonsTrkIsoOfflineBase = cms.EDProducer( "L3MuonCombinedRelativeIsolationProducer",
        printDebug = cms.bool( False ),
        CutsPSet = cms.PSet( 
          applyCutsORmaxNTracks = cms.bool( False ),
          maxNTracks = cms.int32( -1 ),
          Thresholds = cms.vdouble( 1e9 ),
          EtaBounds = cms.vdouble( 2.411 ),
          ComponentName = cms.string( "SimpleCuts" ),
          ConeSizes = cms.vdouble( 1.0 )  ### 0.3
        ),
        OutputMuIsoDeposits = cms.bool( True ),
        TrackPt_Min = cms.double( -1.0 ),
        CaloDepositsLabel = cms.InputTag( "notUsed" ),
        CaloExtractorPSet = cms.PSet( 
          DR_Veto_H = cms.double( 0.1 ),
          Vertex_Constraint_Z = cms.bool( False ),
          DR_Veto_E = cms.double( 0.07 ),
          Weight_H = cms.double( 1.0 ),
          CaloTowerCollectionLabel = cms.InputTag( "hltTowerMakerForAll" ),
          DR_Max = cms.double( 0.3 ),
          DepositLabel = cms.untracked.string( "EcalPlusHcal" ),
          Vertex_Constraint_XY = cms.bool( False ),
          Threshold_H = cms.double( 0.5 ),
          Threshold_E = cms.double( 0.2 ),
          ComponentName = cms.string( "CaloExtractor" ),
          Weight_E = cms.double( 1.0 )
        ),
        inputMuonCollection = cms.InputTag( "hltPhase2L3MuonCandidates" ),
        TrkExtractorPSet = cms.PSet( 
          Diff_z = cms.double( 0.2 ),
          inputTrackCollection = cms.InputTag( "generalTracks" ),
          Chi2Ndof_Max = cms.double( 1.0E64 ),
          BeamSpotLabel = cms.InputTag( "offlineBeamSpot" ),
          DR_Veto = cms.double( 0.01 ),
          Pt_Min = cms.double( -1.0 ),
          VetoLeadingTrack = cms.bool( True ),
          DR_Max = cms.double( 1.0 ),  ### 0.3
          DepositLabel = cms.untracked.string( "PXLS" ),
          PtVeto_Min = cms.double( 2.0 ),
          NHits_Min = cms.uint32( 0 ),
          PropagateTracksToRadius = cms.bool( True ),
          ReferenceRadius = cms.double( 6.0 ),
          Chi2Prob_Min = cms.double( -1.0 ),
          Diff_r = cms.double( 0.1 ),
          BeamlineOption = cms.string( "BeamSpotFromEvent" ),
          ComponentName = cms.string( "PixelTrackExtractor" ),
          DR_VetoPt = cms.double( 0.025 )
        ),
        UseRhoCorrectedCaloDeposits = cms.bool( False ),
        UseCaloIso = cms.bool( False )
    )

    return process

def customizePhase2MuonHLTIsolation(process, processName = "MYHLT"):

    # -- Calo Isolations -- #
    process = customizePhase2MuonHLTCaloLocalReco(process)
    process = customizePhase2MuonHLTFixedGridRho(process)
    process = customizePhase2MuonHLTEcalIsolation(process)
    process = customizePhase2MuonHLTHcalIsolation(process)
    process = customizePhase2MuonHLTHgcalLCIsolation(process)

    process.hltPhase2L3MuonsEcalIsodR0p3dRVeto0p000 = cms.EDProducer("MuonHLTEcalPFClusterIsolationProducer",
        absEtaLowEdges = cms.vdouble(0.0, 1.479),
        doRhoCorrection = cms.bool(False),
        drMax = cms.double(0.3),
        drVetoBarrel = cms.double(0.0),
        drVetoEndcap = cms.double(0.0),
        effectiveAreas = cms.vdouble(0.35, 0.193),
        energyBarrel = cms.double(0.0),
        energyEndcap = cms.double(0.0),
        etaStripBarrel = cms.double(0.0),
        etaStripEndcap = cms.double(0.0),
        pfClusterProducer = cms.InputTag("hltParticleFlowClusterECALUnseeded"),
        recoCandidateProducer = cms.InputTag("hltPhase2L3MuonCandidates"),
        rhoMax = cms.double(99999999.0),
        rhoProducer = cms.InputTag("hltFixedGridRhoFastjetAllCaloForMuons"),
        rhoScale = cms.double(1.0)
    )

    process.hltPhase2L3MuonsHcalIsodR0p3dRVeto0p000 = cms.EDProducer("MuonHLTHcalPFClusterIsolationProducer",
        absEtaLowEdges = cms.vdouble(0.0, 1.479),
        doRhoCorrection = cms.bool(False),
        drMax = cms.double(0.3),
        drVetoBarrel = cms.double(0.0),
        drVetoEndcap = cms.double(0.0),
        effectiveAreas = cms.vdouble(0.227, 0.372),
        energyBarrel = cms.double(0.0),
        energyEndcap = cms.double(0.0),
        etaStripBarrel = cms.double(0.0),
        etaStripEndcap = cms.double(0.0),
        pfClusterProducerHCAL = cms.InputTag("hltParticleFlowClusterHCALForMuons"),
        pfClusterProducerHFEM = cms.InputTag(""),
        pfClusterProducerHFHAD = cms.InputTag(""),
        recoCandidateProducer = cms.InputTag("hltPhase2L3MuonCandidates"),
        rhoMax = cms.double(99999999.0),
        rhoProducer = cms.InputTag("hltFixedGridRhoFastjetAllCaloForMuons"),
        rhoScale = cms.double(1.0),
        useEt = cms.bool(True),
        useHF = cms.bool(False)
    )

    process.hltPhase2L3MuonsHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00 = cms.EDProducer("MuonHLTHGCalLayerClusterIsolationProducer",
        doRhoCorrection = cms.bool(False),
        drMax = cms.double(0.2),
        drVetoEM = cms.double(0.0),
        drVetoHad = cms.double(0.02),
        effectiveAreas = cms.vdouble(0.0, 0.0),
        layerClusterProducer = cms.InputTag("hgcalLayerClusters"),
        minEnergyEM = cms.double(0.0),
        minEnergyHad = cms.double(0.0),
        minEtEM = cms.double(0.0),
        minEtHad = cms.double(0.0),
        recoCandidateProducer = cms.InputTag("hltPhase2L3MuonCandidates"),
        rhoMax = cms.double(99999999.0),
        rhoProducer = cms.InputTag("hltFixedGridRhoFastjetAllCaloForMuons"),
        rhoScale = cms.double(1.0),
        useEt = cms.bool(False)
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

    process.HLTPhase2L3MuonCaloIsoSequenceIgnoreFilters = cms.Sequence(
        process.HLTDoFullUnpackingEgammaEcalSequence +
        process.HLTDoLocalHcalSequence +
        process.HLTFastJetForMuons +

        process.HLTPFClusteringForMuonsUnseeded +
        process.hltPhase2L3MuonsEcalIsodR0p3dRVeto0p000 +
        cms.ignore( process.hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41 ) +

        process.HLTPFHcalClusteringForMuons +
        process.hltPhase2L3MuonsHcalIsodR0p3dRVeto0p000 +
        cms.ignore( process.hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40 ) +

        process.HLTHgcalTiclLayerClusteringForMuons + 
        process.hltPhase2L3MuonsHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00 +
        cms.ignore( process.hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70 )
    )

    process.HLTPhase2L3MuonCaloIsoSequence = cms.Sequence(
        process.HLTDoFullUnpackingEgammaEcalSequence +
        process.HLTDoLocalHcalSequence +
        process.HLTFastJetForMuons +

        process.HLTPFClusteringForMuonsUnseeded +
        process.hltPhase2L3MuonsEcalIsodR0p3dRVeto0p000 +
        process.hltL3crIsoL1TkSingleMu22L3f24QL3pfecalIsoFiltered0p41 +

        process.HLTPFHcalClusteringForMuons +
        process.hltPhase2L3MuonsHcalIsodR0p3dRVeto0p000 +
        process.hltL3crIsoL1TkSingleMu22L3f24QL3pfhcalIsoFiltered0p40 +

        process.HLTHgcalTiclLayerClusteringForMuons + 
        process.hltPhase2L3MuonsHgcalLCIsodR0p2dRVetoEM0p00dRVetoHad0p02minEEM0p00minEHad0p00 +
        process.hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70
    )


    # -- Tracker Isolation -- #
    process = customizePhase2MuonHLTTrkIsolationRegionalNew(process)

    process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0 = cms.EDProducer("L3MuonCombinedRelativeIsolationProducer",
        CaloDepositsLabel = cms.InputTag("notUsed"),
        CaloExtractorPSet = cms.PSet(
            CaloTowerCollectionLabel = cms.InputTag("hltTowerMakerForAll"),
            ComponentName = cms.string('CaloExtractor'),
            DR_Max = cms.double(0.3),
            DR_Veto_E = cms.double(0.07),
            DR_Veto_H = cms.double(0.1),
            DepositLabel = cms.untracked.string('EcalPlusHcal'),
            Threshold_E = cms.double(0.2),
            Threshold_H = cms.double(0.5),
            Vertex_Constraint_XY = cms.bool(False),
            Vertex_Constraint_Z = cms.bool(False),
            Weight_E = cms.double(1.0),
            Weight_H = cms.double(1.0)
        ),
        CutsPSet = cms.PSet(
            ComponentName = cms.string('SimpleCuts'),
            ConeSizes = cms.vdouble(0.3),
            EtaBounds = cms.vdouble(2.411),
            Thresholds = cms.vdouble(1000000000.0),
            applyCutsORmaxNTracks = cms.bool(False),
            maxNTracks = cms.int32(-1)
        ),
        OutputMuIsoDeposits = cms.bool(True),
        TrackPt_Min = cms.double(-1.0),
        TrkExtractorPSet = cms.PSet(
            BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
            BeamlineOption = cms.string('BeamSpotFromEvent'),
            Chi2Ndof_Max = cms.double(1e+64),
            Chi2Prob_Min = cms.double(-1.0),
            ComponentName = cms.string('PixelTrackExtractor'),
            DR_Max = cms.double(0.3),
            DR_Veto = cms.double(0.005),
            DR_VetoPt = cms.double(0.025),
            DepositLabel = cms.untracked.string('PXLS'),
            Diff_r = cms.double(0.2),
            Diff_z = cms.double(0.25),
            NHits_Min = cms.uint32(0),
            PropagateTracksToRadius = cms.bool(True),
            PtVeto_Min = cms.double(2.0),
            Pt_Min = cms.double(-1.0),
            ReferenceRadius = cms.double(6.0),
            VetoLeadingTrack = cms.bool(True),
            inputTrackCollection = cms.InputTag("hltPhase2L3MuonGeneralTracks")
        ),
        UseCaloIso = cms.bool(False),
        UseRhoCorrectedCaloDeposits = cms.bool(False),
        inputMuonCollection = cms.InputTag("hltPhase2L3MuonCandidates"),
        printDebug = cms.bool(False)
    )

    process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p07 = process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0.clone()
    process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p07.CutsPSet.Thresholds = cms.vdouble( 0.07 )

    process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p4 = process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0.clone()
    process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p4.CutsPSet.Thresholds = cms.vdouble( 0.4 )


    process.hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk = cms.EDFilter( "HLTMuonIsoFilter",
        saveTags = cms.bool( True ),
        PreviousCandTag = cms.InputTag( "hltL3crIsoL1TkSingleMu22L3f24QL3pfhgcalIsoFiltered4p70" ),
        MinN = cms.int32( 1 ),
        IsolatorPSet = cms.PSet(  ),
        CandTag = cms.InputTag( "hltPhase2L3MuonCandidates" ),
        DepTag = cms.VInputTag( 'hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p07' )
    )

    process.HLTPhase2L3MuonIsoSequenceTrkRegionalNewIgnoreFilters = cms.Sequence(
        process.HLTPhase2L3MuonCaloIsoSequenceIgnoreFilters + 

        process.tracking_v6_1_L3Muon_NoVertexReco +
        process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p07 +
        cms.ignore( process.hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk )
    )

    process.HLTPhase2L3MuonIsoSequenceTrkRegionalNew = cms.Sequence(
        process.HLTPhase2L3MuonCaloIsoSequence + 

        process.tracking_v6_1_L3Muon_NoVertexReco +
        process.hltPhase2L3MuonsTrkIsoRegionalNewdR0p3dRVeto0p005dz0p25dr0p20ChisqInfPtMin0p0Cut0p07 +
        process.hltL3crIsoL1TkSingleMu22L3f24QL3trkIsoRegionalNewFiltered0p07EcalHcalHgcalTrk
    )

    return process

def customizePhase2MuonHLTIsolationForOpt(process, processName = "MYHLT"):

    process = customizePhase2MuonHLTCaloLocalReco(process)
    process = customizePhase2MuonHLTFixedGridRho(process)
    process = customizePhase2MuonHLTEcalIsolation(process)
    process = customizePhase2MuonHLTHcalIsolation(process)
    process = customizePhase2MuonHLTHgcalLCIsolation(process)
    process = customizePhase2MuonHLTTrkIsolationRegional(process)
    process = customizePhase2MuonHLTTrkIsolationRegionalNew(process)
    process = customizePhase2MuonHLTTrkIsolationFull(process)
    process = customizePhase2MuonHLTTrkIsolationOffline(process)

    process.HLTPhase2L3MuonBaseIsoSequence = cms.Sequence(

        process.HLTDoFullUnpackingEgammaEcalSequence +
        process.HLTDoLocalHcalSequence +
        process.HLTFastJetForMuons +

        process.HLTPFClusteringForMuonsUnseeded +
        process.hltPhase2L3MuonsEcalIsoBase +

        process.HLTPFHcalClusteringForMuons +
        process.hltPhase2L3MuonsHcalIsoBase +

        process.HLTHgcalTiclLayerClusteringForMuons +
        process.hltPhase2L3MuonsHgcalLCIsoBase +

        process.HLTTrackReconstructionForIsoPhase2L3MuonIter02 +
        process.hltPhase2L3MuonsTrkIsoBase +

        process.tracking_v6_1_NoVertexReco +
        process.hltPhase2L3MuonsTrkIsoFullBase+

        process.tracking_v6_1_L3Muon_NoVertexReco +
        process.hltPhase2L3MuonsTrkIsoRegionalNewBase +

        process.hltPhase2L3MuonsTrkIsoOfflineBase
    )


    pfIsodRs     = [ 0.3 ]  #, 0.4 ]  # , 0.5 ]
    pfIsodRVetos = [ 0.00, 0.03, 0.05 ]  # , 0.07, 0.1 ]

    pfIsoTags = []
    pfIsoLabels = []
    pfIsoMods = []
    for pfIsodR in pfIsodRs:
        for pfIsodRVeto in pfIsodRVetos:

            EcalIsoTag = "EcalIsodR%.1fdRVeto%.3f" % (pfIsodR, pfIsodRVeto)
            EcalIsoTag = EcalIsoTag.replace(".", "p")
            EcalIsoModName = "hltPhase2L3Muons"+EcalIsoTag
            EcalIsoMod = process.hltPhase2L3MuonsEcalIsoBase.clone(
                drMax = cms.double( pfIsodR ),
                drVetoBarrel = cms.double( pfIsodRVeto ),
                drVetoEndcap = cms.double( pfIsodRVeto )
            )

            setattr(process, EcalIsoModName, EcalIsoMod)
            pfIsoTags.append( "pf"+EcalIsoTag )
            pfIsoLabels.append( EcalIsoModName )
            pfIsoMods.append( EcalIsoMod )

            HcalIsoTag  = "HcalIsodR%.1fdRVeto%.3f" % (pfIsodR, pfIsodRVeto)
            HcalIsoTag = HcalIsoTag.replace(".", "p")
            HcalIsoModName = "hltPhase2L3Muons"+HcalIsoTag
            HcalIsoMod = process.hltPhase2L3MuonsHcalIsoBase.clone(
                drMax = cms.double( pfIsodR ),
                drVetoBarrel = cms.double( pfIsodRVeto ),
                drVetoEndcap = cms.double( pfIsodRVeto )
            )

            setattr(process, HcalIsoModName, HcalIsoMod)
            pfIsoTags.append( "pf"+HcalIsoTag )
            pfIsoLabels.append( HcalIsoModName )
            pfIsoMods.append( HcalIsoMod )

            # HgcalIsoTag = "HgcalIsodR%.1fdRVeto%.3f" % (pfIsodR, pfIsodRVeto)
            # HgcalIsoTag = HgcalIsoTag.replace(".", "p")
            # HgcalIsoModName = "hltPhase2L3Muons"+HgcalIsoTag
            # HgcalIsoMod = process.hltPhase2L3MuonsHgcalIsoBase.clone(
            #     drMax = cms.double( pfIsodR ),
            #     drVetoBarrel = cms.double( pfIsodRVeto ),
            #     drVetoEndcap = cms.double( pfIsodRVeto )
            # )

            # setattr(process, HgcalIsoModName, HgcalIsoMod)
            # pfIsoTags.append( "pf"+HgcalIsoTag )
            # pfIsoLabels.append( HgcalIsoModName )
            # pfIsoMods.append( HgcalIsoMod )


    lcIsodRs        = [ 0.2 ]  # , 0.3, 0.4 ]
    lcIsodRVetoEMs  = [ 0.00, 0.02, 0.04 ]
    lcIsodRVetoHads = [ 0.00, 0.02, 0.04 ]
    lcIsoMinEEMs    = [ 0.00 ]  # , 0.02, 0.05, 0.07 ]
    lcIsoMinEHads   = [ 0.00 ]  # , 0.02, 0.05, 0.07 ]
    for lcIsodR in lcIsodRs:
        for lcIsodRVetoEM in lcIsodRVetoEMs:
            for lcIsodRVetoHad in lcIsodRVetoHads:
                for lcIsoMinEEM in lcIsoMinEEMs:
                    for lcIsoMinEHad in lcIsoMinEHads:
                        HgcalLCIsoTag = "HgcalLCIsodR%.1fdRVetoEM%.2fdRVetoHad%.2fminEEM%.2fminEHad%.2f" % (lcIsodR, lcIsodRVetoEM, lcIsodRVetoHad, lcIsoMinEEM, lcIsoMinEHad)
                        HgcalLCIsoTag = HgcalLCIsoTag.replace(".", "p")
                        HgcalLCIsoModName = "hltPhase2L3Muons"+HgcalLCIsoTag
                        HgcalLCIsoMod = process.hltPhase2L3MuonsHgcalLCIsoBase.clone(
                            drMax = cms.double(lcIsodR),
                            drVetoEM = cms.double(lcIsodRVetoEM),
                            drVetoHad = cms.double(lcIsodRVetoHad),
                            minEnergyEM = cms.double(lcIsoMinEEM),
                            minEnergyHad = cms.double(lcIsoMinEHad)
                        )

                        setattr(process, HgcalLCIsoModName, HgcalLCIsoMod)
                        pfIsoTags.append( "pf"+HgcalLCIsoTag )
                        pfIsoLabels.append( HgcalLCIsoModName )
                        pfIsoMods.append( HgcalLCIsoMod )

    pfIsoSeq = reduce(lambda x,y: x+y, pfIsoMods)

    trkIsodRs     = [ 0.3 ]  # , 0.4 ]
    trkIsodRVetos = [ 0.005 ]  # , 0.01 ]  # , 0.02 ]
    trkIsodzs     = [ 0.1, 0.2, 0.25 ]  # [ 0.1, 0.15, 0.2, 0.25 ]
    trkIsodrs     = [ 0.1, 0.2 ]
    trkIsoChi2s   = [ 1.0E64 ]
    trkIsoPtMins  = [ -1.0 ]
    trkIsoTypes   = [ "Regional", "Full", "RegionalNew", "Offline" ]

    trkIsoTags = []
    trkIsoLabels = []
    trkIsoMods = []
    for trkIsodR in trkIsodRs:
        for trkIsodRVeto in trkIsodRVetos:
            for trkIsodz in trkIsodzs:
                for trkIsodr in trkIsodrs:
                    for trkIsoChi2 in trkIsoChi2s:
                        for trkIsoPtMin in trkIsoPtMins:
                            for trkIsoType in trkIsoTypes:

                                trkIsoTag = "trkIso%sdR%.1fdRVeto%.3fdz%.2fdr%.2fChisq%.0fPtMin%.1f" % (trkIsoType, trkIsodR, trkIsodRVeto, trkIsodz, trkIsodr, trkIsoChi2, trkIsoPtMin)
                                trkIsoTag = trkIsoTag.replace(".", "p")
                                trkIsoTag = trkIsoTag.replace("Chisq%.0f"%1.0E64, "ChisqInf")
                                trkIsoTag = trkIsoTag.replace("PtMin-1p0", "PtMin0p0")

                                trkIsoModName = "hltPhase2L3Muons"+trkIsoTag
                                trkIsoModName = trkIsoModName.replace("trkIso", "TrkIso")

                                if trkIsoType == "Regional":
                                    trkIsoMod = process.hltPhase2L3MuonsTrkIsoBase.clone()
                                elif trkIsoType == "Full":
                                    trkIsoMod = process.hltPhase2L3MuonsTrkIsoFullBase.clone()
                                elif trkIsoType == "RegionalNew":
                                    trkIsoMod = process.hltPhase2L3MuonsTrkIsoRegionalNewBase.clone()
                                elif trkIsoType == "Offline":
                                    trkIsoMod = process.hltPhase2L3MuonsTrkIsoOfflineBase.clone()
                                else:
                                    trkIsoMod = None

                                trkIsoMod.CutsPSet.ConeSizes            = cms.vdouble( trkIsodR )
                                trkIsoMod.TrkExtractorPSet.DR_Max       = cms.double( trkIsodR )
                                trkIsoMod.TrkExtractorPSet.DR_Veto      = cms.double( trkIsodRVeto )
                                trkIsoMod.TrkExtractorPSet.Diff_z       = cms.double( trkIsodz )
                                trkIsoMod.TrkExtractorPSet.Diff_r       = cms.double( trkIsodr )
                                trkIsoMod.TrkExtractorPSet.Chi2Ndof_Max = cms.double( trkIsoChi2 )
                                trkIsoMod.TrkExtractorPSet.Pt_Min       = cms.double( trkIsoPtMin )

                                setattr(process, trkIsoModName, trkIsoMod)
                                trkIsoTags.append( trkIsoTag )
                                trkIsoLabels.append( trkIsoModName+":trkIsoDeposits" )
                                trkIsoMods.append( trkIsoMod )

    trkIsoSeq = reduce(lambda x,y: x+y, trkIsoMods)

    process.IsolationStudySequence = cms.Sequence(
        pfIsoSeq+
        trkIsoSeq
    )

    return process, (pfIsoTags, pfIsoLabels, pfIsoMods), (trkIsoTags, trkIsoLabels, trkIsoMods)



