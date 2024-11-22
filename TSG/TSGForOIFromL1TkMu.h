#ifndef RecoMuon_TrackerSeedGenerator_TSGForOIFromL1TkMu_H
#define RecoMuon_TrackerSeedGenerator_TSGForOIFromL1TkMu_H

/**
 \class    TSGForOIFromL1TkMu
 \brief    Create L3MuonTrajectorySeeds from L1 Tracker Muons in an outside-in manner
  \author   Arnab Purohit, Dmitry Kondratyev (Purdue University, West Lafayette, USA)
 */

#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTrackerEvent.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/MeasurementDet/interface/MeasurementDet.h"
#include "TrackingTools/PatternTools/interface/TrajMeasLessEstim.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/StateOnTrackerBound.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "DataFormats/L1TCorrelator/interface/TkMuon.h"
#include "DataFormats/L1TCorrelator/interface/TkMuonFwd.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"

class TSGForOIFromL1TkMu : public edm::global::EDProducer<> {
public:
  explicit TSGForOIFromL1TkMu(const edm::ParameterSet& iConfig);
  ~TSGForOIFromL1TkMu() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void produce(edm::StreamID sid, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

private:

  const edm::EDGetTokenT<l1t::TrackerMuonCollection> src_;
//const edm::EDGetTokenT<l1t::TkMuonCollection> src_;

  /// Minimum Pt of the L1TkMuons to remove soft fakes and as also the L1TkMuon pt resolution is very close to the L3 hence we can apply a cut on pt close to the pt cut we apply on the L3 muons.
  const double minPtOfL1TKMuons_;

  /// Maximum number of seeds for each L1 TkMu
  const unsigned int maxSeeds_;

  /// Maximum number of hitless seeds for each L1 TkMu
  const unsigned int maxHitlessSeeds_;

  /// Maximum number of hitbased seeds for each L1 TkMu
  const unsigned int maxHitSeeds_;

  /// How many layers to try
  const unsigned int numOfLayersToTry_;

  /// How many hits to try per layer
  const unsigned int numOfHitsToTry_;

  /// Rescale L1 TkMu parameter uncertainties (fixed error vs pT, eta)
  const double fixedErrorRescalingForHits_;
  const double fixedErrorRescalingForHitless_;

  /// Whether or not to use an automatically calculated scale-factor value
  const bool adjustErrorsDynamicallyForHits_;
  const bool adjustErrorsDynamicallyForHitless_;

  /// Estimator used to find dets and TrajectoryMeasurements
  const std::string estimatorName_;

  /// Minimum eta value to activate searching in the TEC
  const double minEtaForTEC_;

  /// Maximum eta value to activate searching in the TOB
  const double maxEtaForTOB_;

  /// Switch ON  (True) : use additional hits for seeds depending on the L1 TkMu properties (ignores numOfMaxSeeds_)
  /// Switch OFF (False): the numOfMaxSeeds_ defines if we will use hitless (numOfMaxSeeds_==1) or hitless+hits (numOfMaxSeeds_>1)
  const bool useHitLessSeeds_;

  /// KFUpdator defined in constructor
  const std::unique_ptr<TrajectoryStateUpdator> updator_;

  const edm::EDGetTokenT<MeasurementTrackerEvent> measurementTrackerTag_;

  /// pT, eta ranges and scale factor values
  const double pT1_, pT2_, pT3_;
  const double eta1_, eta2_, eta3_, eta4_, eta5_, eta6_, eta7_;
  const double SF1_, SF2_, SF3_, SF4_, SF5_, SF6_;

  /// Distance of L1 TkMu TSOSs before and after updated with vertex
  const double tsosDiff1_;
  const double tsosDiff2_;

  /// Counters and flags for the implementation
  const std::string propagatorName_;
  const std::string theCategory_;


  const edm::ESGetToken<Chi2MeasurementEstimatorBase, TrackingComponentsRecord> estimatorToken_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magfieldToken_;
  const edm::ESGetToken<Propagator, TrackingComponentsRecord> propagatorToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tmpTkGeometryToken_;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geometryToken_;
  const edm::ESGetToken<Propagator, TrackingComponentsRecord> sHPOppositeToken_;






  /// Create seeds without hits on a given layer (TOB or TEC)
  void makeSeedsWithoutHits(const GeometricSearchDet& layer,
                            const TrajectoryStateOnSurface& tsos,
                            const Propagator& propagatorAlong,
                            const edm::ESHandle<Chi2MeasurementEstimatorBase>& estimator,
                            double errorSF,
                            unsigned int& hitlessSeedsMade,
                            unsigned int& numSeedsMade,
                            std::vector<TrajectorySeed>& out) const;

  /// Find hits on a given layer (TOB or TEC) and create seeds from updated TSOS with hit
  void makeSeedsFromHits(const GeometricSearchDet& layer,
                         const TrajectoryStateOnSurface& tsos,
                         const Propagator& propagatorAlong,
                         const edm::ESHandle<Chi2MeasurementEstimatorBase>& estimator,
                         const edm::Handle<MeasurementTrackerEvent>& measurementTracker,
                         double errorSF,
                         unsigned int& hitSeedsMade,
                         unsigned int& numSeedsMade,
                         unsigned int& layerCount,
                         std::vector<TrajectorySeed>& out) const;

  /// Calculate the dynamic error SF by analysing the L1 TkMu
  double calculateSFFromL1TkMu(const TTTrack<Ref_Phase2TrackerDigi_>& track) const;

  /// Find compatability between two TSOSs
  double match_Chi2(const TrajectoryStateOnSurface& tsos1, const TrajectoryStateOnSurface& tsos2) const;

  FreeTrajectoryState initialFreeStateL1TTrack(const TTTrack<Ref_Phase2TrackerDigi_>& tk,
                                               const MagneticField* field,
                                               bool withErr = false) const;

};

#endif
