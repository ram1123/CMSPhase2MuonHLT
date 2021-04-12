/**
  \class    TSGForOIFromL1TkMu
  \brief    Create L3MuonTrajectorySeeds from L1 Tracker Muons in an outside-in manner
  \author   Arnab Purohit (Purdue University, West Lafayette, USA)
 */

#include "RecoMuon/TrackerSeedGenerator/plugins/TSGForOIFromL1TkMu.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/L1TCorrelator/interface/TkMuon.h"
#include "DataFormats/L1TCorrelator/interface/TkMuonFwd.h"

#include <memory>

TSGForOIFromL1TkMu::TSGForOIFromL1TkMu(const edm::ParameterSet& iConfig)
  : src_(consumes<l1t::TkMuonCollection>(iConfig.getParameter<edm::InputTag>("src"))),
    minPtOfL1TKMuons_(iConfig.getParameter<double>("minPtOfL1TKMuons")),
    maxSeeds_(iConfig.getParameter<uint32_t>("maxSeeds")),
      maxHitlessSeeds_(iConfig.getParameter<uint32_t>("maxHitlessSeeds")),
      maxHitSeeds_(iConfig.getParameter<uint32_t>("maxHitSeeds")),
      numOfLayersToTry_(iConfig.getParameter<int32_t>("layersToTry")),
      numOfHitsToTry_(iConfig.getParameter<int32_t>("hitsToTry")),
      fixedErrorRescalingForHits_(iConfig.getParameter<double>("fixedErrorRescaleFactorForHits")),
      fixedErrorRescalingForHitless_(iConfig.getParameter<double>("fixedErrorRescaleFactorForHitless")),
      adjustErrorsDynamicallyForHits_(iConfig.getParameter<bool>("adjustErrorsDynamicallyForHits")),
      adjustErrorsDynamicallyForHitless_(iConfig.getParameter<bool>("adjustErrorsDynamicallyForHitless")),
      estimatorName_(iConfig.getParameter<std::string>("estimator")),
      minEtaForTEC_(iConfig.getParameter<double>("minEtaForTEC")),
      maxEtaForTOB_(iConfig.getParameter<double>("maxEtaForTOB")),
      useHitLessSeeds_(iConfig.getParameter<bool>("UseHitLessSeeds")),
      updator_(new KFUpdator()),
      measurementTrackerTag_(
          consumes<MeasurementTrackerEvent>(iConfig.getParameter<edm::InputTag>("MeasurementTrackerEvent"))),
      pT1_(iConfig.getParameter<double>("pT1")),
      pT2_(iConfig.getParameter<double>("pT2")),
      pT3_(iConfig.getParameter<double>("pT3")),
      eta1_(iConfig.getParameter<double>("eta1")),
      eta2_(iConfig.getParameter<double>("eta2")),
      eta3_(iConfig.getParameter<double>("eta3")),
      eta4_(iConfig.getParameter<double>("eta4")),
      eta5_(iConfig.getParameter<double>("eta5")),
      eta6_(iConfig.getParameter<double>("eta6")),
      eta7_(iConfig.getParameter<double>("eta7")),
      SF1_(iConfig.getParameter<double>("SF1")),
      SF2_(iConfig.getParameter<double>("SF2")),
      SF3_(iConfig.getParameter<double>("SF3")),
      SF4_(iConfig.getParameter<double>("SF4")),
      SF5_(iConfig.getParameter<double>("SF5")),
      SF6_(iConfig.getParameter<double>("SF6")),
      tsosDiff1_(iConfig.getParameter<double>("tsosDiff1")),
      tsosDiff2_(iConfig.getParameter<double>("tsosDiff2")),
      propagatorName_(iConfig.getParameter<std::string>("propagatorName")),
      theCategory_(std::string("Muon|RecoMuon|TSGForOIFromL1TkMu")) {
  produces<std::vector<TrajectorySeed> >();
}

TSGForOIFromL1TkMu::~TSGForOIFromL1TkMu() {}

//
// Produce seeds
//
void TSGForOIFromL1TkMu::produce(edm::StreamID sid, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  //  std::cout<<"Event number is "<<iEvent.id().event()<<std::endl;

  // Initialize variables
  unsigned int numSeedsMade = 0;
  unsigned int layerCount = 0;
  unsigned int hitlessSeedsMadeIP = 0;
  unsigned int hitSeedsMade = 0;
  int l1tkmu_idx = 0;
  // Surface used to make a TSOS at the PCA to the beamline
  Plane::PlanePointer dummyPlane = Plane::build(Plane::PositionType(), Plane::RotationType());

  // Read ESHandles
  edm::Handle<MeasurementTrackerEvent> measurementTrackerH;
  edm::ESHandle<Chi2MeasurementEstimatorBase> estimatorH;
  edm::ESHandle<MagneticField> magfieldH;
  edm::ESHandle<Propagator> propagatorAlongH;
  edm::ESHandle<Propagator> propagatorOppositeH;
  edm::ESHandle<TrackerGeometry> tmpTkGeometryH;
  edm::ESHandle<GlobalTrackingGeometry> geometryH;

  iSetup.get<IdealMagneticFieldRecord>().get(magfieldH);
  iSetup.get<TrackingComponentsRecord>().get(propagatorName_, propagatorOppositeH);
  iSetup.get<TrackingComponentsRecord>().get(propagatorName_, propagatorAlongH);
  iSetup.get<GlobalTrackingGeometryRecord>().get(geometryH);
  iSetup.get<TrackerDigiGeometryRecord>().get(tmpTkGeometryH);
  iSetup.get<TrackingComponentsRecord>().get(estimatorName_, estimatorH);
  iEvent.getByToken(measurementTrackerTag_, measurementTrackerH);

  // Read L1 TkMuon collection
  edm::Handle<l1t::TkMuonCollection> l1TkMuCol;

  iEvent.getByToken(src_, l1TkMuCol); 

  // The product
  std::unique_ptr<std::vector<TrajectorySeed> > result(new std::vector<TrajectorySeed>());

  // Get vector of Detector layers
  std::vector<BarrelDetLayer const*> const& tob = measurementTrackerH->geometricSearchTracker()->tobLayers();
  std::vector<ForwardDetLayer const*> const& tecPositive =
      tmpTkGeometryH->isThere(GeomDetEnumerators::P2OTEC)
          ? measurementTrackerH->geometricSearchTracker()->posTidLayers()
          : measurementTrackerH->geometricSearchTracker()->posTecLayers();
  std::vector<ForwardDetLayer const*> const& tecNegative =
      tmpTkGeometryH->isThere(GeomDetEnumerators::P2OTEC)
          ? measurementTrackerH->geometricSearchTracker()->negTidLayers()
          : measurementTrackerH->geometricSearchTracker()->negTecLayers();

  // Get suitable propagators
  std::unique_ptr<Propagator> propagatorAlong = SetPropagationDirection(*propagatorAlongH, alongMomentum);
  std::unique_ptr<Propagator> propagatorOpposite = SetPropagationDirection(*propagatorOppositeH, oppositeToMomentum);

  // Stepping Helix Propagator for propogation from muon system to tracker
  edm::ESHandle<Propagator> SHPOpposite;
  iSetup.get<TrackingComponentsRecord>().get("hltESPSteppingHelixPropagatorOpposite", SHPOpposite);

  // Loop over the L1 Tracker Muons and make seeds for all of them
  LogTrace(theCategory_) << "TSGForOIFromL1TkMu::produce: Number of L1TkMu's: " << l1TkMuCol->size();
  for (auto ittkmu = l1TkMuCol->begin(); ittkmu != l1TkMuCol->end(); ittkmu++) {

    auto l1TkMu_trkPtr = ittkmu->trkPtr();

    // Container of Seeds
    std::vector<TrajectorySeed> out;

    auto p3 = l1TkMu_trkPtr->momentum();

    LogTrace("TSGForOIFromL1TkMu") << "TSGForOIFromL1TkMu::produce: L1TKMU_TRKPTR muon pT, eta, phi --> " << p3.perp() << " , " << p3.eta() << " , " << p3.phi() << std::endl;
    if(p3.perp()<minPtOfL1TKMuons_) continue;

    FreeTrajectoryState fts = initialFreeStateL1TTrack(*l1TkMu_trkPtr, magfieldH.product(), true);

    dummyPlane->move(fts.position() - dummyPlane->position());
    TrajectoryStateOnSurface tsosAtPCA = TrajectoryStateOnSurface(fts, *dummyPlane);
    LogTrace("TSGForOIFromL1TkMu") << "TSGForOIFromL1TkMu::produce: Created TSOSatPCA: " << tsosAtPCA << std::endl;

    LogTrace("TSGForOIFromL1TkMu")
        << "TSGForOIFromL1TkMu::produce: Check the error of the L1TkMu_trkPtr parameter and use hit seeds if big errors" << std::endl;

    float tk_eta = p3.eta();
    float tk_aeta = std::abs(tk_eta);

    numSeedsMade = 0;
    hitlessSeedsMadeIP = 0;
    hitSeedsMade = 0;

    // calculate scale factors
    double errorSFHits = 
        (adjustErrorsDynamicallyForHits_ ? calculateSFFromL1TkMu(*l1TkMu_trkPtr) : fixedErrorRescalingForHits_);
    double errorSFHitless =
        (adjustErrorsDynamicallyForHitless_ ? calculateSFFromL1TkMu(*l1TkMu_trkPtr) : fixedErrorRescalingForHitless_);

    // BARREL
    if (tk_aeta < maxEtaForTOB_) {
      layerCount = 0;
      for (auto it = tob.rbegin(); it != tob.rend(); ++it) {
        LogTrace("TSGForOIFromL1TkMu") << "TSGForOIFromL1TkMu::produce: looping in TOB layer " << layerCount << std::endl;
        if (useHitLessSeeds_ && hitlessSeedsMadeIP < maxHitlessSeeds_ && numSeedsMade < maxSeeds_){
          makeSeedsWithoutHits(**it,
                               tsosAtPCA,
                               *(propagatorAlong.get()),
                               estimatorH,
                               errorSFHitless,
                               hitlessSeedsMadeIP,
                               numSeedsMade,
                               out);
        }

        if (hitSeedsMade < maxHitSeeds_ && numSeedsMade < maxSeeds_){
          makeSeedsFromHits(**it,
                            tsosAtPCA,
                            *(propagatorAlong.get()),
                            estimatorH,
                            measurementTrackerH,
                            errorSFHits,
                            hitSeedsMade,
                            numSeedsMade,
                            layerCount,
                            out);
        }

      }
      LogTrace("TSGForOIFromL1TkMu") << "TSGForOIFromL1TkMu:::produce: NumSeedsMade = " << numSeedsMade
                      << " , layerCount = " << layerCount << std::endl;
    }

    // Reset number of seeds if in overlap region
    if (tk_aeta > minEtaForTEC_ && tk_aeta < maxEtaForTOB_) {
      numSeedsMade = 0;
      hitlessSeedsMadeIP = 0;
      hitSeedsMade = 0;
    }

    // ENDCAP+
    if (tk_eta > minEtaForTEC_) {
      layerCount = 0;
      for (auto it = tecPositive.rbegin(); it != tecPositive.rend(); ++it) {
        LogTrace("TSGForOIFromL1TkMu") << "TSGForOIFromL1TkMu::produce: looping in TEC+ layer " << layerCount << std::endl;
        if (useHitLessSeeds_ && hitlessSeedsMadeIP < maxHitlessSeeds_ && numSeedsMade < maxSeeds_){
          makeSeedsWithoutHits(**it,
                               tsosAtPCA,
                               *(propagatorAlong.get()),
                               estimatorH,
                               errorSFHitless,
                               hitlessSeedsMadeIP,
                               numSeedsMade,
                               out);
        }

        if (hitSeedsMade < maxHitSeeds_ && numSeedsMade < maxSeeds_){
          makeSeedsFromHits(**it,
                            tsosAtPCA,
                            *(propagatorAlong.get()),
                            estimatorH,
                            measurementTrackerH,
                            errorSFHits,
                            hitSeedsMade,
                            numSeedsMade,
                            layerCount,
                            out);
        }

      }
      LogTrace("TSGForOIFromL1TkMu") << "TSGForOIFromL1TkMu:::produce: NumSeedsMade = " << numSeedsMade
                                 << " , layerCount = " << layerCount << std::endl;
    }

    // ENDCAP-
    if (tk_eta < -minEtaForTEC_) {
      layerCount = 0;
      for (auto it = tecNegative.rbegin(); it != tecNegative.rend(); ++it) {
        LogTrace("TSGForOIFromL1TkMu") << "TSGForOIFromL1TkMu::produce: looping in TEC- layer " << layerCount << std::endl;
        if (useHitLessSeeds_ && hitlessSeedsMadeIP < maxHitlessSeeds_ && numSeedsMade < maxSeeds_){
          makeSeedsWithoutHits(**it,
                               tsosAtPCA,
                               *(propagatorAlong.get()),
                               estimatorH,
                               errorSFHitless,
                               hitlessSeedsMadeIP,
                               numSeedsMade,
                               out);
        }

        if (hitSeedsMade < maxHitSeeds_ && numSeedsMade < maxSeeds_){
          makeSeedsFromHits(**it,
                            tsosAtPCA,
                            *(propagatorAlong.get()),
                            estimatorH,
                            measurementTrackerH,
                            errorSFHits,
                            hitSeedsMade,
                            numSeedsMade,
                            layerCount,
                            out);
        }

      }
      LogTrace("TSGForOIFromL1TkMu") << "TSGForOIFromL1TkMu:::produce: NumSeedsMade = " << numSeedsMade
                                 << " , layerCount = " << layerCount << std::endl;
    }

    for (std::vector<TrajectorySeed>::iterator it = out.begin(); it != out.end(); ++it) {
      result->push_back(*it);
    }
    l1tkmu_idx++;  
  }  // L1TkMu Collection

  edm::LogInfo(theCategory_) << "TSGForOIFromL1TkMu::produce: number of seeds made: " << result->size();

  iEvent.put(std::move(result));
}

//
// Create seeds without hits on a given layer (TOB or TEC)
//
void TSGForOIFromL1TkMu::makeSeedsWithoutHits(const GeometricSearchDet& layer,
                                              const TrajectoryStateOnSurface& tsos,
                                              const Propagator& propagatorAlong,
                                              edm::ESHandle<Chi2MeasurementEstimatorBase>& estimator,
                                              double errorSF,
                                              unsigned int& hitlessSeedsMade,
                                              unsigned int& numSeedsMade,
                                              std::vector<TrajectorySeed>& out) const {
  // create hitless seeds
  LogTrace("TSGForOIFromL1TkMu") << "TSGForOIFromL1TkMu::makeSeedsWithoutHits: Start hitless" << std::endl;
  std::vector<GeometricSearchDet::DetWithState> dets;
  layer.compatibleDetsV(tsos, propagatorAlong, *estimator, dets);
  if (!dets.empty()) {
    auto const& detOnLayer = dets.front().first;
    auto const& tsosOnLayer = dets.front().second;
    LogTrace("TSGForOIFromL1TkMu") << "TSGForOIFromL1TkMu::makeSeedsWithoutHits: tsosOnLayer " << tsosOnLayer << std::endl;
    if (!tsosOnLayer.isValid()) {
      edm::LogInfo(theCategory_) << "ERROR!: Hitless TSOS is not valid!";
    } else {
      dets.front().second.rescaleError(errorSF);
      PTrajectoryStateOnDet const& ptsod =
          trajectoryStateTransform::persistentState(tsosOnLayer, detOnLayer->geographicalId().rawId());
      TrajectorySeed::RecHitContainer rHC;
      out.push_back(TrajectorySeed(ptsod, rHC, oppositeToMomentum));
      hitlessSeedsMade++;
      numSeedsMade++;
      LogTrace("TSGForOIFromL1TkMu") << "TSGForOIFromL1TkMu::makeSeedsWithoutHits: number of hitless seeds found: " << hitlessSeedsMade<<" and number of seeds made: "<<numSeedsMade<<std::endl;
    }
  }
}

//
// Find hits on a given layer (TOB or TEC) and create seeds from updated TSOS with hit
//
void TSGForOIFromL1TkMu::makeSeedsFromHits(const GeometricSearchDet& layer,
                                           const TrajectoryStateOnSurface& tsos,
                                           const Propagator& propagatorAlong,
                                           edm::ESHandle<Chi2MeasurementEstimatorBase>& estimator,
                                           edm::Handle<MeasurementTrackerEvent>& measurementTracker,
                                           double errorSF,
                                           unsigned int& hitSeedsMade,
                                           unsigned int& numSeedsMade,
                                           unsigned int& layerCount,
                                           std::vector<TrajectorySeed>& out) const {
  if (layerCount > numOfLayersToTry_)
    return;
  //else std::cout<<"layerCount has reached numOfLayersToTry_"<<std::endl;
  // Error Rescaling
  TrajectoryStateOnSurface onLayer(tsos);
  onLayer.rescaleError(errorSF);

  std::vector<GeometricSearchDet::DetWithState> dets;
  layer.compatibleDetsV(onLayer, propagatorAlong, *estimator, dets);

  // Find Measurements on each DetWithState
  LogTrace("TSGForOIFromL1TkMu") << "TSGForOIFromL1TkMu::makeSeedsFromHits: Find measurements on each detWithState  " << dets.size() << std::endl;
  std::vector<TrajectoryMeasurement> meas;
  for (std::vector<GeometricSearchDet::DetWithState>::iterator it = dets.begin(); it != dets.end(); ++it) {
    MeasurementDetWithData det = measurementTracker->idToDet(it->first->geographicalId());
    if (det.isNull())
      continue;
    if (!it->second.isValid())
      continue;  // Skip if TSOS is not valid

    std::vector<TrajectoryMeasurement> mymeas =
        det.fastMeasurements(it->second, onLayer, propagatorAlong, *estimator);  // Second TSOS is not used
    for (std::vector<TrajectoryMeasurement>::const_iterator it2 = mymeas.begin(), ed2 = mymeas.end(); it2 != ed2;
         ++it2) {
      if (it2->recHit()->isValid())
        meas.push_back(*it2);  // Only save those which are valid
    }
  }

  // Update TSOS using TMs after sorting, then create Trajectory Seed and put into vector
  //LogTrace("TSGForOIFromL1TkMu") << "TSGForOIFromL1TkMu::makeSeedsFromHits: Update TSOS using TMs after sorting, then create "
  //                              "Trajectory Seed, number of TM = "
  //                           << meas.size() << std::endl;
  std::sort(meas.begin(), meas.end(), TrajMeasLessEstim());

  unsigned int found = 0;
  for (std::vector<TrajectoryMeasurement>::const_iterator it = meas.begin(); it != meas.end(); ++it) {
    TrajectoryStateOnSurface updatedTSOS = updator_->update(it->forwardPredictedState(), *it->recHit());
    LogTrace("TSGForOIFromL1TkMu") << "TSGForOIFromL1TkMu::makeSeedsFromHits: TSOS for TM " << found << std::endl;
    if (not updatedTSOS.isValid())
      continue;

    edm::OwnVector<TrackingRecHit> seedHits;
    seedHits.push_back(*it->recHit()->hit());
    PTrajectoryStateOnDet const& pstate =
        trajectoryStateTransform::persistentState(updatedTSOS, it->recHit()->geographicalId().rawId());
    LogTrace("TSGForOIFromL1TkMu") << "TSGForOIFromL1TkMu::makeSeedsFromHits: Number of seedHits: " << seedHits.size() << std::endl;
    TrajectorySeed seed(pstate, std::move(seedHits), oppositeToMomentum);
    out.push_back(seed);
    found++;
    numSeedsMade++;
    hitSeedsMade++;
    if (found == numOfHitsToTry_)
      break;
    if (hitSeedsMade > maxHitSeeds_)
      return;
  }

  if (found)
    layerCount++;
}

//
// Calculate the dynamic error SF by analysing the L2
//

double TSGForOIFromL1TkMu::calculateSFFromL1TkMu(const TTTrack<Ref_Phase2TrackerDigi_>& track) const {
  double theSF = 1.0;
  // L1 TkMu direction vs pT blowup - as was previously done:
  // Split into 4 pT ranges: <pT1_, pT1_<pT2_, pT2_<pT3_, <pT4_: 13,30,70
  // Split into different eta ranges depending in pT
  auto p3 = track.momentum();
  auto pt = p3.perp();
  double abseta = std::abs(p3.eta());
  if (pt <= pT1_)
    theSF = SF1_;
  else if (pt > pT1_ && pt <= pT2_) {
    if (abseta <= eta3_)
      theSF = SF3_;
    else if (abseta > eta3_ && abseta <= eta6_)
      theSF = SF2_;
    else if (abseta > eta6_)
      theSF = SF3_;
  } else if (pt > pT2_ && pt <= pT3_) {
    if (abseta <= eta1_)
      theSF = SF6_;
    else if (abseta > eta1_ && abseta <= eta2_)
      theSF = SF4_;
    else if (abseta > eta2_ && abseta <= eta3_)
      theSF = SF6_;
    else if (abseta > eta3_ && abseta <= eta4_)
      theSF = SF1_;
    else if (abseta > eta4_ && abseta <= eta5_)
      theSF = SF1_;
    else if (abseta > eta5_)
      theSF = SF5_;
  } else if (pt > pT3_) {
    if (abseta <= eta3_)
      theSF = SF5_;
    else if (abseta > eta3_ && abseta <= eta4_)
      theSF = SF4_;
    else if (abseta > eta4_ && abseta <= eta5_)
      theSF = SF4_;
    else if (abseta > eta5_)
      theSF = SF5_;
  }

  LogTrace(theCategory_) << "TSGForOIFromL1TkMu::calculateSFFromL2: SF has been calculated as: " << theSF;

  return theSF;
}

//
// calculate Chi^2 of two trajectory states
//
double TSGForOIFromL1TkMu::match_Chi2(const TrajectoryStateOnSurface& tsos1, const TrajectoryStateOnSurface& tsos2) const {
  if (!tsos1.isValid() || !tsos2.isValid())
    return -1.;

  AlgebraicVector5 v(tsos1.localParameters().vector() - tsos2.localParameters().vector());
  AlgebraicSymMatrix55 m(tsos1.localError().matrix() + tsos2.localError().matrix());

  bool ierr = !m.Invert();

  if (ierr) {
    edm::LogInfo("TSGForOIFromL1TkMu") << "Error inverting covariance matrix";
    return -1;
  }

  double est = ROOT::Math::Similarity(v, m);

  return est;
}

//
FreeTrajectoryState TSGForOIFromL1TkMu::initialFreeStateL1TTrack(
    const TTTrack<Ref_Phase2TrackerDigi_>& tk,
    const MagneticField* field,
    bool withErr) const {
  Basic3DVector<float> pos(tk.POCA());
  GlobalPoint gpos(pos);
  GlobalVector gmom = tk.momentum();
  int charge = tk.rInv() > 0.f ? 1 : -1;
  
  auto p3 = tk.momentum();
  float tk_pt = p3.perp();
  //float tk_p = p3.mag();
  float tk_eta = p3.eta();
  //float tk_aeta = std::abs(tk_eta);
  //float tk_phi = p3.phi();
  //float tk_q = tk->rInv() > 0 ? 1. : -1.;
  //float tk_phi = p3.phi();
  //float tk_q = tk.rInv()>0? 1.: -1.;
  //float tk_z  = tk.POCA().z();

  //bool barrel = tk_aeta < 1.1 ? true : false;
  
  bool barrel = true;
  if(fabs(tk_eta)>1.1 ) barrel = false;
  
  GlobalTrajectoryParameters par(gpos, gmom, charge, field);
  if (!withErr)
    return FreeTrajectoryState(par);
  //AlgebraicSymMatrix55 mat = AlgebraicMatrixID();
  AlgebraicSymMatrix55 mat;
  // mat *= 1e-8;
  
  mat[0][0] = (0.25 / tk_pt) * (0.25 / tk_pt);  // sigma^2(charge/abs_momentum)
  if (!barrel)
    mat[0][0] = (0.4 / tk_pt) * (0.4 / tk_pt);
  
  mat[1][1] = 0.05 * 0.05;  // sigma^2(lambda)
  mat[2][2] = 0.2 * 0.2;    // sigma^2(phi)
  //mat[3][3] = 2. * 2.;    // sigma^2(x_transverse))
  //mat[4][4] = 2. * 2.;    // sigma^2(y_transverse))
  //Phase 2 values
  if(tk_pt>20){
    if(fabs(tk_eta)<0.6){
      mat[3][3] = 0.006367*0.006367;    // sigma^2(x_transverse))
      mat[4][4] = 0.08525*0.08525;    // sigma^2(y_transverse))
    }
    else if(fabs(tk_eta)<1.6){
      mat[3][3] = 0.0119825*0.0119825;    // sigma^2(x_transverse))
      mat[4][4] = 0.07075*0.07075;    // sigma^2(y_transverse))
    }
    else if(fabs(tk_eta)<2.4){
      mat[3][3] = 0.017295*0.017295;    // sigma^2(x_transverse))
      mat[4][4] = 0.06295*0.06295;    // sigma^2(y_transverse))
    }
  }
  CurvilinearTrajectoryError error(mat);
  //std::cout<<"Here is the freestate error matrix: "<<error.matrix()<<std::endl;
  return FreeTrajectoryState(par, error);
}

//
//
void TSGForOIFromL1TkMu::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("hltL2Muons", "UpdatedAtVtx"));
  desc.add<double>("minPtOfL1TKMuons", 23.0);
  desc.add<int>("layersToTry", 2);
  desc.add<double>("fixedErrorRescaleFactorForHitless", 2.0);
  desc.add<int>("hitsToTry", 1);
  desc.add<bool>("adjustErrorsDynamicallyForHits", false);
  desc.add<bool>("adjustErrorsDynamicallyForHitless", true);
  desc.add<edm::InputTag>("MeasurementTrackerEvent", edm::InputTag("hltSiStripClusters"));
  desc.add<bool>("UseHitLessSeeds", true);
  desc.add<std::string>("estimator", "hltESPChi2MeasurementEstimator100");
  desc.add<double>("maxEtaForTOB", 1.8);
  desc.add<double>("minEtaForTEC", 0.7);
  desc.addUntracked<bool>("debug", false);
  desc.add<double>("fixedErrorRescaleFactorForHits", 1.0);
  desc.add<unsigned int>("maxSeeds", 20);
  desc.add<unsigned int>("maxHitlessSeeds", 5);
  desc.add<unsigned int>("maxHitSeeds", 1);
  desc.add<double>("pT1", 13.0);
  desc.add<double>("pT2", 30.0);
  desc.add<double>("pT3", 70.0);
  desc.add<double>("eta1", 0.2);
  desc.add<double>("eta2", 0.3);
  desc.add<double>("eta3", 1.0);
  desc.add<double>("eta4", 1.2);
  desc.add<double>("eta5", 1.6);
  desc.add<double>("eta6", 1.4);
  desc.add<double>("eta7", 2.1);
  desc.add<double>("SF1", 3.0);
  desc.add<double>("SF2", 4.0);
  desc.add<double>("SF3", 5.0);
  desc.add<double>("SF4", 7.0);
  desc.add<double>("SF5", 10.0);
  desc.add<double>("SF6", 2.0);
  desc.add<double>("tsosDiff1", 0.2);
  desc.add<double>("tsosDiff2", 0.02);
  desc.add<std::string>("propagatorName", "PropagatorWithMaterialParabolicMf");
  descriptions.add("TSGForOIFromL1TkMu", desc);
}

DEFINE_FWK_MODULE(TSGForOIFromL1TkMu);
