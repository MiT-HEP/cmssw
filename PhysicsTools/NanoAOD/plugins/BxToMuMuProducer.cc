#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/PointingKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include <TLorentzVector.h>
#include <TVector.h>
#include <TMatrix.h>

// 
// BxToMuMuProducer is designed for Bs/d->mumu analysis and it stores the following 
// signatures:
// - Vertexed dimuon pairs
// - BtoKmm 
// 

// TODO
// - add pointing constraint
// - alpha3D
// - 3D impact parameter and its significance
// - 3D flight length and its significance
// - isolation
// - close tracks
// - refitted primary vertex excluding B cand tracks
// - Chi2/ndof

typedef reco::Candidate::LorentzVector LorentzVector;

struct KinematicFitResult{
  bool treeIsValid;
  bool vertexIsValid;
  RefCountedKinematicVertex      refitVertex;
  RefCountedKinematicParticle    refitMother;
  RefCountedKinematicTree        refitTree;
  std::vector<RefCountedKinematicParticle> refitDaughters;
  float lxy, lxyErr, sigLxy, cosAlpha;
  KinematicFitResult():treeIsValid(false),vertexIsValid(false),
		       lxy(-1.0), sigLxy(-1.0), cosAlpha(-999.)
  {}

  bool valid() const {
    return treeIsValid and vertexIsValid;
  }

  void postprocess(const reco::BeamSpot& beamSpot)
  {
    if ( not valid() ) return;
    // displacement information
    TVector v(2);
    v[0] = refitVertex->position().x()-beamSpot.position().x();
    v[1] = refitVertex->position().y()-beamSpot.position().y();

    TMatrix errVtx(2,2);
    errVtx(0,0) = refitVertex->error().cxx();
    errVtx(0,1) = refitVertex->error().matrix()(0,1);
    errVtx(1,0) = errVtx(0,1);
    errVtx(1,1) = refitVertex->error().cyy();

    TMatrix errBS(2,2);
    errBS(0,0) = beamSpot.covariance()(0,0);
    errBS(0,1) = beamSpot.covariance()(0,1);
    errBS(1,0) = beamSpot.covariance()(1,0);
    errBS(1,1) = beamSpot.covariance()(1,1);
    
    lxy = sqrt(v.Norm2Sqr());
    lxyErr = sqrt( v*(errVtx*v) + v*(errBS*v) ) / lxy;
    if (lxyErr > 0) sigLxy = lxy/lxyErr;
    
    // compute cosAlpha
    v[0] = refitVertex->position().x()-beamSpot.position().x();
    v[1] = refitVertex->position().y()-beamSpot.position().y();
    TVector w(2);
    w[0] = refitMother->currentState().globalMomentum().x();
    w[1] = refitMother->currentState().globalMomentum().y();
    cosAlpha = v*w/sqrt(v.Norm2Sqr()*w.Norm2Sqr());
  }
  
  float mass() const
  {
    if ( not valid() ) return -1.0;
    return refitMother->currentState().mass();
  }

  float massErr() const
  {
    if ( not valid() ) return -1.0;
    return sqrt(refitMother->currentState().kinematicParametersError().matrix()(6,6));
  }

  float vtxProb() const
  {
    if ( not valid() ) return -1.0;
    return TMath::Prob((double)refitVertex->chiSquared(), int(rint(refitVertex->degreesOfFreedom())));
  }
  
};

struct KalmanVertexFitResult{
  float vtxProb;
  bool  valid;
  std::vector<LorentzVector> refitVectors;
  GlobalPoint position;
  GlobalError err;
  float lxy, lxyErr, sigLxy;

  KalmanVertexFitResult():vtxProb(-1.0),valid(false),lxy(-1.0),lxyErr(-1.0),sigLxy(-1.0){}

  float mass() const
  {
    if (not valid) return -1.0;
    LorentzVector p4;
    for (auto v: refitVectors)
      p4 += v;
    return p4.mass();
  }
  
  void postprocess(const reco::BeamSpot& bs)
  {
    if (not valid) return;
    // position of the beam spot at a given z value (it takes into account the dxdz and dydz slopes)
    reco::BeamSpot::Point bs_at_z(bs.position(position.z()));
    GlobalPoint xy_displacement(position.x() - bs_at_z.x(),
				position.y() - bs_at_z.y(),
				0);
    lxy = xy_displacement.perp();
    lxyErr = sqrt(err.rerr(xy_displacement));
    if (lxyErr > 0) sigLxy = lxy/lxyErr;
  }
};

LorentzVector makeLorentzVectorFromPxPyPzM(double px, double py, double pz, double m){
  double p2 = px*px+py*py+pz*pz;
  return LorentzVector(px,py,pz,sqrt(p2+m*m));
}

struct GenMatchInfo{
  int mu1_pdgId, mu1_motherPdgId, mu2_pdgId, mu2_motherPdgId, kaon_pdgId, kaon_motherPdgId,
    mm_pdgId, mm_motherPdgId, kmm_pdgId;
  float mu1_pt, mu2_pt, kaon_pt, mm_mass, mm_pt, kmm_mass, kmm_pt;
  GenMatchInfo():mu1_pdgId(0), mu1_motherPdgId(0), mu2_pdgId(0), mu2_motherPdgId(0), 
		 kaon_pdgId(0), kaon_motherPdgId(0), mm_pdgId(0), mm_motherPdgId(0), 
		 kmm_pdgId(0), mu1_pt(0), mu2_pt(0), kaon_pt(0), mm_mass(0), mm_pt(0), 
		 kmm_mass(0), kmm_pt(0){}
};
struct GenEventInfo{};


using namespace std;

class BxToMuMuProducer : public edm::EDProducer {
    
public:
    
  explicit BxToMuMuProducer(const edm::ParameterSet &iConfig);
    
  ~BxToMuMuProducer() override {};
    
    
private:
    
  virtual void produce(edm::Event&, const edm::EventSetup&);

  bool isGoodMuon(const pat::Muon& muon);
    
  KalmanVertexFitResult 
  vertexWithKalmanFitter(edm::ESHandle<TransientTrackBuilder>& theTTBuilder,
			 std::vector<const reco::Track*> trks, 
			 std::vector<float> masses);

  KalmanVertexFitResult 
  vertexMuonsWithKalmanFitter(edm::ESHandle<TransientTrackBuilder>& theTTBuilder,
			      const pat::Muon& muon1,
			      const pat::Muon& muon2);

  KinematicFitResult 
  vertexWithKinematicFitter(edm::ESHandle<TransientTrackBuilder>& theTTBuilder,
			    std::vector<const reco::Track*> trks,
			    std::vector<float> masses);

  KinematicFitResult 
  vertexMuonsWithKinematicFitter(edm::ESHandle<TransientTrackBuilder>& theTTBuilder,
				 const pat::Muon& muon1,
				 const pat::Muon& muon2);

  /// BToKJPsiMuMuFitResult
  KinematicFitResult
  fitBToKJPsiMuMu( RefCountedKinematicParticle jpsi,
		   const pat::PackedCandidate& kaon,
		   edm::ESHandle<TransientTrackBuilder> theTTBuilder,
		   bool applyJpsiMassConstraint);

  KinematicFitResult
  fitBToKJPsiMuMuNew( RefCountedKinematicTree jpsi,
		      const pat::PackedCandidate& kaon,
		      edm::ESHandle<TransientTrackBuilder> theTTBuilder,
		      bool applyJpsiMassConstraint);

  KinematicFitResult
  refitWithPointingConstraint( RefCountedKinematicTree tree,
			       const reco::Vertex& primaryVertex);

  pair<double,double> computeDCA(const pat::PackedCandidate &kaon,
   				 edm::ESHandle<MagneticField> bFieldHandle,
				 reco::BeamSpot beamSpot);
  GenMatchInfo getGenMatchInfo( const edm::View<reco::GenParticle>& genParticles,
				const pat::Muon& muon1,
				const pat::Muon& muon2,
				const pat::PackedCandidate* kaon = 0 );
  // Two track DOCA
  float distanceOfClosestApproach( const reco::Track* track1,
				   const reco::Track* track2,
				   edm::ESHandle<TransientTrackBuilder> theTTBuilder);

  // ----------member data ---------------------------
    
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfCandToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> >   prunedGenToken_;

  bool isMC_;

  double ptMinMu_;
  double etaMaxMu_;
  double ptMinKaon_;

  double etaMaxKaon_;
  double DCASigMinKaon_;
  bool   diMuonCharge_;
  double minBKmmMass_;
  double maxBKmmMass_;
  double maxTwoTrackDOCA_;
    
    float MuonMass_ = 0.10565837;
    float MuonMassErr_ = 3.5*1e-9;
    float KaonMass_ = 0.493677;
    float KaonMassErr_ = 1.6e-5;
    float JPsiMass_ = 3.0969;
    float JPsiMassErr_ = 92.9e-6;

};

BxToMuMuProducer::BxToMuMuProducer(const edm::ParameterSet &iConfig):
beamSpotToken_( consumes<reco::BeamSpot> ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
vertexToken_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) ),
muonToken_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
pfCandToken_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
prunedGenToken_( consumes<edm::View<reco::GenParticle>> ( iConfig.getParameter<edm::InputTag>( "prunedGenParticleCollection" ) ) ),
isMC_(            iConfig.getParameter<bool>( "isMC" ) ),
ptMinMu_(         iConfig.getParameter<double>( "MuonMinPt" ) ),
etaMaxMu_(        iConfig.getParameter<double>( "MuonMaxEta" ) ),
ptMinKaon_(       iConfig.getParameter<double>( "KaonMinPt" ) ),
etaMaxKaon_(      iConfig.getParameter<double>( "KaonMaxEta" ) ),
DCASigMinKaon_(   iConfig.getParameter<double>( "KaonMinDCASig" ) ),
diMuonCharge_(    iConfig.getParameter<bool>(   "DiMuonChargeCheck" ) ),
minBKmmMass_(     iConfig.getParameter<double>( "minBKmmMass" ) ),
maxBKmmMass_(     iConfig.getParameter<double>( "maxBKmmMass" ) ),
maxTwoTrackDOCA_( iConfig.getParameter<double>( "maxTwoTrackDOCA" ) )
{
    produces<pat::CompositeCandidateCollection>("DiMuon");
    produces<pat::CompositeCandidateCollection>("BToKmumu");
}

bool BxToMuMuProducer::isGoodMuon(const pat::Muon& muon){
  if ( not muon.isLooseMuon() ) return false;
  if ( not muon.isTrackerMuon() ) return false;
  if ( not muon.innerTrack()->quality(reco::Track::highPurity) ) return false; 
  if ( muon.pt() < ptMinMu_ || fabs(muon.eta()) > etaMaxMu_ ) return false;
  return true;
}

void addFitInfo(pat::CompositeCandidate& cand, const KinematicFitResult& fit, std::string name){
  cand.addUserInt(   name+"_valid",       fit.valid() );
  cand.addUserFloat( name+"_vtx_prob",    fit.vtxProb() );
  cand.addUserFloat( name+"_mass",        fit.mass() );
  cand.addUserFloat( name+"_massErr",     fit.massErr() );
  cand.addUserFloat( name+"_lxy",         fit.lxy );
  cand.addUserFloat( name+"_sigLxy",      fit.sigLxy );
  cand.addUserFloat( name+"_cosAlpha",    fit.cosAlpha );
}


void BxToMuMuProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    edm::ESHandle<MagneticField> bFieldHandle;
    edm::ESHandle<TransientTrackBuilder> theTTBuilder;

    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    edm::Handle<reco::VertexCollection> vertexHandle;
    
    iEvent.getByToken(beamSpotToken_, beamSpotHandle);
    
    if ( ! beamSpotHandle.isValid() ) {
        edm::LogError("BxToMuMuProducer") << "No beam spot available from EventSetup" ;
    }
    
    const reco::BeamSpot& beamSpot(*beamSpotHandle);
    
    // iEvent.getByToken(vertexToken_, vertexHandle);
    // const reco::Vertex & primaryVertex = vertexHandle->front();

    edm::Handle<std::vector<pat::Muon>> muonHandle;
    edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle;
    edm::Handle<edm::View<pat::PackedCandidate>> lostTrackHandle;
    
    iEvent.getByToken(muonToken_, muonHandle);
    iEvent.getByToken(pfCandToken_, pfCandHandle);
    
    edm::Handle<edm::View<reco::GenParticle> > prunedGenParticles;
    if ( isMC_ ) iEvent.getByToken(prunedGenToken_,prunedGenParticles);
    
    auto nMuons = muonHandle->size();
    auto nPFCands = pfCandHandle->size();
    // unsigned int lostTrackNumber = useLostTracks_ ? lostTrackHandle->size() : 0;
    
    // Output collection
    auto dimuon = std::make_unique<pat::CompositeCandidateCollection>();
    auto btokmm = std::make_unique<pat::CompositeCandidateCollection>();
    AddFourMomenta addP4;

    // Build dimuon candidates first

    if ( nMuons > 1 ){
      // Ensure that muon1.pt > muon2.pt
      for (unsigned int i = 0; i < nMuons; ++i) {
	const pat::Muon & muon1 = muonHandle->at(i);
        if (not isGoodMuon(muon1)) continue;
	for (unsigned int j = 0; j < nMuons; ++j) {
	  if (i==j) continue;
	  const pat::Muon & muon2 = muonHandle->at(j);
	  if (muon2.pt() > muon1.pt()) continue;
	  if (not isGoodMuon(muon2)) continue;
	  
	  if (maxTwoTrackDOCA_>0 and 
	      distanceOfClosestApproach(muon1.innerTrack().get(),
					muon2.innerTrack().get(),
					theTTBuilder) > maxTwoTrackDOCA_) continue;
	  if (diMuonCharge_ && muon1.charge()*muon2.charge()>0) continue;

	  pat::CompositeCandidate dimuonCand;
	  dimuonCand.addDaughter( muon1 , "muon1");
	  dimuonCand.addDaughter( muon2 , "muon2");
	  addP4.set( dimuonCand );

	  dimuonCand.addUserInt("mu1_index", i);
	  dimuonCand.addUserInt("mu2_index", j);

	  // Kalman Vertex Fit
	  auto kalmanMuMuVertexFit = vertexMuonsWithKalmanFitter(theTTBuilder, muon1, muon2);
	  kalmanMuMuVertexFit.postprocess(beamSpot);
	  dimuonCand.addUserInt(   "kalman_valid",    kalmanMuMuVertexFit.valid);
	  dimuonCand.addUserFloat( "kalman_vtx_prob", kalmanMuMuVertexFit.vtxProb);
	  dimuonCand.addUserFloat( "kalman_mass",     kalmanMuMuVertexFit.mass() );
	  dimuonCand.addUserFloat( "kalman_lxy",      kalmanMuMuVertexFit.lxy );
	  dimuonCand.addUserFloat( "kalman_sigLxy",   kalmanMuMuVertexFit.sigLxy );
	  
	  // Kinematic Fits
	  auto kinematicMuMuVertexFit = vertexMuonsWithKinematicFitter(theTTBuilder, muon1, muon2);
	  kinematicMuMuVertexFit.postprocess(beamSpot);
	  addFitInfo(dimuonCand, kinematicMuMuVertexFit, "kin");
	  if (isMC_){
	    auto gen_mm = getGenMatchInfo(*prunedGenParticles.product(),muon1,muon2);
	    // int mu1_pdgId, mu1_motherPdgId, mu2_pdgId, mu2_motherPdgId, kaon_pdgId, kaon_motherPdgId,
	    // mm_pdgId, mm_motherPdgId, kmm_pdgId;
	    // float mu1_pt, mu2_pt, kaon_pt, mm_mass, mm_pt, kmm_mass, kmm_pt;
	    dimuonCand.addUserInt(  "gen_mu1_pdgId",  gen_mm.mu1_pdgId);
	    dimuonCand.addUserInt(  "gen_mu1_mpdgId", gen_mm.mu1_motherPdgId);
	    dimuonCand.addUserFloat("gen_mu1_pt",     gen_mm.mu1_pt);
	    dimuonCand.addUserInt(  "gen_mu2_pdgId",  gen_mm.mu2_pdgId);
	    dimuonCand.addUserInt(  "gen_mu2_mpdgId", gen_mm.mu2_motherPdgId);
	    dimuonCand.addUserFloat("gen_mu2_pt",     gen_mm.mu2_pt);
	    dimuonCand.addUserFloat("gen_mass",       gen_mm.mm_mass);
	    dimuonCand.addUserFloat("gen_pt",         gen_mm.mm_pt);
	    dimuonCand.addUserInt(  "gen_pdgId",      gen_mm.mm_pdgId);
	    dimuonCand.addUserInt(  "gen_mpdgId",     gen_mm.mm_motherPdgId);
	  }

	  auto imm = dimuon->size();
	  // BtoJpsiK
	  if (fabs(kinematicMuMuVertexFit.mass()-3.1)<0.2){
	    for (unsigned int k = 0; k < nPFCands; ++k) {
	      const pat::PackedCandidate & pfCand = (*pfCandHandle)[k];
	      if(abs(pfCand.pdgId())!=211) continue; //Charged hadrons
	      if(!pfCand.hasTrackDetails()) continue;
	      if(pfCand.pt()<ptMinKaon_ || abs(pfCand.eta())>etaMaxKaon_) continue;
	      if(deltaR(muon1, pfCand) < 0.01 || deltaR(muon2, pfCand) < 0.01) continue;
	    
	      if (maxTwoTrackDOCA_>0 and 
		  distanceOfClosestApproach(muon1.innerTrack().get(),
					    pfCand.bestTrack(),
					    theTTBuilder) > maxTwoTrackDOCA_) continue;	      
	      if (maxTwoTrackDOCA_>0 and 
		  distanceOfClosestApproach(muon2.innerTrack().get(),
					    pfCand.bestTrack(),
					    theTTBuilder) > maxTwoTrackDOCA_) continue;

	      double kmm_mass = (muon1.p4()+muon2.p4()+pfCand.p4()).mass();
	      if ( kmm_mass<minBKmmMass_ || kmm_mass>maxBKmmMass_ ) continue;

	      pat::CompositeCandidate btokmmCand;
	      btokmmCand.addUserInt("mm_index", imm);
	      btokmmCand.addUserFloat("kaon_pt", pfCand.pt());
	      btokmmCand.addUserFloat("kaon_eta", pfCand.eta());
	      btokmmCand.addUserFloat("kaon_phi", pfCand.phi());
	      btokmmCand.addUserInt("kaon_charge", pfCand.charge());

	      if (isMC_){
		auto gen_kmm = getGenMatchInfo(*prunedGenParticles.product(),muon1,muon2,&pfCand);
		btokmmCand.addUserInt(  "gen_kaon_pdgId",  gen_kmm.kaon_pdgId);
		btokmmCand.addUserInt(  "gen_kaon_mpdgId", gen_kmm.kaon_motherPdgId);
		btokmmCand.addUserFloat("gen_kaon_pt",     gen_kmm.kaon_pt);
		btokmmCand.addUserFloat("gen_mass",        gen_kmm.kmm_mass);
		btokmmCand.addUserFloat("gen_pt",          gen_kmm.kmm_pt);
		btokmmCand.addUserInt(  "gen_pdgId",       gen_kmm.kmm_pdgId);
	      }
	      // if (pfCand.genParticle()){
	      // 	btokmmCand.addUserInt("kaon_mc_pdgId", pfCand.genParticle().pdgId());
	      // } else {
	      // 	btokmmCand.addUserInt("kaon_mc_pdgId", 0);
	      // }
	      
	      auto bToKJPsiMuMuNoMassConstraint = fitBToKJPsiMuMu(kinematicMuMuVertexFit.refitMother, pfCand, theTTBuilder, false);
	      bToKJPsiMuMuNoMassConstraint.postprocess(beamSpot);
	      addFitInfo(btokmmCand, bToKJPsiMuMuNoMassConstraint, "nomc");
	      
	      // worse performing option
	      // auto bToKJPsiMuMuWithMassConstraint = fitBToKJPsiMuMu(kinematicMuMuVertexFit.refitMother, pfCand, theTTBuilder, true);
	      // bToKJPsiMuMuWithMassConstraint.postprocess(beamSpot);
	      // addFitInfo(btokmmCand, bToKJPsiMuMuWithMassConstraint, "jpsimc");

	      auto bToKJPsiMuMu_MC = fitBToKJPsiMuMuNew(kinematicMuMuVertexFit.refitTree, pfCand, theTTBuilder, true);
	      bToKJPsiMuMu_MC.postprocess(beamSpot);
	      addFitInfo(btokmmCand, bToKJPsiMuMu_MC, "jpsimc");
	      
	      // broken pointing constraint
	      // auto bToKJPsiMuMu_MC_PC = refitWithPointingConstraint(bToKJPsiMuMu_MC.refitTree, primaryVertex);
	      // bToKJPsiMuMu_MC_PC.postprocess(beamSpot);
	      // addFitInfo(btokmmCand, bToKJPsiMuMu_MC_PC, "mcpc");

	      btokmm->push_back(btokmmCand);
	    }                    
	  }
        
	  dimuon->push_back(dimuonCand);
	}
      }
    }
    
    iEvent.put(std::move(dimuon),"DiMuon");
    iEvent.put(std::move(btokmm),"BToKmumu");
    
}

KalmanVertexFitResult 
BxToMuMuProducer::vertexWithKalmanFitter(edm::ESHandle<TransientTrackBuilder>& theTTBuilder,
					 std::vector<const reco::Track*> trks, 
					 std::vector<float> masses){
  if (trks.size()!=masses.size()) 
    throw cms::Exception("Error") << "number of tracks and number of masses should match";
  KalmanVertexFitResult results;
  std::vector<reco::TransientTrack> transTrks;
  for (auto trk: trks){
    transTrks.push_back((*theTTBuilder).build(trk));
  }
  KalmanVertexFitter kvf(true);
  TransientVertex tv = kvf.vertex(transTrks);

  if ( tv.isValid() ){
    results.vtxProb = TMath::Prob(tv.totalChiSquared(), (int)tv.degreesOfFreedom());
    results.valid = true;
    results.position = tv.position();
    results.err = tv.positionError();
    if (tv.hasRefittedTracks()){
      assert(tv.refittedTracks().size()==transTrks.size());
      for (unsigned int i=0; i<transTrks.size(); ++i){
	// Is it safe to assume that the order hasn't changed?
	GlobalVector gvP = tv.refittedTracks()[i].trajectoryStateClosestToPoint(tv.position()).momentum();
	// GlobalVector gvP = tv.originalTracks()[i].trajectoryStateClosestToPoint(tv.position()).momentum();
	results.refitVectors.push_back(makeLorentzVectorFromPxPyPzM(gvP.x(), gvP.y(), gvP.z(), masses[i]));
      }
    }
  }
  return results;
}

KalmanVertexFitResult 
BxToMuMuProducer::vertexMuonsWithKalmanFitter(edm::ESHandle<TransientTrackBuilder>& theTTBuilder,
					      const pat::Muon& muon1,
					      const pat::Muon& muon2)
{
  std::vector<const reco::Track*> trks;
  std::vector<float> masses;
  trks.push_back( muon1.innerTrack().get() );
  masses.push_back(MuonMass_);
  trks.push_back( muon2.innerTrack().get() );
  masses.push_back(MuonMass_);
  return vertexWithKalmanFitter(theTTBuilder,trks,masses);
}


KinematicFitResult 
BxToMuMuProducer::vertexWithKinematicFitter(edm::ESHandle<TransientTrackBuilder>& theTTBuilder,
					    std::vector<const reco::Track*> trks,
					    std::vector<float> masses)
{
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideKinematicVertexFit
  if ( trks.size() != masses.size() ) 
    throw cms::Exception("Error") << "number of tracks and number of masses should match";

  std::vector<reco::TransientTrack> transTrks;

  KinematicParticleFactoryFromTransientTrack factory;
  KinematicParticleVertexFitter fitter;
    
  std::vector<RefCountedKinematicParticle> particles;

  double chi = 0.;
  double ndf = 0.;

  for (auto trk: trks){
    transTrks.push_back((*theTTBuilder).build(trk));
    particles.push_back(factory.particle(transTrks.back(),MuonMass_,chi,ndf,MuonMassErr_));
  }

  RefCountedKinematicTree vertexFitTree = fitter.fit(particles);
  KinematicFitResult result;
    
  if ( !vertexFitTree->isValid()) return result;
  
  result.treeIsValid = true;

  vertexFitTree->movePointerToTheTop();
  result.refitVertex = vertexFitTree->currentDecayVertex();
  result.refitMother = vertexFitTree->currentParticle();
  result.refitTree   = vertexFitTree;
  if ( !result.refitVertex->vertexIsValid()) return result;

  result.vertexIsValid = true;

  // extract the re-fitted tracks
  vertexFitTree->movePointerToTheTop();
  
  if ( vertexFitTree->movePointerToTheFirstChild() ){
    do {
      result.refitDaughters.push_back(vertexFitTree->currentParticle());
    } while (vertexFitTree->movePointerToTheNextChild());
  }
  return result;
}

KinematicFitResult
BxToMuMuProducer::vertexMuonsWithKinematicFitter(edm::ESHandle<TransientTrackBuilder>& theTTBuilder,
						 const pat::Muon& muon1,
						 const pat::Muon& muon2)
{
  std::vector<const reco::Track*> trks;
  std::vector<float> masses;
  trks.push_back( muon1.innerTrack().get() );
  masses.push_back(MuonMass_);
  trks.push_back( muon2.innerTrack().get() );
  masses.push_back(MuonMass_);
  return vertexWithKinematicFitter(theTTBuilder,trks,masses);
}

KinematicFitResult
BxToMuMuProducer::fitBToKJPsiMuMu( RefCountedKinematicParticle refitMuMu,
				   const pat::PackedCandidate &kaon,
				   edm::ESHandle<TransientTrackBuilder> theTTBuilder,
				   bool applyJpsiMassConstraint)
{
  const reco::TransientTrack mmTT = refitMuMu->refittedTransientTrack();
  const reco::TransientTrack kaonTT = theTTBuilder->build(kaon.bestTrack());

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter fitter;

  std::vector<RefCountedKinematicParticle> BToKMuMuParticles;
  double chi = 0.;
  double ndf = 0.;

  float MuMu_mass = refitMuMu->currentState().mass();
  float MuMu_mass_err = sqrt(refitMuMu->currentState().kinematicParametersError().matrix()(6,6));

  if ( applyJpsiMassConstraint ){
    MuMu_mass = JPsiMass_;
    MuMu_mass_err = JPsiMassErr_;
  }

  BToKMuMuParticles.push_back(partFactory.particle(mmTT,MuMu_mass,chi,ndf,MuMu_mass_err));
  BToKMuMuParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,KaonMassErr_));

  RefCountedKinematicTree vertexFitTree = fitter.fit(BToKMuMuParticles);
  KinematicFitResult result; 

  if ( !vertexFitTree->isValid()) return result;

  result.treeIsValid = true;

  vertexFitTree->movePointerToTheTop();
  result.refitVertex = vertexFitTree->currentDecayVertex();
  result.refitMother = vertexFitTree->currentParticle();
  result.refitTree   = vertexFitTree;

  if ( !result.refitVertex->vertexIsValid()) return result;

  result.vertexIsValid = true;

  // extract the re-fitted tracks
  vertexFitTree->movePointerToTheTop();

  if ( vertexFitTree->movePointerToTheFirstChild() ){
    do {
      result.refitDaughters.push_back(vertexFitTree->currentParticle());
    } while (vertexFitTree->movePointerToTheNextChild());
  }
  return result;
}

KinematicFitResult
BxToMuMuProducer::fitBToKJPsiMuMuNew( RefCountedKinematicTree jpsiTree,
				      const pat::PackedCandidate& kaon,
				      edm::ESHandle<TransientTrackBuilder> theTTBuilder,
				      bool applyJpsiMassConstraint)
{
  KinematicFitResult result; 
  if ( !jpsiTree->isValid()) return result;

  KinematicConstraint* jpsi_mc(0);
  if (applyJpsiMassConstraint){
    ParticleMass jpsi = JPsiMass_;
    // jpsi mass constraint fit
    KinematicParticleFitter csFitter;
    float jp_m_sigma = JPsiMassErr_;
    // FIXME: memory leak
    jpsi_mc = new MassKinematicConstraint(jpsi, jp_m_sigma);
    jpsiTree = csFitter.fit(jpsi_mc, jpsiTree);
  }

  const reco::TransientTrack kaonTT = theTTBuilder->build(kaon.bestTrack());

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter fitter;

  std::vector<RefCountedKinematicParticle> BToKMuMuParticles;
  double chi = 0.;
  double ndf = 0.;

  jpsiTree->movePointerToTheTop();
  BToKMuMuParticles.push_back(jpsiTree->currentParticle());
  BToKMuMuParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,KaonMassErr_));

  RefCountedKinematicTree vertexFitTree = fitter.fit(BToKMuMuParticles);

  if ( !vertexFitTree->isValid()) return result;

  result.treeIsValid = true;

  vertexFitTree->movePointerToTheTop();
  result.refitVertex = vertexFitTree->currentDecayVertex();
  result.refitMother = vertexFitTree->currentParticle();
  result.refitTree   = vertexFitTree;

  if ( !result.refitVertex->vertexIsValid()) return result;

  result.vertexIsValid = true;

  // extract the re-fitted tracks
  vertexFitTree->movePointerToTheTop();

  if ( vertexFitTree->movePointerToTheFirstChild() ){
    do {
      result.refitDaughters.push_back(vertexFitTree->currentParticle());
    } while (vertexFitTree->movePointerToTheNextChild());
  }
  return result;
}

KinematicFitResult
BxToMuMuProducer::refitWithPointingConstraint( RefCountedKinematicTree tree,
					       const reco::Vertex& primaryVertex)
{
  KinematicFitResult result; 
  if ( !tree->isValid()) return result;

  GlobalPoint pv(primaryVertex.position().x(), 
		 primaryVertex.position().y(), 
		 primaryVertex.position().z());
  // FIXIT: potential memory leak
  PointingKinematicConstraint* pointing_constraint = new PointingKinematicConstraint(pv);

  KinematicParticleFitter fitter;

  tree->movePointerToTheTop(); // not sure it's needed

  RefCountedKinematicTree refittedTree = fitter.fit(pointing_constraint,tree);

  if ( !refittedTree->isValid()) return result;

  result.treeIsValid = true;

  refittedTree->movePointerToTheTop();
  result.refitVertex = refittedTree->currentDecayVertex();
  result.refitMother = refittedTree->currentParticle();
  result.refitTree   = refittedTree;

  if ( !result.refitVertex->vertexIsValid()) return result;

  result.vertexIsValid = true;

  // extract the re-fitted tracks
  refittedTree->movePointerToTheTop();

  if ( refittedTree->movePointerToTheFirstChild() ){
    do {
      result.refitDaughters.push_back(refittedTree->currentParticle());
    } while (refittedTree->movePointerToTheNextChild());
  }
  return result;
}

pair<double,double> BxToMuMuProducer::computeDCA(const pat::PackedCandidate &kaon,
                                                 edm::ESHandle<MagneticField> bFieldHandle,
                                                 reco::BeamSpot beamSpot){

  const reco::TransientTrack trackTT((*(kaon.bestTrack())), &(*bFieldHandle));

  TrajectoryStateClosestToPoint theDCAXBS = trackTT.trajectoryStateClosestToPoint( GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()) );  
  
  double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
  double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
  pair<double,double> DCA = make_pair(DCABS,DCABSErr);
    
  return DCA;
}

bool dr_match(const LorentzVector& reco , const LorentzVector& gen){
  if (fabs(reco.pt()-gen.pt())/gen.pt()<0.1 and deltaR(reco,gen)<0.02)
    return true;
  return false;
}

GenMatchInfo BxToMuMuProducer::getGenMatchInfo( const edm::View<reco::GenParticle>& genParticles,
						const pat::Muon& muon1,
						const pat::Muon& muon2,
						const pat::PackedCandidate* kaon )
{
  auto result = GenMatchInfo();
  const reco::GenParticle* mc_mu1(0);
  const reco::GenParticle* mc_mu2(0);
  const reco::GenParticle* mc_kaon(0);
  const reco::Candidate*   mm(0);
  if (muon1.genParticle()){
    // should we use sim instead for more information?
    mc_mu1 = muon1.genParticle();
    result.mu1_pdgId = muon1.genParticle()->pdgId();
    result.mu1_pt    = muon1.genParticle()->pt();
    if (muon1.genParticle()->mother()){
      result.mu1_motherPdgId = muon1.genParticle()->mother()->pdgId();
    }
  }
  if (muon2.genParticle()){
    mc_mu2 = muon2.genParticle();
    result.mu2_pdgId = muon2.genParticle()->pdgId();
    result.mu2_pt    = muon2.genParticle()->pt();
    if (muon2.genParticle()->mother()){
      result.mu2_motherPdgId = muon2.genParticle()->mother()->pdgId();
    }
  }
  if ( mc_mu1 and mc_mu2 ){
    mm = mc_mu1->mother();
    if ( mm and mm == mc_mu2->mother() ){
      result.mm_mass  = mm->mass();
      result.mm_pt    = mm->pt();
      result.mm_pdgId = mm->pdgId();
      if (mm->mother())
	result.mm_motherPdgId = mm->mother()->pdgId();
    }
  }
  if (kaon){
    for (auto const & genParticle: genParticles){
      if (dr_match(kaon->p4(),genParticle.p4())){
	mc_kaon = &genParticle;
	  result.kaon_pdgId = genParticle.pdgId();
	  result.kaon_pt    = genParticle.pt();
	  if (genParticle.mother()){
	    result.kaon_motherPdgId = genParticle.mother()->pdgId();
	  }
	  break;
	}
    }
    if (mm and mc_kaon and mc_kaon->mother()){
      if (mm == mc_kaon->mother() or mm->mother() == mc_kaon->mother()){
	result.kmm_pdgId = mc_kaon->mother()->pdgId();
	result.kmm_mass  = mc_kaon->mother()->mass();
	result.kmm_pt = mc_kaon->mother()->pt();
      }
    }
  }

  return result;
}

float BxToMuMuProducer::distanceOfClosestApproach( const reco::Track* track1,
						   const reco::Track* track2,
						   edm::ESHandle<TransientTrackBuilder> theTTBuilder)
{
  TwoTrackMinimumDistance md;
  const reco::TransientTrack tt1 = theTTBuilder->build(track1);
  const reco::TransientTrack tt2 = theTTBuilder->build(track2);
  if ( not md.calculate( tt1.initialFreeState(), tt2.initialFreeState() ) ) return -1.0;
  return md.distance();
}

// GenInfoSummary HeavyFlavDileptonNtupleMakerMiniAOD::getGenInfoSummary(const edm::Event& iEvent){
//   GenInfoSummary summary;
  
//   edm::Handle<edm::View<reco::GenParticle> > pruned;
//   iEvent.getByToken(pruned_gen_token,pruned);

//   // Packed particles are all the status 1, so usable to remake jets
//   // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
//   edm::Handle<edm::View<pat::PackedGenParticle> > packed;
//   iEvent.getByToken(packed_gen_token,packed);

//   //let's try to find all status1 originating directly from a B meson decay
//   for (auto const & bMeson: *pruned){
//     if(abs(bMeson.pdgId()) != 521 and abs(bMeson.pdgId()) != 511) continue;
//     std::vector<const reco::Candidate*> b_muons;
//     std::vector<const reco::Candidate*> b_hadrons;
//     std::vector<const reco::Candidate*> other_muons;
//     for (auto const& pa: *packed){ 
//       //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
//       const reco::Candidate * mother = pa.mother(0) ;
//       // FIXME: may need to check particle status to avoid intermediate gen particles
//       if(mother != nullptr && isAncestor( &bMeson , mother)){
// 	if (abs(pa.pdgId())==13)
// 	  b_muons.push_back(&pa);
// 	if (abs(pa.pdgId())==211 or abs(pa.pdgId())==321)
// 	  b_hadrons.push_back(&pa);
//       }else{
// 	if (abs(pa.pdgId())==13)
// 	  other_muons.push_back(&pa);
//       }
//     }
//     if (b_muons.size()==2 && b_hadrons.size()>0){
//       summary.b = &bMeson;
//       summary.b_muons = b_muons;
//       summary.b_hadrons = b_hadrons;
//       std::sort(other_muons.begin(),other_muons.end(), [](const reco::Candidate* a,const reco::Candidate* b){
// 	  return a->pt() > b->pt();
// 	});
//       summary.other_muons = other_muons;
//     }
//   }

//   return summary;
// }


DEFINE_FWK_MODULE(BxToMuMuProducer);
