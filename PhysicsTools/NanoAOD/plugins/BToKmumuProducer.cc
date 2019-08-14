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
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include <TLorentzVector.h>
#include <TVector.h>
#include <TMatrix.h>


//
// class declaration
//

using namespace std;

class BToKmumuProducer : public edm::EDProducer {
    
public:
    
    explicit BToKmumuProducer(const edm::ParameterSet &iConfig);
    
    ~BToKmumuProducer() override {};
    
    
private:
    
    virtual void produce(edm::Event&, const edm::EventSetup&);
    
    bool MuMuVertexRefitting(const pat::Muon & muon1,
                             const pat::Muon & muon2,
                             edm::ESHandle<TransientTrackBuilder> theTTBuilder,
                             RefCountedKinematicVertex &refitVertex,
                             RefCountedKinematicParticle &refitMuMu,
                             RefCountedKinematicParticle &refitMu1,
                             RefCountedKinematicParticle &refitMu2);
    
    bool BToKMuMuVertexRefitting(const pat::Muon &muon1,
                                 const pat::Muon &muon2,
                                 const pat::PackedCandidate &kaon,
                                 edm::ESHandle<TransientTrackBuilder> theTTBuilder,
                                 RefCountedKinematicVertex &refitVertex,
                                 RefCountedKinematicParticle &refitBToKMuMu,
                                 RefCountedKinematicParticle &refitMu1,
                                 RefCountedKinematicParticle &refitMu2,
                                 RefCountedKinematicParticle &refitKaon);

    bool BToKJPsiMuMuVertexRefitting(const RefCountedKinematicParticle refitMuMu,
				     const pat::PackedCandidate &kaon,
				     edm::ESHandle<TransientTrackBuilder> theTTBuilder,
				     RefCountedKinematicVertex &refitVertex,
				     RefCountedKinematicParticle &refitBToKJPsiMuMu,
				     RefCountedKinematicParticle &refitJPsi,
				     RefCountedKinematicParticle &refitKaon);
    
    pair<double,double> computeLS(RefCountedKinematicVertex refitVertex,
				  reco::BeamSpot beamSpot);
    
    double computeCosAlpha(RefCountedKinematicParticle refitBToKMuMu,
                           RefCountedKinematicVertex vertexFitTree,
                           reco::BeamSpot beamSpot);

    pair<double,double> computeDCA(const pat::PackedCandidate &kaon,
				   edm::ESHandle<MagneticField> bFieldHandle,
				   reco::BeamSpot beamSpot);
    
    // ----------member data ---------------------------
    
    edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
    edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
    edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate>> PFCandSrc_;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate>> lostTrackSrc_;
    
    double ptMinMu_;
    double etaMaxMu_;
    double ptMinKaon_;
    double etaMaxKaon_;
    double DCASigMinKaon_;
    bool diMuonCharge_;
    double JPsiMassConstraint_;
    bool save2TrkRefit_;
    bool useLostTracks_;
    
    float MuonMass_ = 0.10565837;
    float MuonMassErr_ = 3.5*1e-9;
    float KaonMass_ = 0.493677;
    float KaonMassErr_ = 1.6e-5;
    //float JPsiMass_ = 3.096916;  //Configurable parameter
    float JPsiMassErr_ = 0.011;

};



BToKmumuProducer::BToKmumuProducer(const edm::ParameterSet &iConfig):
beamSpotSrc_( consumes<reco::BeamSpot> ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
vertexSrc_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) ),
muonSrc_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
PFCandSrc_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
lostTrackSrc_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "lostTrackCollection" ) ) ),
ptMinMu_( iConfig.getParameter<double>( "MuonMinPt" ) ),
etaMaxMu_( iConfig.getParameter<double>( "MuonMaxEta" ) ),
ptMinKaon_( iConfig.getParameter<double>( "KaonMinPt" ) ),
etaMaxKaon_( iConfig.getParameter<double>( "KaonMaxEta" ) ),
DCASigMinKaon_( iConfig.getParameter<double>( "KaonMinDCASig" ) ),
diMuonCharge_( iConfig.getParameter<bool>( "DiMuonChargeCheck" ) ),
JPsiMassConstraint_( iConfig.getParameter<double>( "JPsiMassConstraint" ) ),
save2TrkRefit_( iConfig.getParameter<bool>( "save2TrackRefit" ) ),
useLostTracks_( iConfig.getParameter<bool>( "useLostTracks" ) )
{
    produces<pat::CompositeCandidateCollection>();
}


void BToKmumuProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    edm::ESHandle<MagneticField> bFieldHandle;
    edm::ESHandle<TransientTrackBuilder> theTTBuilder;

    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    edm::Handle<reco::VertexCollection> vertexHandle;
    
    iEvent.getByToken(beamSpotSrc_, beamSpotHandle);
    
    if ( ! beamSpotHandle.isValid() ) {
        edm::LogError("BToKmumuProducer") << "No beam spot available from EventSetup" ;
    }
    
    reco::BeamSpot beamSpot = *beamSpotHandle;
    
    iEvent.getByToken(vertexSrc_, vertexHandle);
    const reco::Vertex & PV = vertexHandle->front();

    edm::Handle<std::vector<pat::Muon>> muonHandle;
    edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle;
    edm::Handle<edm::View<pat::PackedCandidate>> lostTrackHandle;
    
    iEvent.getByToken(muonSrc_, muonHandle);
    iEvent.getByToken(PFCandSrc_, pfCandHandle);
    if(useLostTracks_) iEvent.getByToken(lostTrackSrc_, lostTrackHandle);
    
    unsigned int muonNumber = muonHandle->size();
    unsigned int pfCandNumber = pfCandHandle->size();
    unsigned int lostTrackNumber = useLostTracks_ ? lostTrackHandle->size() : 0;
    
    // Output collection
    std::unique_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );
    
    
    if(muonNumber>1){
        
        // loop on all the mumuK triplets
        for (unsigned int i = 0; i < muonNumber; ++i) {

            const pat::Muon & muon1 = (*muonHandle)[i];
            
            if(!(muon1.isLooseMuon() && muon1.isSoftMuon(PV))) continue;
            if(muon1.pt()<ptMinMu_ || abs(muon1.eta())>etaMaxMu_) continue;
            
            for (unsigned int j = 0; j < muonNumber; ++j) {
                
	        if(i==j) continue;

                const pat::Muon & muon2 = (*muonHandle)[j];

                if(muon1.pt()<muon2.pt()) continue; //Muon 1 is always saved as the leading one
                if(!(muon2.isLooseMuon() && muon2.isSoftMuon(PV))) continue;
                if(muon2.pt()<ptMinMu_ || abs(muon2.eta())>etaMaxMu_) continue;
                
                if(diMuonCharge_ && muon1.charge()*muon2.charge()>0) continue;

                bool passedDiMuon = false;

                double MuMuLSBS = -1.;
                double MuMuLSBSErr = -1.;
                double MuMuVtx_Chi2 = -1.;
                double MuMuVtx_CL = -1.;
                double MuMu_mass_err = -1.;

                RefCountedKinematicParticle refitMuMu;
                math::XYZVector refitMuMuV3D;

                if(save2TrkRefit_){

                  RefCountedKinematicVertex refitVertexMuMu;
                  RefCountedKinematicParticle refitMu1_MuMu;
                  RefCountedKinematicParticle refitMu2_MuMu;
                  
                  passedDiMuon = MuMuVertexRefitting(muon1, muon2,
                                                     theTTBuilder,
                                                     refitVertexMuMu,
                                                     refitMuMu,
                                                     refitMu1_MuMu,
                                                     refitMu2_MuMu);
                  
                  if (passedDiMuon){
                    
                     math::XYZVector refitMu1V3D_MuMu = refitMu1_MuMu->refittedTransientTrack().track().momentum();
                     math::XYZVector refitMu2V3D_MuMu = refitMu2_MuMu->refittedTransientTrack().track().momentum();
                     refitMuMuV3D = refitMu1V3D_MuMu + refitMu2V3D_MuMu;
                                        
                     MuMu_mass_err = sqrt(refitMuMu->currentState().kinematicParametersError().matrix()(6,6));

                     pair<double,double> MuMuLS = computeLS(refitVertexMuMu,beamSpot);
                     MuMuLSBS = MuMuLS.first;
                     MuMuLSBSErr = MuMuLS.second;

                     MuMuVtx_Chi2 = (double)refitVertexMuMu->chiSquared();
                     MuMuVtx_CL = TMath::Prob((double)refitVertexMuMu->chiSquared(),
					      int(rint(refitVertexMuMu->degreesOfFreedom())));

                  }

                }

                //Kaon
                for (unsigned int k = 0; k < (pfCandNumber+lostTrackNumber); ++k) {
                    bool isPFCand = k<pfCandNumber;
                    const pat::PackedCandidate & pfCand = isPFCand ? (*pfCandHandle)[k] : (*lostTrackHandle)[k-pfCandNumber];
                    if(abs(pfCand.pdgId())!=211) continue; //Charged hadrons
                    if(!pfCand.hasTrackDetails()) continue;
                    if(pfCand.pt()<ptMinKaon_ || abs(pfCand.eta())>etaMaxKaon_) continue;
		    if(deltaR(muon1, pfCand) < 0.01 || deltaR(muon2, pfCand) < 0.01) continue;

                    pair<double,double> DCA = computeDCA(pfCand,
                                                         bFieldHandle,
                                                         beamSpot);
                    double DCABS = DCA.first;
                    double DCABSErr = DCA.second;

                    if(fabs(DCABS/DCABSErr)<DCASigMinKaon_) continue;
                    
                    RefCountedKinematicVertex refitVertexBToKMuMu;
                    RefCountedKinematicParticle refitBToKMuMu;
                    RefCountedKinematicParticle refitMuon1;
                    RefCountedKinematicParticle refitMuon2;
                    RefCountedKinematicParticle refitKaon;

                    bool passed = BToKMuMuVertexRefitting(muon1, muon2, pfCand,
                                                          theTTBuilder,
                                                          refitVertexBToKMuMu,
                                                          refitBToKMuMu,
                                                          refitMuon1,
                                                          refitMuon2,
                                                          refitKaon);
                    
                    if (!passed) continue;
                    
                    pair<double,double> BToKMuMuLS = computeLS(refitVertexBToKMuMu,beamSpot);
                    double LSBS = BToKMuMuLS.first;
                    double LSBSErr = BToKMuMuLS.second;
                    
                    double BToKMuMuVtx_Chi2 = (double)refitVertexBToKMuMu->chiSquared();
                    double BToKMuMuVtx_CL = TMath::Prob((double)refitVertexBToKMuMu->chiSquared(),
                                                    int(rint(refitVertexBToKMuMu->degreesOfFreedom())));

                    
                    double cosAlpha = computeCosAlpha(refitBToKMuMu,refitVertexBToKMuMu,beamSpot);
                    
                    double mass_err = sqrt(refitBToKMuMu->currentState().kinematicParametersError().matrix()(6,6));

                    
                    pat::CompositeCandidate BToKMuMuCand;
                    BToKMuMuCand.addDaughter( muon1 , "muon1");
                    BToKMuMuCand.addDaughter( muon2 , "muon2");
                    BToKMuMuCand.addDaughter( pfCand, "kaon");
                    BToKMuMuCand.addUserInt("mu1_index", i);
                    BToKMuMuCand.addUserInt("mu2_index", j);
                    BToKMuMuCand.addUserInt("kaon_index", isPFCand ? k : -1);
                    BToKMuMuCand.addUserInt("kaon_lostTrack_index", isPFCand ? -1 : k-pfCandNumber);
                    BToKMuMuCand.addUserInt("kaon_isPFCand", (int)isPFCand);

                    math::XYZVector refitMu1V3D = refitMuon1->refittedTransientTrack().track().momentum();
                    BToKMuMuCand.addUserFloat("mu1_pt",     sqrt(refitMu1V3D.perp2()));
                    BToKMuMuCand.addUserFloat("mu1_eta",    refitMu1V3D.eta());
                    BToKMuMuCand.addUserFloat("mu1_phi",    refitMu1V3D.phi());
                    BToKMuMuCand.addUserInt("mu1_charge", refitMuon1->currentState().particleCharge());

                    math::XYZVector refitMu2V3D = refitMuon2->refittedTransientTrack().track().momentum();
                    BToKMuMuCand.addUserFloat("mu2_pt",     sqrt(refitMu2V3D.perp2()));
                    BToKMuMuCand.addUserFloat("mu2_eta",    refitMu2V3D.eta());
                    BToKMuMuCand.addUserFloat("mu2_phi",    refitMu2V3D.phi());
                    BToKMuMuCand.addUserInt("mu2_charge", refitMuon2->currentState().particleCharge());

                    TLorentzVector muon1cand;
                    muon1cand.SetPtEtaPhiM(sqrt(refitMu1V3D.perp2()), refitMu1V3D.eta(), refitMu1V3D.phi(), MuonMass_);
                    TLorentzVector muon2cand;
                    muon2cand.SetPtEtaPhiM(sqrt(refitMu2V3D.perp2()), refitMu2V3D.eta(), refitMu2V3D.phi(), MuonMass_);
                    BToKMuMuCand.addUserFloat("mumuKFit_mumu_mass", (muon1cand+muon2cand).Mag());

                    math::XYZVector refitKaonV3D = refitKaon->refittedTransientTrack().track().momentum();
                    BToKMuMuCand.addUserFloat("kaon_pt",    sqrt(refitKaonV3D.perp2()));
                    BToKMuMuCand.addUserFloat("kaon_eta",   refitKaonV3D.eta());
                    BToKMuMuCand.addUserFloat("kaon_phi",   refitKaonV3D.phi());
                    BToKMuMuCand.addUserInt("kaon_charge",refitKaon->currentState().particleCharge());
                    BToKMuMuCand.addUserFloat("kaon_DCASig", DCABS/DCABSErr);

                    BToKMuMuCand.addUserInt("mumuRefit", (int)passedDiMuon);
                    BToKMuMuCand.addUserFloat("mumu_pt",    (passedDiMuon)? sqrt(refitMuMuV3D.perp2()) : -1.);
                    BToKMuMuCand.addUserFloat("mumu_eta",   (passedDiMuon)? refitMuMuV3D.eta() : -9.);
                    BToKMuMuCand.addUserFloat("mumu_phi",   (passedDiMuon)? refitMuMuV3D.phi() : -9.);
                    BToKMuMuCand.addUserFloat("mumu_mass",  (passedDiMuon)? refitMuMu->currentState().mass() : -1.);
                    BToKMuMuCand.addUserFloat("mumu_mass_err", (passedDiMuon)? MuMu_mass_err : -1.);
                    BToKMuMuCand.addUserFloat("mumu_Lxy", (passedDiMuon)? (float) MuMuLSBS/MuMuLSBSErr : -1.);
                    BToKMuMuCand.addUserFloat("mumu_ctxy", (passedDiMuon)? (float) MuMuLSBS/sqrt(refitMuMuV3D.perp2()) : -1.);
                    BToKMuMuCand.addUserFloat("mumu_Chi2_vtx", (passedDiMuon)? (float) MuMuVtx_Chi2 : -1.);
                    BToKMuMuCand.addUserFloat("mumu_CL_vtx", (passedDiMuon)? (float) MuMuVtx_CL : -1.);
                    
                    math::XYZVector refitBToKMuMuV3D = refitMu1V3D + refitMu2V3D + refitKaonV3D;
                    BToKMuMuCand.addUserFloat("pt",     sqrt(refitBToKMuMuV3D.perp2()));
                    BToKMuMuCand.addUserFloat("eta",    refitBToKMuMuV3D.eta());
                    BToKMuMuCand.addUserFloat("phi",    refitBToKMuMuV3D.phi());
                    BToKMuMuCand.addUserFloat("mass",   refitBToKMuMu->currentState().mass());
                    BToKMuMuCand.addUserFloat("mass_err", mass_err);
                    
                    BToKMuMuCand.addUserFloat("Lxy", (float) LSBS/LSBSErr);
                    BToKMuMuCand.addUserFloat("ctxy", (float) LSBS/sqrt(refitBToKMuMuV3D.perp2()));
                    BToKMuMuCand.addUserFloat("Chi2_vtx", (float) BToKMuMuVtx_Chi2);
                    BToKMuMuCand.addUserFloat("CL_vtx", (float) BToKMuMuVtx_CL);
                    BToKMuMuCand.addUserFloat("cosAlpha", (float) cosAlpha);

                    float pt_2trk = -9999.;
                    float eta_2trk = -9999.;
                    float phi_2trk = -9999.;
                    float mass_2trk = -9999.;
                    float mass_err_2trk = -9999.;
                    float Lxy_2trk = -9999.;
                    float ctxy_2trk = -9999.;
                    float Chi2_vtx_2trk = -9999.;
                    float CL_vtx_2trk = -9999.;
                    float cosAlpha_2trk = -9999.;

                    //This is crashing for some unexplained reason
                    //https://hypernews.cern.ch/HyperNews/CMS/get/physTools/2746/1/2/1/1.html
                    //Disabled for now

                    //if(save2TrkRefit_ && passed_MuMuRefit){

                    if(0){

                       RefCountedKinematicVertex refitVertexBToKJPsiMuMu;
                       RefCountedKinematicParticle refitBToKJPsiMuMu;
                       RefCountedKinematicParticle refitJPsiMuMu;
                       RefCountedKinematicParticle refitKaon_KJPsi;

                       passed = BToKJPsiMuMuVertexRefitting(refitMuMu, pfCand,
                                                            theTTBuilder,
                                                            refitVertexBToKJPsiMuMu,
                                                            refitBToKJPsiMuMu,
                                                            refitJPsiMuMu,
                                                            refitKaon_KJPsi);

                       if(passed){

                          math::XYZVector refitJPsiMuMuV3D = refitJPsiMuMu->refittedTransientTrack().track().momentum();
                          math::XYZVector refitKaonV3D_KJPsi = refitKaon_KJPsi->refittedTransientTrack().track().momentum();
                          math::XYZVector refitBToKJPsiMuMuV3D = refitJPsiMuMuV3D + refitKaonV3D_KJPsi;

                          pt_2trk = sqrt(refitBToKJPsiMuMuV3D.perp2());
                          eta_2trk = refitBToKJPsiMuMuV3D.eta();
                          phi_2trk = refitBToKJPsiMuMuV3D.phi();
                          mass_2trk = refitBToKJPsiMuMu->currentState().mass();
                          mass_err_2trk = sqrt(refitBToKJPsiMuMu->currentState().kinematicParametersError().matrix()(6,6));

                          pair<double,double> BToKJPsiMuMuLS = computeLS(refitVertexBToKJPsiMuMu,beamSpot);
                          double LSBS_2trk = BToKJPsiMuMuLS.first;
                          double LSBSErr_2trk = BToKJPsiMuMuLS.second;
                          Lxy_2trk = LSBS_2trk/LSBSErr_2trk;
                          ctxy_2trk = LSBS_2trk/pt_2trk;
                          Chi2_vtx_2trk = (double)refitVertexBToKJPsiMuMu->chiSquared();
                          CL_vtx_2trk = TMath::Prob((double)refitVertexBToKJPsiMuMu->chiSquared(),
                                                    int(rint(refitVertexBToKJPsiMuMu->degreesOfFreedom())));
                          cosAlpha_2trk = computeCosAlpha(refitBToKJPsiMuMu,refitVertexBToKJPsiMuMu,beamSpot);

                          }

                    }

                    BToKMuMuCand.addUserFloat("pt_2trk", pt_2trk);
                    BToKMuMuCand.addUserFloat("eta_2trk", eta_2trk);
                    BToKMuMuCand.addUserFloat("phi_2trk", phi_2trk);
                    BToKMuMuCand.addUserFloat("mass_2trk", mass_2trk);
                    BToKMuMuCand.addUserFloat("mass_err_2trk", mass_err_2trk);
                    BToKMuMuCand.addUserFloat("Lxy_2trk", Lxy_2trk);
                    BToKMuMuCand.addUserFloat("ctxy_2trk", ctxy_2trk);
                    BToKMuMuCand.addUserFloat("Chi2_vtx_2trk", Chi2_vtx_2trk);
                    BToKMuMuCand.addUserFloat("CL_vtx_2trk", CL_vtx_2trk);
                    BToKMuMuCand.addUserFloat("cosAlpha_2trk", cosAlpha_2trk);
                    
                    result->push_back(BToKMuMuCand);
                    
                    
                }
                
            }
            
        }
        
    }
    
    iEvent.put(std::move(result));
    
}



bool BToKmumuProducer::MuMuVertexRefitting(const pat::Muon & muon1,
                                           const pat::Muon & muon2,
                                           edm::ESHandle<TransientTrackBuilder> theTTBuilder,
                                           RefCountedKinematicVertex &refitVertex,
                                           RefCountedKinematicParticle &refitMuMu,
                                           RefCountedKinematicParticle &refitMu1,
                                           RefCountedKinematicParticle &refitMu2){
    
    const reco::TransientTrack muon1TT = theTTBuilder->build(muon1.innerTrack());
    const reco::TransientTrack muon2TT = theTTBuilder->build(muon2.innerTrack());
    
    KinematicParticleFactoryFromTransientTrack partFactory;
    KinematicParticleVertexFitter PartVtxFitter;
    
    std::vector<RefCountedKinematicParticle> muonParticles;
    double chi = 0.;
    double ndf = 0.;
    muonParticles.push_back(partFactory.particle(muon1TT,MuonMass_,chi,ndf,MuonMassErr_));
    muonParticles.push_back(partFactory.particle(muon2TT,MuonMass_,chi,ndf,MuonMassErr_));
    RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles);
    
    if ( !mumuVertexFitTree->isValid()) return false;
    
    mumuVertexFitTree->movePointerToTheTop();
    refitVertex = mumuVertexFitTree->currentDecayVertex();
    refitMuMu = mumuVertexFitTree->currentParticle();
    
    if ( !refitVertex->vertexIsValid()) return false;

    // extract the re-fitted tracks
    mumuVertexFitTree->movePointerToTheTop();
    
    mumuVertexFitTree->movePointerToTheFirstChild();
    refitMu1 = mumuVertexFitTree->currentParticle();
    
    mumuVertexFitTree->movePointerToTheNextChild();
    refitMu2 = mumuVertexFitTree->currentParticle();
    
    return true;
    
}




bool BToKmumuProducer::BToKMuMuVertexRefitting(const pat::Muon &muon1,
                                               const pat::Muon &muon2,
                                               const pat::PackedCandidate &kaon,
                                               edm::ESHandle<TransientTrackBuilder> theTTBuilder,
                                               RefCountedKinematicVertex &refitVertex,
                                               RefCountedKinematicParticle &refitBToKMuMu,
                                               RefCountedKinematicParticle &refitMu1,
                                               RefCountedKinematicParticle &refitMu2,
                                               RefCountedKinematicParticle &refitKaon){

    const reco::TransientTrack muon1TT = theTTBuilder->build(muon1.innerTrack());
    const reco::TransientTrack muon2TT = theTTBuilder->build(muon2.innerTrack());
    const reco::TransientTrack kaonTT = theTTBuilder->build(kaon.bestTrack());

    KinematicParticleFactoryFromTransientTrack partFactory;
    KinematicParticleVertexFitter PartVtxFitter;

    std::vector<RefCountedKinematicParticle> BToKMuMuParticles;
    double chi = 0.;
    double ndf = 0.;
    BToKMuMuParticles.push_back(partFactory.particle(muon1TT,MuonMass_,chi,ndf,MuonMassErr_));
    BToKMuMuParticles.push_back(partFactory.particle(muon2TT,MuonMass_,chi,ndf,MuonMassErr_));
    BToKMuMuParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,KaonMassErr_));

    RefCountedKinematicTree BToKMuMuVertexFitTree = PartVtxFitter.fit(BToKMuMuParticles);

    if ( !BToKMuMuVertexFitTree->isValid()) return false;

    BToKMuMuVertexFitTree->movePointerToTheTop();
    refitVertex = BToKMuMuVertexFitTree->currentDecayVertex();
    refitBToKMuMu = BToKMuMuVertexFitTree->currentParticle();

    if ( !refitVertex->vertexIsValid()) return false;

    // extract the re-fitted tracks
    BToKMuMuVertexFitTree->movePointerToTheTop();

    BToKMuMuVertexFitTree->movePointerToTheFirstChild();
    refitMu1 = BToKMuMuVertexFitTree->currentParticle();

    BToKMuMuVertexFitTree->movePointerToTheNextChild();
    refitMu2 = BToKMuMuVertexFitTree->currentParticle();

    BToKMuMuVertexFitTree->movePointerToTheNextChild();
    refitKaon = BToKMuMuVertexFitTree->currentParticle();

    return true;

}





bool BToKmumuProducer::BToKJPsiMuMuVertexRefitting(const RefCountedKinematicParticle refitMuMu,
                                                   const pat::PackedCandidate &kaon,
                                                   edm::ESHandle<TransientTrackBuilder> theTTBuilder,
                                                   RefCountedKinematicVertex &refitVertex,
                                                   RefCountedKinematicParticle &refitBToKJPsiMuMu,
                                                   RefCountedKinematicParticle &refitJPsi,
                                                   RefCountedKinematicParticle &refitKaon){

  const reco::TransientTrack MuMuTT = refitMuMu->refittedTransientTrack();
  const reco::TransientTrack kaonTT = theTTBuilder->build(kaon.bestTrack());

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;

  std::vector<RefCountedKinematicParticle> BToKMuMuParticles;
  double chi = 0.;
  double ndf = 0.;

  float MuMu_mass = refitMuMu->currentState().mass();
  float MuMu_mass_err = sqrt(refitMuMu->currentState().kinematicParametersError().matrix()(6,6));

  if(JPsiMassConstraint_ > 0){
    MuMu_mass = JPsiMassConstraint_;
    MuMu_mass_err = JPsiMassErr_;
  }

  BToKMuMuParticles.push_back(partFactory.particle(MuMuTT,MuMu_mass,chi,ndf,MuMu_mass_err));
  BToKMuMuParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,KaonMassErr_));

  RefCountedKinematicTree BToKMuMuVertexFitTree = PartVtxFitter.fit(BToKMuMuParticles);

  if ( !BToKMuMuVertexFitTree->isValid()) return false;

  BToKMuMuVertexFitTree->movePointerToTheTop();
  refitVertex = BToKMuMuVertexFitTree->currentDecayVertex();
  refitBToKJPsiMuMu = BToKMuMuVertexFitTree->currentParticle();

  if ( !refitVertex->vertexIsValid()) return false;

  // extract the re-fitted tracks
  BToKMuMuVertexFitTree->movePointerToTheTop();

  BToKMuMuVertexFitTree->movePointerToTheFirstChild();
  refitJPsi = BToKMuMuVertexFitTree->currentParticle();

  BToKMuMuVertexFitTree->movePointerToTheNextChild();
  refitKaon = BToKMuMuVertexFitTree->currentParticle();

  return true;



}







pair<double,double> BToKmumuProducer::computeLS(RefCountedKinematicVertex refitVertex,
                                                reco::BeamSpot beamSpot){
    
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
    
  double LSBS = sqrt(v.Norm2Sqr());
  double LSBSErr = sqrt( v*(errVtx*v) + v*(errBS*v) ) / LSBS;
    
  pair<double,double> LS = make_pair(LSBS,LSBSErr);
    
  return LS;
}



double BToKmumuProducer::computeCosAlpha(RefCountedKinematicParticle refitBToKMuMu,
                                         RefCountedKinematicVertex refitVertex,
                                         reco::BeamSpot beamSpot){
    
  TVector v(2);
  v[0] = refitVertex->position().x()-beamSpot.position().x();
  v[1] = refitVertex->position().y()-beamSpot.position().y();

  TVector w(2);
  w[0] = refitBToKMuMu->currentState().globalMomentum().x();
  w[1] = refitBToKMuMu->currentState().globalMomentum().y();
    
  double cosAlpha = v*w/sqrt(v.Norm2Sqr()*w.Norm2Sqr());
  return cosAlpha;
}






pair<double,double> BToKmumuProducer::computeDCA(const pat::PackedCandidate &kaon,
                                                 edm::ESHandle<MagneticField> bFieldHandle,
                                                 reco::BeamSpot beamSpot){

  const reco::TransientTrack trackTT((*(kaon.bestTrack())), &(*bFieldHandle));

  TrajectoryStateClosestToPoint theDCAXBS = trackTT.trajectoryStateClosestToPoint( GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()) );  
  
  double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
  double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
  pair<double,double> DCA = make_pair(DCABS,DCABSErr);
    
  return DCA;
}




DEFINE_FWK_MODULE(BToKmumuProducer);
