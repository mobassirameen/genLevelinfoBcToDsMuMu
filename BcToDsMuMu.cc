

// -*- C++ -*-
//
// Package:    BcToDsMuMu
// Class:      BcToDsMuMu
// 
//=================================================
// Original by: Chandiprasad Kar                  |
//<chandiprasad.kar@cern.ch>                      |
//Major variables are added 05 Sept, 2020         |
//=================================================

// system include files
#include <memory>

// user include files
#include "bctodsmumu_analysis/BcToDsMuMuPAT/src/BcToDsMuMu.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"

#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include <vector>
#include <utility>
#include <string>
#include <iostream>
//
// constants, enums and typedefs
//

typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//



//
// constructors and destructor
//
BcToDsMuMu::BcToDsMuMu(const edm::ParameterSet& iConfig):

  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  triggerPrescales_ (consumes<pat::PackedTriggerPrescales> (iConfig.getParameter<edm::InputTag>("prescales"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),  
  trigTable_( iConfig.getParameter<std::vector<std::string> >("TriggerNames")),   
  puToken_(consumes<std::vector< PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PuInfoTag"))), 

  v0PtrCollection_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secondaryVerticesPtr"))),	       

  genParticles_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))),
  packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter <edm::InputTag> ("packedGenParticles"))),
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  doMC_ ( iConfig.getUntrackedParameter<bool>("doMC",false) ),

  tree_(0), 

  mumC2(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupNHits(0), mupNPHits(0),
  mumdxy(0), mupdxy(0), mumdz(0), mupdz(0),
  muon_dca(0),

  tri_Dim25(0), tri_JpsiTk(0), tri_JpsiTkTk(0),
  dr0(0),dr1(0),
  dpt0(0),dpt1(0),
  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0), 
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),

  nVtx(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),
  
  //************************ ****************************************************
  
  //*******************************************************
  nB(0), nMu(0),
  B_mass(0), B_px(0), B_py(0), B_pz(0),
  
  B_Ds_mass(0), B_Ds_px(0), B_Ds_py(0), B_Ds_pz(0),
  B_Ds_pt1(0), B_Ds_px1(0), B_Ds_py1(0), B_Ds_pz1(0), 
  B_Ds_pt2(0), B_Ds_px2(0), B_Ds_py2(0), B_Ds_pz2(0), 

  B_Ds_px1_track(0), B_Ds_py1_track(0), B_Ds_pz1_track(0), 
  B_Ds_px2_track(0), B_Ds_py2_track(0), B_Ds_pz2_track(0), 
  
  pi1dxy(0), pi2dxy(0), pi1dz(0), pi2dz(0),
  pi1dxy_e(0), pi2dxy_e(0), pi1dz_e(0), pi2dz_e(0),
  B_Ds_charge1(0), B_Ds_charge2(0),

  B_mumu_mass(0), B_mumu_px(0), B_mumu_py(0), B_mumu_pz(0),
  B_mumu_pt1(0), B_mumu_px1(0), B_mumu_py1(0), B_mumu_pz1(0), 
  B_mumu_pt2(0), B_mumu_px2(0), B_mumu_py2(0), B_mumu_pz2(0), 
  B_mumu_charge1(0), B_mumu_charge2(0),

  B_Ds_chi2(0), B_mumu_chi2(0), B_chi2(0), B_chi2dof(0),
  B_Prob(0), B_mumu_Prob(0), B_Ds_Prob(0),
  bDecayVtxX(0), bDecayVtxY(0), bDecayVtxZ(0), bDecayVtxXE(0), bDecayVtxYE(0), bDecayVtxZE(0),
  bDecayVtxXYE(0), bDecayVtxXZE(0), bDecayVtxYZE(0),
  
  VDecayVtxX(0), VDecayVtxY(0), VDecayVtxZ(0), VDecayVtxXE(0), VDecayVtxYE(0), VDecayVtxZE(0),
  VDecayVtxXYE(0), VDecayVtxXZE(0), VDecayVtxYZE(0), 
  pVtxIPX(0),  pVtxIPY(0),  pVtxIPZ(0),
  pVtxIPXE(0),  pVtxIPYE(0),  pVtxIPZE(0),  pVtxIPCL(0),
  pVtxIPXYE(0),  pVtxIPXZE(0),  pVtxIPYZE(0),
  
  B_l3d(0),  B_l3dE(0),  B_lxy(0), B_lxyE(0),
  B_cosalpha(0),   B_cosalphaxy(0), alpha(0),  B_treco(0),   B_trecoe(0),  B_trecoxy(0), B_trecoxye(0),
  B_pvip(0), B_pviperr(0), B_pvips(0), B_pvlzip(0), B_pvlziperr(0), B_pvlzips(0),
  B_pv2ip(0), B_pv2iperr(0), B_pv2ips(0), B_pv2lzip(0), B_pv2lziperr(0), B_pv2lzips(0),

  B_l3d_pv2(0),  B_l3dE_pv2(0),
  B_iso(0), B_mum_iso(0), B_mup_iso(0), B_pi1_iso(0),B_pi2_iso(0),
  
  istruemum(0), istruemup(0), istruekp(0), istruekm(0), istruebs(0),
  bunchXingMC(0), numInteractionsMC(0), trueNumInteractionsMC(0),
  run(0), event(0),
  lumiblock(0)
  
{
   //now do what ever initialization is needed
}

BcToDsMuMu::~BcToDsMuMu()
{

}

// ------------ method called to for each event  ------------
void BcToDsMuMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  //*********************************
  // Get event content information
  //*********************************  
  
  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);

  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genParticles_, pruned);

  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packedGenToken_,packed);

  //***************************************
  // Get Gen Level information
  //***************************************
  gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_ks_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_b_vtx.SetXYZ(0.,0.,0.);
  gen_jpsi_vtx.SetXYZ(0.,0.,0.);
  gen_b_ct = -9999.;
  
  if ( (isMC_ || OnlyGen_) && pruned.isValid() ) {
    int foundit = 0;
    for (size_t i=0; i<pruned->size(); i++) {
      foundit = 0;
      const reco::Candidate *dau = &(*pruned)[i];
      if ( (abs(dau->pdgId()) == 531) ) { //&& (dau->status() == 2) ) 
	foundit++;
	gen_b_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
	gen_b_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
	for (size_t k=0; k<dau->numberOfDaughters(); k++) {
	  const reco::Candidate *gdau = dau->daughter(k);
	  if (gdau->pdgId()==443 ) { //&& gdau->status()==2) {
	    foundit++;
	    gen_jpsi_vtx.SetXYZ(gdau->vx(),gdau->vy(),gdau->vz());
	    gen_b_ct = GetLifetime(gen_b_p4,gen_b_vtx,gen_jpsi_vtx);
	    int nm=0;
	    for (size_t l=0; l<gdau->numberOfDaughters(); l++) {
	      const reco::Candidate *mm = gdau->daughter(l);
	      if (mm->pdgId()==13) { foundit++;
		if (mm->status()!=1) {
		  for (size_t m=0; m<mm->numberOfDaughters(); m++) {
		    const reco::Candidate *mu = mm->daughter(m);
		    if (mu->pdgId()==13 ) { //&& mu->status()==1) {
		      nm++;
		      gen_muon1_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
		      break;
		    }
		  }
		} else {
		  gen_muon1_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
		  nm++;
		}
	      }
	      if (mm->pdgId()==-13) { foundit++;
		if (mm->status()!=1) {
		  for (size_t m=0; m<mm->numberOfDaughters(); m++) {
		    const reco::Candidate *mu = mm->daughter(m);
		    if (mu->pdgId()==-13 ) { //&& mu->status()==1) {
		      nm++;
		      gen_muon2_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
		      break;
		    }
		  }
		} else {
		  gen_muon2_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
		  nm++;
		}
	      }
	    }
	    if (nm==2) gen_jpsi_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
	    else foundit-=nm;
	  }
	  for (size_t lk=0; lk<packed->size(); lk++) {
	    const reco::Candidate * dauInPrunedColl = (*packed)[lk].mother(0);
	    int stable_id = (*packed)[lk].pdgId();
	    if (dauInPrunedColl != nullptr && isAncestor(gdau,dauInPrunedColl)) {
	      if(stable_id == 211) {foundit++;
		gen_pion1_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
	      }
	      if(stable_id == -211){ foundit++;
		gen_pion2_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
	      }
	    }
	  }

        } // for (size_t k                                                                                                                                                          
      }   // if (abs(dau->pdgId())==531 )                                                                                                                                           
      if (foundit>=7) break;
    } // for i                                                                                                                                                                      
    if (foundit!=7) {
      gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_ks_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_b_vtx.SetXYZ(0.,0.,0.);
      gen_jpsi_vtx.SetXYZ(0.,0.,0.);
      gen_b_ct = -9999.;
      std::cout << "Does not found the given decay " << run << "," << event << " foundit=" << foundit << std::endl; // sanity check
    }
  }

  nB = 0; nMu = 0;
  if ( OnlyGen_ ) { 
    tree_->Fill();
    return;
  }

  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event(); 
  
		if(nB==0){		    
		  //nMu  = nMu_tmp;
		  //cout<< "*Number of Muons : " << nMu_tmp << endl;
		} // end nB==0		     
		
			      
		nB++;	
		
		//if(isMC_)saveTruthMatch(iEvent);
		
		
		
	     // }
	   // }
	    
	 // }
//	}
   // } ye extra hai ya nai confusion hai
  
  
   //fill the tree and clear the vectors
  if (nB > 0 ) 
    {
      std::cout << "filling tree " << nB<<std::endl;
      tree_->Fill();
    }
  // *********
  
  nB = 0; //nMu = 0;
  
  
}

bool BcToDsMuMu::IsTheSame(const reco::Track& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
bool BcToDsMuMu::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
  if (ancestor == particle ) return true;
  for (size_t i=0; i< particle->numberOfMothers(); i++) {
    if (isAncestor(ancestor,particle->mother(i))) return true;
  }
  return false;
}

double BcToDsMuMu::GetLifetime(TLorentzVector b_p4, TVector3 production_vtx, TVector3 decay_vtx) {
  TVector3 pv_dv = decay_vtx - production_vtx;
  TVector3 b_p3  = b_p4.Vect();
  pv_dv.SetZ(0.);
  b_p3.SetZ(0.);
  Double_t lxy   = pv_dv.Dot(b_p3)/b_p3.Mag();
  return lxy*b_p4.M()/b_p3.Mag();
}

// ------------ method called once each job just before starting event loop  ------------

void 
BcToDsMuMu::beginJob()
{

  //std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;
  std::cout << "Beginning analyzer job with value of isMC_ = " << isMC_ << std::endl;
  edm::Service<TFileService> fs;

  tree_ = fs->make<TTree>("ntuple","Bs->J/psi Ds ntuple");

  tree_->Branch("nB",&nB,"nB/i");
  //tree_->Branch("nMu",&nMu,"nMu/i");


  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");

  
  //gen information
  if (isMC_) {
    tree_->Branch("gen_b_p4",     "TLorentzVector",  &gen_b_p4);
    tree_->Branch("gen_ks_p4",   "TLorentzVector",  &gen_ks_p4);
    tree_->Branch("gen_pion1_p4",  "TLorentzVector",  &gen_pion1_p4);
    tree_->Branch("gen_pion2_p4",  "TLorentzVector",  &gen_pion2_p4);
    tree_->Branch("gen_jpsi_p4",   "TLorentzVector",  &gen_jpsi_p4);
    tree_->Branch("gen_muon1_p4",  "TLorentzVector",  &gen_muon1_p4);
    tree_->Branch("gen_muon2_p4",  "TLorentzVector",  &gen_muon2_p4);
    tree_->Branch("gen_b_vtx",    "TVector3",        &gen_b_vtx);
    tree_->Branch("gen_jpsi_vtx",  "TVector3",        &gen_jpsi_vtx);
    tree_->Branch("gen_b_ct",     &gen_b_ct,        "gen_b_ct/F");
    }


}
// ------------ method called once each job just after ending the event loop  ------------
void BcToDsMuMu::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(BcToDsMuMu);
