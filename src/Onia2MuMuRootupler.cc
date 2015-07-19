// -*- C++ -*-
//
// Package:    Onia2MuMuRootupler
// Class:      Onia2MuMuRootupler
// 
// Description: Dump  Onia(mu+ mu-)  decays
//
// Author:  Alberto Sanchez Hernandez
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"

//
// class declaration
//

class Onia2MuMuRootupler:public edm::EDAnalyzer {
      public:
	explicit Onia2MuMuRootupler(const edm::ParameterSet &);
	~Onia2MuMuRootupler();

	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

      private:
        UInt_t getTriggerBits(const edm::Event &);
        bool   isAncestor(const reco::Candidate *, const reco::Candidate *);
        const  reco::Candidate* GetAncestor(const reco::Candidate *);

	virtual void beginJob();
	virtual void analyze(const edm::Event &, const edm::EventSetup &);
	virtual void endJob();

	virtual void beginRun(edm::Run const &, edm::EventSetup const &);
	virtual void endRun(edm::Run const &, edm::EventSetup const &);
	virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
	virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);

	// ----------member data ---------------------------
	std::string file_name;
	edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_Label;
        edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
        edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
        int  pdgid_;
        std::vector<double> OniaMassCuts_;
	bool isMC_;
        bool OnlyBest_;
        bool OnlyGen_;

	UInt_t run;
	UInt_t event;
        Int_t  irank;
        UInt_t trigger;
        Int_t  charge; 

	TLorentzVector dimuon_p4;
	TLorentzVector muonP_p4;
	TLorentzVector muonN_p4;

        Float_t MassErr;
        Float_t vProb;
        Float_t DCA;
        Float_t ppdlPV;
        Float_t ppdlErrPV;
        Float_t ppdlBS;
        Float_t ppdlErrBS;
        Float_t cosAlpha;
        Float_t lxyPV;
        Float_t lxyBS;

	Int_t numPrimaryVertices;

	TTree *onia_tree;

        Int_t mother_pdgId;
        Int_t dimuon_pdgId;
	TLorentzVector gen_dimuon_p4;
	TLorentzVector gen_muonP_p4;
	TLorentzVector gen_muonM_p4;
          
        edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
        edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;

};

//
// constructors and destructor
//

Onia2MuMuRootupler::Onia2MuMuRootupler(const edm::ParameterSet & iConfig):
dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuons"))),
primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
pdgid_(iConfig.getParameter<uint32_t>("onia_pdgid")),
OniaMassCuts_(iConfig.getParameter<std::vector<double>>("onia_mass_cuts")),
isMC_(iConfig.getParameter<bool>("isMC")),
OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
OnlyGen_(iConfig.getParameter<bool>("OnlyGen"))
{
  edm::Service < TFileService > fs;
  onia_tree = fs->make < TTree > ("oniaTree", "Tree of Onia2MuMu");

  if (!OnlyGen_) {
    onia_tree->Branch("run",     &run,     "run/i");
    onia_tree->Branch("event",   &event,   "event/i");
    onia_tree->Branch("irank",   &irank,   "irank/I");
    onia_tree->Branch("trigger", &trigger, "trigger/i");
    onia_tree->Branch("charge",  &charge,  "charge/I");

    onia_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
    onia_tree->Branch("muonP_p4",  "TLorentzVector", &muonP_p4);
    onia_tree->Branch("muonN_p4",  "TLorentzVector", &muonN_p4);

    onia_tree->Branch("MassErr",   &MassErr,    "MassErr/F");
    onia_tree->Branch("vProb",     &vProb,      "vProb/F");
    onia_tree->Branch("DCA",       &DCA,        "DCA/F");
    onia_tree->Branch("ppdlPV",    &ppdlPV,     "ppdlPV/F");
    onia_tree->Branch("ppdlErrPV", &ppdlErrPV,  "ppdlErrPV/F");
    onia_tree->Branch("ppdlBS",    &ppdlBS,     "ppdlBS/F");
    onia_tree->Branch("ppdlErrBS", &ppdlErrBS,  "ppdlErrBS/F");
    onia_tree->Branch("cosAlpha",  &cosAlpha,   "cosAlpha/F");
    onia_tree->Branch("lxyPV",     &lxyPV,      "lxyPV/F");
    onia_tree->Branch("lxyBS",     &lxyBS,      "lxyBS/F");

    onia_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
  }

  if (isMC_ || OnlyGen_) {
     onia_tree->Branch("mother_pdgId",  &mother_pdgId,     "mother_pdgId/I");
     onia_tree->Branch("dimuon_pdgId",  &dimuon_pdgId,     "dimuon_pdgId/I");
     onia_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
     onia_tree->Branch("gen_muonP_p4",  "TLorentzVector",  &gen_muonP_p4);
     onia_tree->Branch("gen_muonN_p4",  "TLorentzVector",  &gen_muonM_p4);
  }
  genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
  packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
}

Onia2MuMuRootupler::~Onia2MuMuRootupler() {}

//
// member functions
//

const reco::Candidate* Onia2MuMuRootupler::GetAncestor(const reco::Candidate* p) {
   if (p->numberOfMothers()) {
      if  ((p->mother(0))->pdgId() == p->pdgId()) return GetAncestor(p->mother(0));
      else return p->mother(0);
   }
   //std::cout << "GetAncestor: Inconsistet ancestor, particle does not have a mother " << std::endl;
   return p;
}

//Check recursively if any ancestor of particle is the given one
bool Onia2MuMuRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}

/* Grab Trigger information. Save it in variable trigger, trigger is an int between 0 and 127, in binary it is:
   (pass 2)(pass 1)(pass 0)
   ex. 7 = pass 0, 1 and 2
   ex. 6 = pass 2, 3
   ex. 1 = pass 0
*/

UInt_t Onia2MuMuRootupler::getTriggerBits(const edm::Event& iEvent ) {
   UInt_t itrigger = 0;
   edm::Handle<edm::TriggerResults> triggerResults_handle;
   iEvent.getByToken(triggerResults_Label, triggerResults_handle);
   if ( triggerResults_handle.isValid() ) { 
      const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
      std::vector <unsigned int> bits_0, bits_1, bits_2, bits_3, bits_4, bits_5, bits_6, bits_7;
      for ( int version = 1; version<3; version ++ ) {
         std::stringstream ss0,ss1,ss2,ss3,ss4,ss5,ss6,ss7;
         ss0<<"HLT_Dimuon16_Jpsi_v"<<version;
         bits_0.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss0.str()).label().c_str()));
         ss1<<"HLT_Dimuon13_PsiPrime_v"<<version;
         bits_1.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss1.str()).label().c_str()));
         ss2<<"HLT_Dimuon13_Upsilon_v"<<version;
         bits_2.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss2.str()).label().c_str()));
         ss3<<"HLT_Dimuon10_Jpsi_Barrel_v"<<version;
         bits_3.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss3.str()).label().c_str()));
         ss4<<"HLT_Dimuon8_PsiPrime_Barrel_v"<<version;
         bits_4.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss4.str()).label().c_str()));
         ss5<<"HLT_Dimuon8_Upsilon_Barrel_v"<<version;
         bits_5.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss5.str()).label().c_str()));
         ss6<<"HLT_Dimuon20_Jpsi_v"<<version;
         bits_6.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss6.str()).label().c_str()));
         ss7<<"HLT_Dimuon0_Phi_Barrel_v"<<version;
         bits_7.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss7.str()).label().c_str()));
      }
      for (unsigned int i=0; i<bits_0.size(); i++) {
         unsigned int bit = bits_0[i];
         if ( bit < triggerResults_handle->size() ){
	   if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 1;
             break;
           }
         }
      }
      for (unsigned int i=0; i<bits_1.size(); i++) {
         unsigned int bit = bits_1[i];
         if ( bit < triggerResults_handle->size() ){
           if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 2;
             break;
           }
         }
      }
      for (unsigned int i=0; i<bits_2.size(); i++) {
         unsigned int bit = bits_2[i];
         if ( bit < triggerResults_handle->size() ){
           if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 4;
             break;
           }
         }
      }
      for (unsigned int i=0; i<bits_3.size(); i++) {
         unsigned int bit = bits_3[i];
         if ( bit < triggerResults_handle->size() ){
           if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 8;
             break;
           }
         }
      }
      for (unsigned int i=0; i<bits_4.size(); i++) {
         unsigned int bit = bits_4[i];
         if ( bit < triggerResults_handle->size() ){
           if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 16;
             break;
           }
         }
      }
      for (unsigned int i=0; i<bits_5.size(); i++) {
         unsigned int bit = bits_5[i];
         if ( bit < triggerResults_handle->size() ){
           if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 32;
             break;
           }
         }
      }
      for (unsigned int i=0; i<bits_6.size(); i++) {
         unsigned int bit = bits_6[i];
         if ( bit < triggerResults_handle->size() ){
           if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 64;
             break;
           }
         }
      }
      for (unsigned int i=0; i<bits_7.size(); i++) {
         unsigned int bit = bits_7[i];
         if ( bit < triggerResults_handle->size() ){
           if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
             itrigger += 128;
             break;
           }
         }
      }
   }
   return itrigger;
}

// ------------ method called for each event  ------------
void Onia2MuMuRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  iEvent.getByToken(dimuon_Label,dimuons);

  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);
 
  if (!OnlyGen_) {
    numPrimaryVertices = -1;
    if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();

    trigger = getTriggerBits(iEvent);
    run     = iEvent.id().run();
    event   = iEvent.id().event();
  }

  dimuon_pdgId = 0;
  mother_pdgId = 0;
  irank = 0;

  // Pruned particles are the one containing "important" stuff
  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genCands_, pruned);

  // Packed particles are all the status 1, so usable to remake jets
  // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packCands_,  packed);

  if ((isMC_ || OnlyGen_) && packed.isValid() && pruned.isValid()) {
     dimuon_pdgId  = 0;
     gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
     int foundit   = 0;

     for (size_t i=0; i<pruned->size(); i++) {
         int p_id = abs((*pruned)[i].pdgId());
         const reco::Candidate *aonia = &(*pruned)[i];
         if (( p_id == pdgid_ ) && (aonia->status() == 2)) {
            dimuon_pdgId = p_id;
            foundit++;
            for (size_t j=0; j<packed->size(); j++) { //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
               const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
               const reco::Candidate * d = &(*packed)[j];
               if ( motherInPrunedCollection != nullptr && (d->pdgId() == 13 ) && isAncestor(aonia , motherInPrunedCollection) ){
                  gen_muonM_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                  foundit++;
               } 
               if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(aonia , motherInPrunedCollection) ) {
                  gen_muonP_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
                  foundit++;
               }
               if ( foundit == 3 ) break;               
            }
            if ( foundit == 3 ) {
               gen_dimuon_p4 = gen_muonM_p4 + gen_muonP_p4;   // this should take into account FSR
               mother_pdgId  = GetAncestor(aonia)->pdgId();
               break;
            } else {
               foundit = 0;
               dimuon_pdgId = 0;
               mother_pdgId = 0;
               gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
            }            
         }  // if ( p_id
     } // for (size

     // sanity check
     //if ( ! dimuon_pdgId ) std::cout << "Onia2MuMuRootupler: does not found the given decay " << run << "," << event << std::endl;
  }  // end if isMC

  float OniaMassMax_ = OniaMassCuts_[1];
  float OniaMassMin_ = OniaMassCuts_[0];

  if ( ! OnlyGen_ && dimuons.isValid() && dimuons->size() > 0) {
     for (pat::CompositeCandidateCollection::const_iterator  dimuonCand = dimuons->begin(); dimuonCand!= dimuons->end(); ++dimuonCand) {
        if (dimuonCand->mass() > OniaMassMin_ && dimuonCand->mass() < OniaMassMax_ && dimuonCand->charge() == 0) {
           dimuon_p4.SetPtEtaPhiM(dimuonCand->pt(), dimuonCand->eta(), dimuonCand->phi(), dimuonCand->mass());
           reco::Candidate::LorentzVector vP = dimuonCand->daughter("muon1")->p4();
           reco::Candidate::LorentzVector vM = dimuonCand->daughter("muon2")->p4();
           if ( dimuonCand->daughter("muon1")->charge() < 0) {
              vP = dimuonCand->daughter("muon2")->p4();
              vM = dimuonCand->daughter("muon1")->p4();
           }
           muonP_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
           muonN_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
           MassErr = dimuonCand->userFloat("MassErr");
           vProb = dimuonCand->userFloat("vProb");
           DCA = dimuonCand->userFloat("DCA");
           ppdlPV = dimuonCand->userFloat("ppdlPV");
           ppdlErrPV = dimuonCand->userFloat("ppdlErrPV");
           ppdlBS = dimuonCand->userFloat("ppdlBS");
           ppdlErrBS = dimuonCand->userFloat("ppdlErrBS");
           cosAlpha = dimuonCand->userFloat("cosAlpha");
           charge = dimuonCand->charge(); 
           TVector3 pperp(dimuonCand->px(), dimuonCand->py(), 0);
           lxyPV = ppdlPV * pperp.Perp() / dimuonCand->mass();
           lxyBS = ppdlBS * pperp.Perp() / dimuonCand->mass();
           irank++;
           onia_tree->Fill();
           if (OnlyBest_) break;
        } 
     }
  } else {
     if (dimuon_pdgId && (OnlyGen_ || isMC_)) onia_tree->Fill();
     //else std::cout << "Onia2MuMuRootupler: does not find a valid dimuon combination " << run << "," << event << std::endl;
  }

  dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  dimuon_pdgId = 0;
  mother_pdgId = 0;
  vProb=-1;
}

// ------------ method called once each job just before starting event loop  ------------
void Onia2MuMuRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void Onia2MuMuRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void Onia2MuMuRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void Onia2MuMuRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void Onia2MuMuRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void Onia2MuMuRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Onia2MuMuRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Onia2MuMuRootupler);
