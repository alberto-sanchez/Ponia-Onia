// -*- C++ -*-
//
// Package:    OniaMM
// Class:      OniaMM
// 
/*
 Description: Dump generation level (genParticles) content for Onia -> mu+ mu- 
              
 Author: Alberto Sanchez Hernandez (asanchez), September 2014
*/

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include <vector>

//
// class declaration
//

class OniaMM : public edm::EDAnalyzer {
public:
  explicit OniaMM(const edm::ParameterSet&);
  ~OniaMM();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  const reco::Candidate* GetStableParticle(const reco::Candidate*);
  const reco::Candidate* GetAncestor(const reco::Candidate*);

  int  pdgid_;
  
  // Variables
  UInt_t    run;
  ULong64_t event;
  UInt_t    lumiblock;
  TLorentzVector gen_muonP_p4;
  TLorentzVector gen_muonM_p4;
  TLorentzVector gen_dimuon_p4;
  Int_t mother_pdgId,dimuon_pdgId;
  TTree* onia_tree;
  
};

//
// constructors and destructor
//

OniaMM::OniaMM(const edm::ParameterSet& iConfig):
pdgid_(iConfig.getParameter<uint32_t>("onia_pdgid"))
{

  edm::Service < TFileService > fs;
  onia_tree = fs->make < TTree > ("oniaTree", "Tree of Onia2MuMu");
  onia_tree->Branch("run",      &run,      "run/i");
  onia_tree->Branch("event",    &event,    "event/l");
  onia_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

  std::cout << "Onia2MuMuRootupler::Onia2MuMuRootupler: Onia id " << pdgid_ << std::endl;
  onia_tree->Branch("mother_pdgId",  &mother_pdgId,     "mother_pdgId/I");
  onia_tree->Branch("dimuon_pdgId",  &dimuon_pdgId,     "dimuon_pdgId/I");
  onia_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
  onia_tree->Branch("gen_muonP_p4",  "TLorentzVector",  &gen_muonP_p4);
  onia_tree->Branch("gen_muonN_p4",  "TLorentzVector",  &gen_muonM_p4);
}

OniaMM::~OniaMM() {}

//
// member functions
//

// ------------ method called for each event  ------------

const reco::Candidate* OniaMM::GetAncestor(const reco::Candidate* p) {
   if (p->numberOfMothers()) {
      if  ((p->mother(0))->pdgId() == p->pdgId()) return GetAncestor(p->mother(0));
      else return p->mother(0);
   }
   std::cout << "GetAncestor: Inconsistet ancestor, particle does not have a mother " << std::endl;
   return p;
}

const reco::Candidate* OniaMM::GetStableParticle(const reco::Candidate* p) {
   if (p->status() == 1) return p;
   int n = p->numberOfDaughters();
   if ( n > 0 ) {
     for (int j = 0; j < n; ++j) {
       const  reco::Candidate* d = p->daughter(j);
       if (d->pdgId() == p->pdgId()) return GetStableParticle(d);
     }
   } 
   std::cout << "GetStableParticle: Inconsistent state of particle, it has not daugthers, but is unstable " << std::endl;
   return p;  
}

void OniaMM::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::GenParticleCollection> GenParticles;
  iEvent.getByLabel("genParticles",GenParticles);

  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  dimuon_pdgId = 0;
  mother_pdgId = 0;

  gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muonP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muonM_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  if ( GenParticles.isValid() ) {

     for ( reco::GenParticleCollection::const_iterator itParticle = GenParticles->begin(); itParticle != GenParticles->end(); ++itParticle ) {
         Int_t pdgId = itParticle->pdgId();
         dimuon_pdgId = 0;
         int foundit = 0;
         if ( (abs(pdgId) == pdgid_) && (itParticle->status() == 2) ) {
            int n = itParticle->numberOfDaughters();
            if (n < 2) continue;

            const  reco::Candidate* onia = &(*itParticle);
            foundit++;
            dimuon_pdgId = pdgId;

            bool yetM = false;
            bool yetP = false;

            for (int j = 0; j < n; ++ j) {
                const  reco::Candidate* d = itParticle->daughter(j);
                Int_t dauId = d->pdgId();
                if ( dauId == 13 && !yetM) {
                   const  reco::Candidate* mM = GetStableParticle(d);
                   gen_muonM_p4.SetPtEtaPhiM(mM->pt(),mM->eta(),mM->phi(),mM->mass());
                   foundit++;
                   yetM = true;
                } 
                if ( dauId == -13 && !yetP) {
                   const  reco::Candidate* mP = GetStableParticle(d);
                   gen_muonP_p4.SetPtEtaPhiM(mP->pt(),mP->eta(),mP->phi(),mP->mass());
                   foundit++;
                   yetP = true;
                } 
            }
            if ( foundit == 3 ) {
               mother_pdgId = GetAncestor(onia)->pdgId();
               gen_dimuon_p4 = gen_muonM_p4 + gen_muonP_p4;
               break;
            } else {
               foundit = 0;
               dimuon_pdgId = 0;
               mother_pdgId = 0;
            }
         }  // if ( pdg

      }   // for ( reco
  }     // if (GenPar
  // sanity check
  if ( ! dimuon_pdgId ) std::cout << "OniaMM: Decay not found [" << iEvent.id().run() << "," << iEvent.id().event() << "]" << std::endl; 
  else onia_tree ->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void OniaMM::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void OniaMM::endJob() {}

// ------------ method called when starting to processes a run  ------------
void OniaMM::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void OniaMM::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void OniaMM::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void OniaMM::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void OniaMM::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(OniaMM);
