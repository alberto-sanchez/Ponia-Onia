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
  
  // Variables
  TLorentzVector muonP_p4;
  TLorentzVector muonM_p4;
  TLorentzVector dimuon_p4;
  Int_t onia_id, mom_id;
  TTree* gen_tree;
  
};

//
// constructors and destructor
//

OniaMM::OniaMM(const edm::ParameterSet& iConfig) {
	edm::Service<TFileService> fs;
	gen_tree = fs->make<TTree>("GenTree","Tree of genParticles");
        gen_tree->Branch("mother_pdgId",  &mom_id,          "mom_id/I");
        gen_tree->Branch("dimuon_pdgId",  &onia_id,         "onia_id/I");
	gen_tree->Branch("gen_muonP_p4",  "TLorentzVector", &muonP_p4);
	gen_tree->Branch("gen_muonM_p4",  "TLorentzVector", &muonM_p4);
        gen_tree->Branch("gen_dimuon_p4", "TLorentzVector", &dimuon_p4);
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
  
  onia_id = 0;
  mom_id = 0;

  if ( GenParticles.isValid() ) {

     for ( reco::GenParticleCollection::const_iterator itParticle = GenParticles->begin(); itParticle != GenParticles->end(); ++itParticle ) {
         Int_t pdgId = itParticle->pdgId();
         onia_id = 0;
         int foundit = 0;
         if ( (pdgId == 443 || pdgId == 100443 || pdgId == 553 || pdgId == 100553 || pdgId == 200553)
              && itParticle->status() == 2
            ) { 
            int n = itParticle->numberOfDaughters();
            if (n < 2) continue;

            const  reco::Candidate* onia = &(*itParticle);
            foundit++;
            onia_id = pdgId;

            bool yetM = false;
            bool yetP = false;

            for (int j = 0; j < n; ++ j) {
                const  reco::Candidate* d = itParticle->daughter(j);
                Int_t dauId = d->pdgId();
                if ( dauId == 13 && !yetM) {
                   const  reco::Candidate* mM = GetStableParticle(d);
                   muonM_p4.SetPtEtaPhiM(mM->pt(),mM->eta(),mM->phi(),mM->mass());
                   foundit++;
                   yetM = true;
                } 
                if ( dauId == -13 && !yetP) {
                   const  reco::Candidate* mP = GetStableParticle(d);
                   muonP_p4.SetPtEtaPhiM(mP->pt(),mP->eta(),mP->phi(),mP->mass());
                   foundit++;
                   yetP = true;
                } 
            }
            if ( foundit == 3 ) {
               mom_id = GetAncestor(onia)->pdgId();
               dimuon_p4 = muonM_p4 + muonP_p4;
               break;
            } else {
               foundit = 0;
               onia_id = 0;
               mom_id = 0;
            }
         }  // if ( pdg

      }   // for ( reco
  }     // if (GenPar
  // sanity check
  if ( ! onia_id ) std::cout << "OniaMM: Decay not found [" << iEvent.id().run() << "," << iEvent.id().event() << "]" << std::endl; 
  else gen_tree ->Fill();
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
