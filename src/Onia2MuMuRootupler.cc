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
        const reco::Candidate* GetAncestor(const reco::Candidate* p);
        bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);

      private:
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
        int  pdgid_;
        std::vector<double> OniaMassCuts_;
	bool isMC_;
        bool OnlyBest_;

	UInt_t run;
	UInt_t event;
        Int_t  irank;

	TLorentzVector dimuon_p4;
	TLorentzVector muonP_p4;
	TLorentzVector muonN_p4;

        Float_t MassErr;
        Float_t vProb;
        Float_t DCA;
        Float_t ppdlPV;
        Float_t ppdlErrPV;
        Float_t cosAlpha;

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
pdgid_(iConfig.getParameter<uint32_t>("onia_pdgid")),
OniaMassCuts_(iConfig.getParameter<std::vector<double>>("onia_mass_cuts")),
isMC_(iConfig.getParameter<bool>("isMC")),
OnlyBest_(iConfig.getParameter<bool>("OnlyBest"))
{
  edm::Service < TFileService > fs;
  onia_tree = fs->make < TTree > ("oniaTree", "Tree of Onia2MuMu");

  onia_tree->Branch("run",   &run,   "run/I");
  onia_tree->Branch("event", &event, "event/I");
  onia_tree->Branch("irank", &irank, "irank/I");

  onia_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
  onia_tree->Branch("muonP_p4",  "TLorentzVector", &muonP_p4);
  onia_tree->Branch("muonN_p4",  "TLorentzVector", &muonN_p4);

  onia_tree->Branch("MassErr",   &MassErr,    "MassErr/F");
  onia_tree->Branch("vProb",     &vProb,      "vProb/F");
  onia_tree->Branch("DCA",       &DCA,        "DCA/F");
  onia_tree->Branch("ppdlPV",    &ppdlPV,     "ppdlPV/F");
  onia_tree->Branch("ppdlErrPV", &ppdlErrPV,  "ppdlErrPV/F");
  onia_tree->Branch("cosAlpha",  &cosAlpha,   "cosAlpha/F");

  onia_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");

  if (isMC_) {
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
   std::cout << "GetAncestor: Inconsistet ancestor, particle does not have a mother " << std::endl;
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

// ------------ method called for each event  ------------
void Onia2MuMuRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  iEvent.getByToken(dimuon_Label,dimuons);

  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  numPrimaryVertices = -1;
  if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();

  run   = iEvent.id().run();
  event = iEvent.id().event();
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

  if (isMC_ && packed.isValid() && pruned.isValid()) {
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
     if ( ! dimuon_pdgId ) std::cout << "Onia2MuMuRootupler: does not found the given decay " << run << "," << event << std::endl;
  }  // end if isMC

  float OniaMassMax_ = OniaMassCuts_[1];
  float OniaMassMin_ = OniaMassCuts_[0];

  if ( dimuons.isValid() && dimuons->size() > 0) {
     for (pat::CompositeCandidateCollection::const_iterator  dimuonCand = dimuons->begin(); dimuonCand!= dimuons->end(); ++dimuonCand) {
        if (dimuonCand->mass() > OniaMassMin_ && dimuonCand->mass() < OniaMassMax_) {
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
           cosAlpha = dimuonCand->userFloat("cosAlpha");
           irank++;
           onia_tree->Fill();
           if (OnlyBest_) break;
        } 
     }
  } else {
     std::cout << "Onia2MuMuRootupler: does not find a valid dimuon combination " << run << "," << event << std::endl;
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
