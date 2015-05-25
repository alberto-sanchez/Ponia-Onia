// -*- C++ -*-
//
// Package:    OniaRootupler
// Class:      OniaRootupler
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
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include <vector>
#include <sstream>

#include "DataFormats/VertexReco/interface/Vertex.h"

//
// class declaration
//

class OniaRootupler:public edm::EDAnalyzer {
      public:
	explicit OniaRootupler(const edm::ParameterSet &);
	~OniaRootupler();

	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
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
	edm::InputTag dimuon_Label;
	edm::InputTag primaryVertices_Label;
	bool isMC_;

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

        Int_t dimuon_pdgId;
	TLorentzVector gen_dimuon_p4;
	TLorentzVector gen_muonP_p4;
	TLorentzVector gen_muonM_p4;

};

//
// constructors and destructor
//

OniaRootupler::OniaRootupler(const edm::ParameterSet & iConfig):
dimuon_Label(iConfig.getParameter < edm::InputTag > ("dimuons")),
primaryVertices_Label(iConfig.getParameter < edm::InputTag > ("primaryVertices")),
isMC_(iConfig.getParameter < bool > ("isMC")
) {
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
     onia_tree->Branch("dimuon_pdgId",  &dimuon_pdgId,     "dimuon_pdgId/I");
     onia_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
     onia_tree->Branch("gen_muonP_p4",  "TLorentzVector",  &gen_muonP_p4);
     onia_tree->Branch("gen_muonN_p4",  "TLorentzVector",  &gen_muonM_p4);
  }
}

OniaRootupler::~OniaRootupler() {}

//
// member functions
//

//Check recursively if any ancestor of particle is the given one
bool OniaRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}

// ------------ method called for each event  ------------
void OniaRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {
  using namespace edm;
  using namespace std;

  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  iEvent.getByLabel(dimuon_Label,dimuons);

  edm::Handle < std::vector < reco::Vertex > >primaryVertices_handle;
  iEvent.getByLabel(primaryVertices_Label, primaryVertices_handle);

  numPrimaryVertices = -1;
  if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();

  run   = iEvent.id().run();
  event = iEvent.id().event();
  dimuon_pdgId = 0;
  irank = 0;

  // Pruned particles are the one containing "important" stuff
  edm::Handle<edm::View<reco::GenParticle> > pruned;
  iEvent.getByLabel("prunedGenParticles",pruned);

  // Packed particles are all the status 1, so usable to remake jets
  // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
  edm::Handle<edm::View<pat::PackedGenParticle> > packed;
  iEvent.getByLabel("packedGenParticles",packed);

  if (isMC_ && packed.isValid() && pruned.isValid()) {
     dimuon_pdgId  = 0;
     gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
     int foundit   = 0;

     for (size_t i=0; i<pruned->size(); i++) {
         int p_id = abs((*pruned)[i].pdgId());
         const reco::Candidate *aonia = &(*pruned)[i];
         if (( p_id == 443 || p_id == 100443 || p_id == 553 || p_id == 100553 || p_id == 200553) && (aonia->status() == 2)) {
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
               //gen_dimuon_p4.SetPtEtaPhiM(j1->pt(),j1->eta(),j1->phi(),j1->mass());
               gen_dimuon_p4 = gen_muonM_p4 + gen_muonP_p4;   // this should take into account FSR
               break;
            } else {
               foundit = 0;
               dimuon_pdgId = 0;
               gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
            }            
         }  // if ( p_id
     } // for (size

     // sanity check
     if ( ! dimuon_pdgId ) std::cout << "OniaRootupler: does not found the given decay " << run << "," << event << std::endl;
  }  // end if isMC

  bool bestCandidateOnly_ = false;

  if ( dimuons.isValid() && dimuons->size() > 0) {
     for (pat::CompositeCandidateCollection::const_iterator  dimuonCand = dimuons->begin(); dimuonCand!= dimuons->end(); ++dimuonCand) {
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
        if (bestCandidateOnly_) break;
     }
  } else {
     std::cout << "OniaRootupler: does not find a valid dimuon combination " << run << "," << event << std::endl;
  }
  if (!irank  && dimuon_pdgId) {
     std::cout << "OniaRootupler: does not find a reco combination but there is an gen particle " << run << "," << event << std::endl;
     dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
     muonP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
     muonN_p4.SetPtEtaPhiM(0.,0.,0.,0.);
     vProb=-1;
     cosAlpha=-100;
     ppdlPV=-100;
     onia_tree->Fill();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void OniaRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void OniaRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void OniaRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void OniaRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void OniaRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void OniaRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void OniaRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(OniaRootupler);
