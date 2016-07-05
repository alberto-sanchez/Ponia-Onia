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
#include "DataFormats/PatCandidates/interface/Muon.h"
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
        edm::EDGetTokenT<pat::MuonCollection> muon_Label;
        edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
        edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
        int  pdgid_;
        std::vector<double> OniaMassCuts_;
        std::vector<std::string> FilterNames_;
	bool isMC_;
        bool OnlyBest_;
        bool OnlyGen_;

	UInt_t    run;
	ULong64_t event;
        UInt_t    lumiblock;
        UInt_t    nonia;
        UInt_t    nmuons;
        UInt_t    trigger;
        Int_t     charge; 

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

	UInt_t numPrimaryVertices;

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
muon_Label(consumes<pat::MuonCollection>(iConfig.getParameter< edm::InputTag>("muons"))),
primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
pdgid_(iConfig.getParameter<uint32_t>("onia_pdgid")),
OniaMassCuts_(iConfig.getParameter<std::vector<double>>("onia_mass_cuts")),
FilterNames_(iConfig.getParameter<std::vector<std::string>>("FilterNames")),
isMC_(iConfig.getParameter<bool>("isMC")),
OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
OnlyGen_(iConfig.getParameter<bool>("OnlyGen"))
{
  edm::Service < TFileService > fs;
  onia_tree = fs->make < TTree > ("oniaTree", "Tree of Onia2MuMu");

  onia_tree->Branch("run",      &run,      "run/i");
  onia_tree->Branch("event",    &event,    "event/l");
  onia_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

  if (!OnlyGen_) {
    onia_tree->Branch("nonia",    &nonia,    "nonia/i");
    onia_tree->Branch("nmuons",   &nmuons,   "nmuons/i");
    onia_tree->Branch("trigger",  &trigger,  "trigger/i");
    onia_tree->Branch("charge",   &charge,   "charge/I");

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

    onia_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");
  }

  if (isMC_ || OnlyGen_) {
     std::cout << "Onia2MuMuRootupler::Onia2MuMuRootupler: Onia id " << pdgid_ << std::endl;
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

/* Grab Trigger information. Save it in variable trigger, trigger is an uint between 0 and 256, in binary it is:
   (pass 2)(pass 1)(pass 0)
   ex. 7 = pass 0, 1 and 2
   ex. 6 = pass 1, 2
   ex. 1 = pass 0
*/

UInt_t Onia2MuMuRootupler::getTriggerBits(const edm::Event& iEvent ) {
   UInt_t trigger = 0;
   edm::Handle<edm::TriggerResults> triggerResults_handle;
   iEvent.getByToken(triggerResults_Label, triggerResults_handle);
   if (triggerResults_handle.isValid()) {
      const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
      for (unsigned int i = 0; i < FilterNames_.size(); i++) {
         for (int version = 1; version < 9; version++) {
            std::stringstream ss;
            ss << FilterNames_[i] << "_v" << version;
            unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label().c_str());
            if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
               trigger += (1<<i);
               break;
            }
         }
      }
   } else std::cout << "Onia2MuMuRootupler::getTriggerBits: *** NO triggerResults found *** " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;
   return trigger;
}

// ------------ method called for each event  ------------
void Onia2MuMuRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muon_Label,muons);

  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  iEvent.getByToken(dimuon_Label,dimuons);

  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock(); 

  numPrimaryVertices = 0;
  if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();
  trigger = getTriggerBits(iEvent);

  dimuon_pdgId = 0;
  mother_pdgId = 0;
  nonia  = 0;
  nmuons = 0;

  dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  muonP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  muonN_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muonP_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muonM_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  // Pruned particles are the one containing "important" stuff
  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genCands_, pruned);

  // Packed particles are all the status 1. The navigation to pruned is possible (the other direction should be made by hand)
  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packCands_,  packed);

  if ( (isMC_ || OnlyGen_) && packed.isValid() && pruned.isValid() ) {
    for (size_t i=0; i<pruned->size(); i++) {
      const reco::Candidate *aonia = &(*pruned)[i];
      if ( (abs(aonia->pdgId()) == pdgid_) && (aonia->status() == 2) ) {
        int foundit = 1;
        dimuon_pdgId = aonia->pdgId();
        for ( size_t j=0; j<packed->size(); j++ ) { //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
          const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
          const reco::Candidate * d = &(*packed)[j];
          if ( motherInPrunedCollection != nullptr && (d->pdgId() ==  13 ) && isAncestor(aonia , motherInPrunedCollection) ) {
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
        } else dimuon_pdgId = 0;            
      }  // if ( p_id
    } // for (size
    if ( ! dimuon_pdgId ) std::cout << "Onia2MuMuRootupler: does not found the given decay " << run << "," << event << std::endl; // sanity check
  }  // end if isMC

  float OniaMassMax_ = OniaMassCuts_[1];
  float OniaMassMin_ = OniaMassCuts_[0];

  bool already_stored = false;
  if ( ! OnlyGen_ ) { // we will look for dimuons, then for muons
    if ( dimuons.isValid() && dimuons->size() > 0) {
      for ( pat::CompositeCandidateCollection::const_iterator dimuonCand = dimuons->begin(); dimuonCand != dimuons->end(); ++dimuonCand ) {
        if (dimuonCand->mass() > OniaMassMin_ && dimuonCand->mass() < OniaMassMax_ && dimuonCand->charge() == 0) {
          dimuon_p4.SetPtEtaPhiM(dimuonCand->pt(),dimuonCand->eta(),dimuonCand->phi(),dimuonCand->mass());
          reco::Candidate::LorentzVector vP = dimuonCand->daughter("muon1")->p4();
          reco::Candidate::LorentzVector vM = dimuonCand->daughter("muon2")->p4();
          if ( dimuonCand->daughter("muon1")->charge() < 0 ) {
              vP = dimuonCand->daughter("muon2")->p4();
              vM = dimuonCand->daughter("muon1")->p4();
          }
          muonP_p4.SetPtEtaPhiM(vP.pt(),vP.eta(),vP.phi(),vP.mass());
          muonN_p4.SetPtEtaPhiM(vM.pt(),vM.eta(),vM.phi(),vM.mass());
          MassErr = dimuonCand->userFloat("MassErr");
          vProb = dimuonCand->userFloat("vProb");
          DCA = -1.;
          if (dimuonCand->hasUserFloat("DCA"))  DCA = dimuonCand->userFloat("DCA");
          ppdlPV = dimuonCand->userFloat("ppdlPV");
          ppdlErrPV = dimuonCand->userFloat("ppdlErrPV");
          ppdlBS = dimuonCand->userFloat("ppdlBS");
          ppdlErrBS = dimuonCand->userFloat("ppdlErrBS");
          cosAlpha = dimuonCand->userFloat("cosAlpha");
          charge = dimuonCand->charge(); 
          TVector3 pperp(dimuonCand->px(),dimuonCand->py(),0);
          lxyPV = ppdlPV * pperp.Perp() / dimuonCand->mass();
          lxyBS = ppdlBS * pperp.Perp() / dimuonCand->mass();
          nonia++;
          if (OnlyBest_) break;
          else { 
            onia_tree->Fill();   // be aware, we are storing all combinations
            already_stored = true;
          }
        } 
      }
    } // if ( dimuons.isValid()
    if ( nonia == 0 && muons.isValid() && muons->size() > 0 ) {
        int mcharge1 = 0, mcharge2 = 0;
        reco::Candidate::LorentzVector v1, v2;
        for ( pat::MuonCollection::const_iterator muonCand = muons->begin(); muonCand!= muons->end(); ++muonCand ) {
          nmuons++;
          if (nmuons == 1) { 
            mcharge1 = muonCand->charge();
            v1 = muonCand->p4();
          } else {
            if ( mcharge1*muonCand->charge() < 0  && mcharge2 == 0 ) { 
              mcharge2 = muonCand->charge();
              v2 = muonCand->p4();
              nmuons = 2;
              break;    // we store only 2 muons
            } 
          }
        }
        if ( mcharge1 > 0 ) { 
          muonP_p4.SetPtEtaPhiM(v1.pt(),v1.eta(),v1.phi(),v1.mass());
          if (mcharge2 < 0 ) muonN_p4.SetPtEtaPhiM(v2.pt(),v2.eta(),v2.phi(),v2.mass());
        } else {
          muonN_p4.SetPtEtaPhiM(v1.pt(),v1.eta(),v1.phi(),v1.mass());
          if (mcharge2 > 0 ) muonP_p4.SetPtEtaPhiM(v2.pt(),v2.eta(),v2.phi(),v2.mass());
        }
    }
  }  // !OnlyGen_

  if ( !already_stored ) {  // we have to make sure, we are not double storing an combination
    if ( !isMC_ ) {
      if ( nonia > 0 ) onia_tree->Fill();   // if not MC filter out
    } else onia_tree->Fill();
  }
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
