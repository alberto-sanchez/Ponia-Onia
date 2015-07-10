import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('Onia2MuMuRootupler',
                          dimuons = cms.InputTag("onia2MuMuPAT"),
                          primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          onia_pdgid = cms.uint32(443),
                          onia_mass_cuts = cms.vdouble(2.2,4.0),
                          isMC = cms.bool(True),
                          OnlyBest = cms.bool(True),
                          OnlyGen = cms.bool(False)
                          )
