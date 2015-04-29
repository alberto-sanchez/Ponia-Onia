import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('OniaRootupler',
                          dimuons = cms.InputTag("onia2MuMuPAT"),
                          primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                          isMC = cms.bool(True)
                          )
