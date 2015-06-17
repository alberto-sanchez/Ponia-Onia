import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('Onia2MuMuRootupler',
                          dimuons = cms.InputTag("onia2MuMuPAT"),
                          primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                          isMC = cms.bool(True)
                          )
