import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('OniaInSingle',
                          dimuons = cms.InputTag("onia2MuMuPAT"),
                          muons = cms.InputTag("replaceme"),
                          primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          onia_pdgid = cms.uint32(443),
                          onia_mass_cuts = cms.vdouble(2.9,3.3),
                          isMC = cms.bool(True),
                          OnlyBest = cms.bool(True),
                          OnlyGen = cms.bool(False),
                          HLTLastFilters = cms.vstring(),
                          propagatorStation1 = cms.PSet(
                                                        useStation2 = cms.bool(False),
                                                        useTrack = cms.string("tracker"),
                                                        useState = cms.string("atVertex"),  # for AOD
                                                        useSimpleGeometry = cms.bool(True),
                          ),
                          propagatorStation2 = cms.PSet(
                                                        useStation2 = cms.bool(True),
                                                        useTrack = cms.string("tracker"),
                                                        useState = cms.string("atVertex"),  # for AOD
                                                        useSimpleGeometry = cms.bool(True),
                          )                         
)
