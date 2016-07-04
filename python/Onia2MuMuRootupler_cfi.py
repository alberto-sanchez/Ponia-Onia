import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('Onia2MuMuRootupler',
                          dimuons = cms.InputTag("onia2MuMuPAT"),
                          muons = cms.InputTag("replaceme"),
                          primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          onia_pdgid = cms.uint32(443),
                          onia_mass_cuts = cms.vdouble(2.2,4.0),
                          FilterNames = cms.vstring(
                                                    'HLT_Dimuon16_Jpsi',
                                                    'HLT_Dimuon13_PsiPrime',
                                                    'HLT_Dimuon13_Upsilon',
                                                    'HLT_Dimuon10_Jpsi_Barrel',
                                                    'HLT_Dimuon8_PsiPrime_Barrel',
                                                    'HLT_Dimuon8_Upsilon_Barrel',
                                                    'HLT_Dimuon20_Jpsi',
                                                    'HLT_Dimuon0_Phi_Barrel',
                                                    'HLT_HIL1DoubleMu0',
                                                    'HLT_HIL2DoubleMu0',
                                                    'HLT_HIL2Mu3',
                                                    'HLT_HIL3Mu3',
                                                    'HLT_Mu16_TkMu0_dEta18_Onia',
                                                    'HLT_Mu25_TkMu0_dEta18_Onia',
                                                    'HLT_Mu16_TkMu0_dEta18_Phi'
                          ),
                          isMC = cms.bool(True),
                          OnlyBest = cms.bool(True),
                          OnlyGen = cms.bool(False)
                          )
