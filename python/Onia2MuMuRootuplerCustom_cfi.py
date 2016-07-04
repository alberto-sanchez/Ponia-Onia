import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('Onia2MuMuRootuplerCustom',
                          dimuons = cms.InputTag("onia2MuMuPAT"),
                          muons = cms.InputTag("replaceme"),
                          primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          onia_pdgid = cms.uint32(443),
                          onia_mass_cuts = cms.vdouble(2.9,3.3),
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
                          SingleFilterNames = cms.vstring(
                                                    'HLT_Mu20',
                                                    'HLT_Mu24_eta2p1',
                                                    'HLT_Mu27',
                                                    'HLT_Mu45_eta2p1',
                                                    'HLT_Mu50',
                                                    'HLT_Mu55',
                                                    'HLT_Mu8',
                                                    'HLT_Mu17',
                                                    'HLT_Mu24',
                                                    'HLT_Mu34',
                                                    'HLT_Mu7p5_Track2_Jpsi'
                          ),
                          isMC = cms.bool(True),
                          OnlyBest = cms.bool(True),
                          OnlyGen = cms.bool(False),
                          HLTLastFilters = cms.vstring('hltDisplacedmumuFilterDimuon10JpsiBarrel',
                                                       'hltDisplacedmumuFilterDimuon16Jpsi',
                                                       'hltDisplacedmumuFilterDimuon20Jpsi'
                          ),
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
