from FWCore.ParameterSet.VarParsing import VarParsing
opt = VarParsing ('analysis')
opt.parseArguments()

from PhysicsTools.PatAlgos.patTemplate_cfg import *
process.options.allowUnscheduled = cms.untracked.bool(True) # switch to uncheduled mode

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.source.fileNames    = cms.untracked.vstring(opt.inputFiles)
process.maxEvents.input     = -1
process.out.fileName        = opt.outputFile
process.options.wantSummary = True

# make patCandidates, select and clean them
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
process.load("PhysicsTools.PatAlgos.cleaningLayer1.cleanPatCandidates_cff")

process.patMuons.embedTrack  = True
process.selectedPatMuons.cut = cms.string('muonID(\"TMOneStationTight\")'
                    '&& abs(innerTrack.dxy) < 3.0'
                    '&& abs(innerTrack.dz)  < 30'
                    '&& track.hitPattern.trackerLayersWithMeasurement > 5'
                    '&& innerTrack.hitPattern.pixelLayersWithMeasurement > 1'
                    '&& innerTrack.normalizedChi2 < 1.8'
                    )

#make patTracks
from PhysicsTools.PatAlgos.tools.trackTools import makeTrackCandidates
makeTrackCandidates(process,
                        label        = 'TrackCands',                  # output collection
                        tracks       = cms.InputTag('generalTracks'), # input track collection
                        particleType = 'pi+',                         # particle type (for assigning a mass)
                        preselection = 'pt > 0.7',                    # preselection cut on candidates
                        selection    = 'pt > 0.7',                    # selection on PAT Layer 1 objects
                        isolation    = {},                            # isolations to use (set to {} for None)
                        isoDeposits  = [],
                        mcAs         = None                           # replicate MC match as the one used for Muons
                        )

process.patTrackCands.embedTrack = True

# dimuon = Onia2MUMU
process.load('HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi')
process.onia2MuMuPAT.muons=cms.InputTag("cleanPatMuons")
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag("offlinePrimaryVertices")
process.onia2MuMuPAT.beamSpotTag=cms.InputTag("offlineBeamSpot")

process.onia2MuMuPATCounter = cms.EDFilter('CandViewCountFilter',
    src = cms.InputTag('onia2MuMuPAT'),
    minNumber = cms.uint32(1),
    filter = cms.bool(True)
    )

# reduce MC genParticles a la miniAOD
process.load("PhysicsTools.PatAlgos.slimming.genParticles_cff")
process.packedGenParticles.inputVertices = cms.InputTag("offlinePrimaryVertices")

# make photon candidate conversions for P-wave studies
process.load("HeavyFlavorAnalysis.PhotonConversion.PhotonConversionProducer_cfi")

# Pick branches you want to keep
SlimmedEventContent=[
                     'keep *_offlinePrimaryVertices_*_*',
                     'keep *_offlineBeamSpot_*_*',
#
                     'keep *_TriggerResults_*_HLT',
                     'keep *_gtDigis_*_RECO',
#
                     'keep *_cleanPatTrackCands_*_*',
                     'keep *_PhotonCandidates_*_*',
                     'keep *_onia2MuMuPAT_*_*',
#
                     'keep *_generalV0Candidates_*_*',
#
                     'keep patPackedGenParticles_packedGenParticles_*_*',
                     'keep recoGenParticles_prunedGenParticles_*_*',
                     'keep PileupSummaryInfos_*_*_*',
                     'keep GenFilterInfo_*_*_*',
                     'keep GenEventInfoProduct_generator_*_*',
                     'keep GenRunInfoProduct_*_*_*'
]

process.out.outputCommands = cms.untracked.vstring('drop *', *SlimmedEventContent )
process.FilterOutput = cms.Path(process.onia2MuMuPATCounter)
process.out.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('FilterOutput'))

