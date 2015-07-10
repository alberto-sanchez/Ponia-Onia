import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:inputfiles')
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('file:ouputfile'),
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load('Ponia.Onia.Onia2MuMuRootupler_cfi')
process.p = cms.Path(process.rootuple)

process.rootuple.isMC = cms.bool(False)                 # is mc?
process.rootuple.onia_mass_cuts = cms.vdouble(2.,4.)    # you may need to adjust this
