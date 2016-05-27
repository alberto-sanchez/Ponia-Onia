import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'/store/group/phys_bphys/asanchez/MuOnia/BPHSkim-v1-Run2015C_25ns-16Dec2015-v1/160519_173955/0000/BPHSkim_1.root',
'/store/group/phys_bphys/asanchez/MuOnia/BPHSkim-v1-Run2015C_25ns-16Dec2015-v1/160519_173955/0000/BPHSkim_2.root',
'/store/group/phys_bphys/asanchez/MuOnia/BPHSkim-v1-Run2015C_25ns-16Dec2015-v1/160519_173955/0000/BPHSkim_3.root'
)
)

process.TFileService = cms.Service("TFileService",fileName = cms.string('rootuple-upsns.root'))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))
process.load('Ponia.Onia.Onia2MuMuRootuplerCustom_cfi')
process.p = cms.Path(process.rootuple)
process.rootuple.onia_mass_cuts = cms.vdouble(8.,12.)
process.rootuple.isMC = cms.bool(False)
process.rootuple.HLTLastFilters = cms.vstring(
        'hltDisplacedmumuFilterDimuon8UpsilonBarrel',
        'hltDisplacedmumuFilterDimuon13Upsilon'
)
