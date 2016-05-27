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
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_1.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_10.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_11.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_12.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_13.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_14.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_15.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_16.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_17.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_18.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_19.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_2.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_20.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_21.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_22.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_23.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_24.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_3.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_4.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_5.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_6.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_7.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_8.root',
       '/store/group/phys_bphys/asanchez/Charmonium/BPHSkim-v1-Run2016B-PromptReco-v2/160521_011201/0000/BPHSkim_9.root'
)
)

process.TFileService = cms.Service("TFileService",fileName = cms.string('rootuple-psi2s.root'))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))
process.load('Ponia.Onia.Onia2MuMuRootuplerCustom_cfi')
process.p = cms.Path(process.rootuple)
process.rootuple.onia_mass_cuts = cms.vdouble(3.35,4.0)
process.rootuple.isMC = cms.bool(False)

process.rootuple.HLTLastFilters = cms.vstring(
        'hltDisplacedmumuFilterDimuon8PsiPrimeBarrel',
        'hltDisplacedmumuFilterDimuon13PsiPrime'
)
