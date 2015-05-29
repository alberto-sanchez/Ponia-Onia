input_filename = 'skim*.root'
ouput_filename = 'rootuple-Onia-psi-skim-v2.root'

import FWCore.ParameterSet.Config as cms
import os
import fnmatch
process = cms.Process("Rootuple")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

dir = '/higgs/raid/asanchez/cms/CommonAOD/data_41p1/JPsiMM-v2/'
inputfiles = []
for file in os.listdir(dir):
    if fnmatch.fnmatch(file, input_filename):
        filename = 'file:'+dir+str(file)
        inputfiles.append(filename)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputfiles)
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string(ouput_filename),
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load('Ponia.Onia.OniaRootupler_cfi')
process.p = cms.Path(process.rootuple)
