from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'rootuple-onia2mumu-psi2s-bphskim-v1'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runOnia2MuMu-psi2s.py'
config.JobType.outputFiles = ['rootuple-psi2s.root']

config.Data.inputDataset = '/Charmonium/asanchez-BPHSkim-v1-Run2016B-PromptReco-v2-f3371d8706d2a64e31e991df87f57a46/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 12

config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag  = 'rootuple-onia2mumu-psi2s-bphskim-v1'
config.Site.storageSite = 'TX_YYY_ZZZ'
