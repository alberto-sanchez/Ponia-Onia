#Ponia-Onia

The Onia2MuMu Rootupler. This package is mean to be run after the BPH CompactSkim. 

* **Onia2MuMuRootupler** - Rootupler of the Onia2MuMu objects
* **OniaMM**             - Rootupler for genparticles objects, AOD

* Setup: (should run with the same release in CompactSkim, but may run in any other)

```
export SCRAM_ARCH=slc6_amd64_gcc491
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_7_4_6_patch6
cd CMSSW_7_4_6_patch6/src/
cmsenv
git clone https://github.com/alberto-sanchez/Ponia-Onia.git Ponia/Onia
scram b
```

* Run: (use your favorite input sample)

```
vi Ponia/Onia/test/runOnia2MuMuRootupler.py
cmsRun Ponia/Onia/test/runOnia2MuMuRootupler.py
```

#CompactSkim
* Setup: (it is part of CMSSW_7_5_X onwards, but for now if you want to use it in CMSSW_7_4_X, you can do)

```
export SCRAM_ARCH=slc6_amd64_gcc491
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_7_4_6_patch6
cd CMSSW_7_4_6_patch6/src/
cmsenv
git cms-merge-topic alberto-sanchez:onia2mumu-74x
scram b
``` 

* Run: (adjust the inputs)

```
vi HeavyFlavorAnalysis/Skimming/test/runCompactSkim.py
cmsRun HeavyFlavorAnalysis/Skimming/test/runCompactSkim.py
```

