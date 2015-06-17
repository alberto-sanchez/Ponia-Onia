#Ponia-Onia

The Onia2MuMu Rootupler. This package is mean to be run after the BPH CompactSkim. 

* **Onia2MuMuRootupler** - Rootupler of the Onia2MuMu objects

* Setup: (should run with the same release in CompactSkim, but may run in any other)

```
export SCRAM_ARCH=slc6_amd64_gcc491
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_7_4_4_patch4
cd CMSSW_7_4_4_patch4/src/
cmsenv
git clone https://github.com/alberto-sanchez/Ponia-Onia.git Ponia/Onia
scram b
```

* Run: (use your favorite input sample)

```
vi Ponia/Onia/test/runOnia2MuMuRootupler.py
cmsRun Ponia/Onia/test/runOnia2MuMuRootupler.py
```

