```
SCRAM_ARCH=slc6_amd64_gcc600
mkdir phase2-tauPerformance
cd phase2-tauPerformance
cmsrel CMSSW_9_4_2
cd CMSSW_9_4_2/src
cmsenv
git cms-init
git cms-addpkg RecoTauTag/RecoTau
cd RecoTauTag
git clone -b TwoTauCollections git@github.com:uhussain/TauPerformancePhaseII.git
cd ../
scram b -j 8
cd RecoTauTag/phase2Taus/test
cmsRun runTauEfficiency.py

```
