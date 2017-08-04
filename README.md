```
SCRAM_ARCH=slc6_amd64_gcc600
mkdir phase2-tauPerformance
cd phase2-tauPerformance
cmsrel CMSSW_9_1_1_patch3
cd CMSSW_9_1_1_patch3/src
cmsenv
git cms-init
git cms-addpkg RecoTauTag/RecoTau
cd RecoTauTag
git clone https://github.com/uhussain/TauPerformancePhaseII.git
cd ../
scram b -j 8
cd RecoTauTag/phase2Taus/test


```
