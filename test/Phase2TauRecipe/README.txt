#
# M. Bluj, NCBj Warsaw
# Created: 20 Oct. 2017
# 

# Installation recipe (basically as in the CMS object school)
scram project -n CMSSW_9_2_9_patch1_tausAtMiniAOD CMSSW CMSSW_9_2_9_patch1
cd CMSSW_9_2_9_patch1_tausAtMiniAOD/src/
git cms-addpkg DataFormats/TauReco DataFormats/PatCandidates 
PhysicsTools/PatAlgos RecoTauTag/Configuration RecoTauTag/RecoTau
git remote add modifiedRepo https://github.com/steggema/cmssw.git
git fetch modifiedRepo
git checkout modifiedRepo/CMSSW_9_2_X_TauRecoMiniAOD
scram b -r -j 4

# Configuration files for tau rereco at miniAOD (in this dir):
tau_miniaod_2023_tst.py # mail cfg file
modify_tau_miniaod.py # file with customizations, currently it switches dynamic strips to fix size ones 

# to be run as follows
cmsRun tau_miniaod_2023_tst.py


