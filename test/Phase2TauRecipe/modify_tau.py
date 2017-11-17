##########
# Customisation of tau reco&Id for MiniAOD production
# M. Bluj 
# September 2017
##########

#useRun1FixedStrips=False
useRun1FixedStrips=True

print 'Production of MiniAOD+PFTau, with a custom tau setup'

######### Add PFTau ReReco and cloned patTaus with original PFTaus
#Add PFTau modules to patTausTaks 
#A bit complicated logic due to lack of tasks in tauReco sequences 
#add PFTau reco modules to cloned makePatTauTask
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
from PhysicsTools.PatAlgos.tools.helpers import listModules
for module in listModules(process.PFTau):
     process.makePatTausTask.add(module)
# when the sequences migrated to tasks (CMSSW>=940pre1) it will be simple as that (uncomment)
#_makePatTausTaskWithTauReReco
#process.makePatTausTask.add(process.PFTauTask) 

# Modified output to produce and store modified and original PFTaus
# Can be omitted in final prod??
for prod in process.RecoTauTagAOD.outputCommands:
     process.MINIAODSIMoutput.outputCommands.append(prod)

# Clone PAT to use original collections from AOD
# Modified taus (fixed strips) will be stored as slimmedTaus 
# while original ones as slimmedTausOrg
orgProcess='RECO'
process.makePatTausOrgTask = cms.Task()
process.slimmedTausOrg = process.slimmedTaus.clone(
     src = cms.InputTag("selectedPatTausOrg")
)
process.makePatTausOrgTask.add(process.slimmedTausOrg)
process.selectedPatTausOrg = process.selectedPatTaus.clone(
     src = cms.InputTag("patTausOrg")
)
process.makePatTausOrgTask.add(process.selectedPatTausOrg)
process.patTausOrg = process.patTaus.clone(
     genParticleMatch = cms.InputTag("tauMatchOrg"),
     genJetMatch = cms.InputTag("tauGenJetMatchOrg"),
     tauSource = cms.InputTag("hpsPFTauProducer"+"::"+orgProcess),
     tauTransverseImpactParameterSource = cms.InputTag("hpsPFTauTransverseImpactParameters"+"::"+orgProcess),
)
for name in process.patTausOrg.tauIDSources.parameterNames_():
     param = getattr(process.patTausOrg.tauIDSources,name)
     param.setProcessName(orgProcess)
     #print name,param
process.makePatTausOrgTask.add(process.patTausOrg)
process.tauMatchOrg = process.tauMatch.clone(
     src = cms.InputTag("hpsPFTauProducer"+"::"+orgProcess)
)
process.makePatTausOrgTask.add(process.tauMatchOrg)
process.tauGenJetMatchOrg = process.tauGenJetMatch.clone(
     src = cms.InputTag("hpsPFTauProducer"+"::"+orgProcess)
)
process.makePatTausOrgTask.add(process.tauGenJetMatchOrg)
process.makePatTausTask.add(process.makePatTausOrgTask)

#Add original slimmed taus to output
process.MINIAODSIMoutput.outputCommands.append('keep *_slimmedTausOrg_*_*')

######### Setup a custom PFTau (re)reco
# Build Run-1 fixed size strips
if useRun1FixedStrips:
     print '\tUse Run-1 fixed strips' 
     import RecoTauTag.RecoTau.RecoTauPiZeroBuilderPlugins_cfi as piZeroBuilders
     process.ak4PFJetsLegacyHPSPiZeros.builders = cms.VPSet(piZeroBuilders.modStrips)


#End
