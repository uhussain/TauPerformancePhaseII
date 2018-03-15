#cmsRun runTauEfficiency.py inputFiles=/store/relval/CMSSW_8_1_0_pre12/RelValZTT_14TeV/MINIAODSIM/81X_mcRun2_asymptotic_v8_2023D1-v1/00000/4CC5E70D-8B8F-E611-8079-0CC47A4D76D0.root

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

#input cmsRun options
options = VarParsing ('analysis')

#options.inputFiles="/store/relval/CMSSW_9_0_0_pre4/RelValZTT_14TeV/MINIAODSIM/90X_upgrade2023_realistic_v3_2023D4Timing-v1/10000/0CF569B0-B4EC-E611-A970-0025905B85DC.root"
#options.inputFiles="/store/relval/CMSSW_9_0_0_pre4/RelValZTT_14TeV/MINIAODSIM/PU25ns_90X_upgrade2023_realistic_v3_D4TPU200c2-v1/10000/22224E09-53F1-E611-A2E9-0CC47A4D7670.root"
#options.inputFiles="/store/user/adewit/Jun22_MC_91X/GluGluHToTauTau_M125_14TeV_powheg_pythia8/crab_GluGluHToTauTau_M-125-barrel-200PU-miniAOD/170711_074522/0000/out_miniAOD_prod_99.root"
#options.inputFiles="/store/user/agilbert/DYJetsToTauTau_M-50_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/miniaod-prod-020117-0PU/170803_115641/0000/miniaod_1.root"
#options.inputFiles="file:/afs/cern.ch/work/d/dtsiakko/public/miniAOD_TauReco_ak4PFJets_ggH_PU200.root"
#options.inputFiles="file:/afs/cern.ch/work/d/dtsiakko/private/CMSSW_9_4_2_tauAtMiniAOD/src/RecoTauTag/Configuration/test/miniAOD_TauReco_ak4PFJets_ggH_PU0.root"
#options.inputFiles="file:/afs/cern.ch/work/d/dtsiakko/private/CMSSW_9_4_2_tauAtMiniAOD/src/RecoTauTag/Configuration/test/miniAOD_TauReco_ak4PFJets_QCD_PU200.root"
#options.inputFiles="file:/afs/cern.ch/work/d/dtsiakko/private/CMSSW_9_4_2_tauAtMiniAOD/src/RecoTauTag/Configuration/test/miniAOD_TauReco_ak4PFJets_QCD_PU0.root"
#options.inputFiles="file:/afs/cern.ch/work/e/eerodoto/public/miniAOD_TauReco_ak4PFJets_pu200_signal.root"
#options.inputFiles="file:/afs/cern.ch/work/d/dtsiakko/private/phase2-tauPerformance/CMSSW_9_4_2/src/RecoTauTag/TauPerformancePhaseII/test/Plotting/Eleni_store/miniAOD_TauReco_ak4PFJets_QCD_pu.root"
#options.inputFiles="file:/afs/cern.ch/work/d/dtsiakko/private/CMSSW_9_4_2_tauAtMiniAOD/src/RecoTauTag/Configuration/test/miniAOD_TauReco_ak4PFJets_GGTT_PU200.root"
options.inputFiles="file:/afs/cern.ch/work/d/dtsiakko/private/CMSSW_9_4_2_tauAtMiniAOD/src/RecoTauTag/Configuration/test/miniAOD_TauReco_ak4PFJets_ggH_PU0.root"
options.outputFile ="new_store/TwoTaus_ggH_PU0_PUPPI.root"
#options.outputFile ="new_store/TwoTaus_QCD_PU0_PUPPI.root"

options.register('inputFileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.parseArguments()

#name the process
process = cms.Process("TreeProducerFromMiniAOD")
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1;
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '91X_upgrade2023_realistic_v3', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '93X_upgrade2023_realistic_v2', '')

#how many events to run over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


if len(options.inputFileList) > 0 :
    with open(options.inputFileList) as f :
        inputFiles = list((line.strip() for line in f))
else :
    inputFiles = cms.untracked.vstring(options.inputFiles)
    
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
)


##################################################
# Main
process.cutBased = cms.EDAnalyzer("phase2Taus",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#    taus = cms.InputTag("selectedPatTaus"),
    taus = cms.InputTag("slimmedTaus"),
#    tausOrg = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJetsPuppi"),
#   jets = cms.InputTag("slimmedJets"),
    genJets=cms.InputTag("slimmedGenJets"),
    tauID = cms.string("byCombinedIsolationDeltaBetaCorrRaw3Hits"),
    pruned = cms.InputTag("prunedGenParticles")
)

#process.MVA = cms.EDAnalyzer("phase2Taus",
#    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
#    taus = cms.InputTag("slimmedTaus"),
#    tauID = cms.string("byLooseIsolationMVArun2v1DBoldDMwLT"),
#    pruned = cms.InputTag("prunedGenParticles")
#)

###################################################
#Global sequence

process.p = cms.Path(
         process.cutBased
 #        process.MVA
                     )

process.TFileService = cms.Service("TFileService",
        fileName = cms.string(options.outputFile)
)

#print out all processes used when running- useful check to see if module ran
#UNCOMMENT BELOW
#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())
