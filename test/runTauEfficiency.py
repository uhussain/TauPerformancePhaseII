#cmsRun runTauEfficiency.py inputFiles=/store/relval/CMSSW_8_1_0_pre12/RelValZTT_14TeV/MINIAODSIM/81X_mcRun2_asymptotic_v8_2023D1-v1/00000/4CC5E70D-8B8F-E611-8079-0CC47A4D76D0.root

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

#input cmsRun options
options = VarParsing ('analysis')

#options.inputFiles="/store/relval/CMSSW_9_0_0_pre4/RelValZTT_14TeV/MINIAODSIM/90X_upgrade2023_realistic_v3_2023D4Timing-v1/10000/0CF569B0-B4EC-E611-A970-0025905B85DC.root"
#options.inputFiles="/store/relval/CMSSW_9_0_0_pre4/RelValZTT_14TeV/MINIAODSIM/PU25ns_90X_upgrade2023_realistic_v3_D4TPU200c2-v1/10000/22224E09-53F1-E611-A2E9-0CC47A4D7670.root"
#options.inputFiles="/store/user/adewit/Jun22_MC_91X/GluGluHToTauTau_M125_14TeV_powheg_pythia8/crab_GluGluHToTauTau_M-125-barrel-200PU-miniAOD/170711_074522/0000/out_miniAOD_prod_99.root"
#options.inputFiles="/store/relval/CMSSW_9_3_0_pre4/RelValZTT_14TeV/MINIAODSIM/PU25ns_93X_upgrade2023_realistic_v0_D17PU200-v1/00000/0E0A0B4C-AF89-E711-B67F-A4BF0103375E.root"
#options.outputFile ="Ztt_930_pre4_pu200_test.root"

options.register('inputFileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.parseArguments()

#name the process
process = cms.Process("TreeProducerFromMiniAOD")
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1;
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '93X_upgrade2023_realistic_v0', '')

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
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
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
