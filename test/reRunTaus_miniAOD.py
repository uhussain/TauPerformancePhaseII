# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --fileout file:SMP-PhaseIITDRSpring17MiniAOD-00011.root --mc --eventcontent MINIAODSIM --runUnscheduled --datatier MINIAODSIM --conditions 91X_upgrade2023_realistic_v3 --step PAT --nThreads 4 --geometry Extended2023D17 --era Phase2C2_timing --python_filename SMP-PhaseIITDRSpring17MiniAOD-00011_1_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n 2880
import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.StandardSequences.Eras import eras

options = VarParsing ('analysis')

#options.inputFiles="/store/relval/CMSSW_9_0_0_pre4/RelValZTT_14TeV/MINIAODSIM/90X_upgrade2023_realistic_v3_2023D4Timing-v1/10000/0CF569B0-B4EC-E611-A970-0025905B85DC.root"
#options.inputFiles="/store/relval/CMSSW_9_1_1_patch1/RelValZTT_14TeV/GEN-SIM-RECO/91X_upgrade2023_realistic_v1_D17-v1/10000/125A37CA-B84A-E711-A9E6-0CC47A7C3428.root"
options.inputFiles="/store/relval/CMSSW_9_1_1_patch3/RelValTTbar_14TeV/GEN-SIM-RECO/91X_upgrade2023_realistic_v3_D17noPUEA300-v2/00000/1AE3DA01-5169-E711-A1E6-0CC47A7C3572.root"
options.outputFile = "MiniAOD_MichalRecipe_91x_TTbar_200_10event.root"

options.register('inputFileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.parseArguments()
process = cms.Process('PAT',eras.Phase2C2_timing)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(10)
            )

# Input source
process.source = cms.Source("PoolSource",
            fileNames = cms.untracked.vstring(options.inputFiles),
                secondaryFileNames = cms.untracked.vstring()
                )

process.options = cms.untracked.PSet(
            allowUnscheduled = cms.untracked.bool(True)
            )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
            annotation = cms.untracked.string('step1 nevts:2880'),
                name = cms.untracked.string('Applications'),
                    version = cms.untracked.string('$Revision: 1.19 $')
                    )

# Output definition

process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
        compressionAlgorithm = cms.untracked.string('LZMA'),
        compressionLevel = cms.untracked.int32(4),
        dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
            ),
        dropMetaData = cms.untracked.string('ALL'),
        eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
        fastCloning = cms.untracked.bool(False),
        fileName = cms.untracked.string(options.outputFile),
        outputCommands = process.MINIAODSIMEventContent.outputCommands,
        overrideInputFileSplitLevels = cms.untracked.bool(True)
         )

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '91X_upgrade2023_realistic_v3', '')

# Path and EndPath definitions
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path()
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path()
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_HBHENoiseFilter = cms.Path()
process.Flag_trkPOG_toomanystripclus53X = cms.Path()
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path()
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.endjob_step,process.MINIAODSIMoutput_step)


process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(4)
process.options.numberOfStreams=cms.untracked.uint32(0)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# customisation of the process.
process.patJets.discriminatorSources = cms.VInputTag(
        cms.InputTag("pfJetBProbabilityBJetTags"),
        cms.InputTag("pfJetProbabilityBJetTags"),
        cms.InputTag("pfTrackCountingHighEffBJetTags"),
        cms.InputTag("pfSimpleSecondaryVertexHighEffBJetTags"),
        cms.InputTag("pfSimpleInclusiveSecondaryVertexHighEffBJetTags"),
        cms.InputTag("pfCombinedSecondaryVertexV2BJetTags"),
        cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
        cms.InputTag("softPFMuonBJetTags"),
        cms.InputTag("softPFElectronBJetTags"),
        cms.InputTag("pfCombinedMVAV2BJetTags"),
        # CTagging
        cms.InputTag('pfCombinedCvsLJetTags'),
        cms.InputTag('pfCombinedCvsBJetTags')
                                       )

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

execfile('modify_tau.py')

# Customisation from command line
#associatePatAlgosToolsTask(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion


#dump_file = open('dump_july6.py','w')
#dump_file.write(process.dumpPython())
