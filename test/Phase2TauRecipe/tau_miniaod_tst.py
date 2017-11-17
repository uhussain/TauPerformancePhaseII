import FWCore.ParameterSet.Config as cms

#from filenames_aod import filenames_aod

runSignal=True
#runSignal=False
#maxEvents=100
maxEvents=-1

#jetCollection = 'slimmedJets'
jetCollection = 'patAK4PFJets'

print 'Running Tau reco\&id with MiniAOD inputs:'
print '\t runSignal:', runSignal
print '\t Jet collection:', jetCollection
 

def addTauReReco(process):
	process.load('PhysicsTools.PatAlgos.producersLayer1.tauProducer_cff')
	process.load('PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff')
	process.load('PhysicsTools.PatAlgos.slimming.slimmedTaus_cfi')

	process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
	# process.ptau = cms.Path(process.PFTau)
	process.PATTauSequence = cms.Sequence(process.PFTau+process.makePatTaus+process.selectedPatTaus)#+process.slimmedTaus)
	process.TauReco = cms.Path(process.PATTauSequence)
	# process.task.add(process.TauReco)
	process.selectedPatTaus.cut="pt > 18. && tauID(\'decayModeFindingNewDMs\')> 0.5"

process = cms.Process("TAURECO")
process.load("Configuration.StandardSequences.MagneticField_cff") # for CH reco
process.load("Configuration.Geometry.GeometryRecoDB_cff")

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(maxEvents)
)
print '\t maxEvents:', process.maxEvents.input.value()

if runSignal:
	readFiles.extend( [
			'file:/home/mbluj/work/data/MiniAOD/92X/GluGluHToTauTau_M125_13TeV_powheg_pythia8_FlatPU28to62HcalNZSRAW_PhaseIFall16_90X_upgrade2017_realistic_v6_C1-v1/06290054-DD1A-E711-9566-C45444922991.root',
			] )
else:
	readFiles.extend( [
			'file:/home/mbluj/work/data/MiniAOD/92X/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_skimByJetV2/patTuple_1.root',
			] )

addTauReReco(process)

process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string('.'),
        connectionRetrialPeriod = cms.untracked.int32(10),
        connectionRetrialTimeOut = cms.untracked.int32(60),
        connectionTimeOut = cms.untracked.int32(0),
        enableConnectionSharing = cms.untracked.bool(True),
        enablePoolAutomaticCleanUp = cms.untracked.bool(False),
        enableReadOnlySessionOnUpdateConnection = cms.untracked.bool(False),
        idleConnectionCleanupPeriod = cms.untracked.int32(10),
        messageLevel = cms.untracked.int32(0)
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    #globaltag = cms.string('90X_upgrade2017_TSG_Hcal_V2'),
    globaltag = cms.string('92X_upgrade2017_realistic_v10'),
    pfnPostfix = cms.untracked.string('None'),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet()
)

process.load('Configuration.EventContent.EventContent_cff')
process.output = cms.OutputModule(
	'PoolOutputModule',
	fileName=cms.untracked.string('outputFULL.root'),
	fastCloning=cms.untracked.bool(False),
	dataset=cms.untracked.PSet(
		dataTier=cms.untracked.string('MINIAODSIM'),
		filterName=cms.untracked.string('')
	),
	#outputCommands=cms.untracked.vstring('keep *'),
	outputCommands = process.MINIAODSIMEventContent.outputCommands,
	SelectEvents=cms.untracked.PSet(
		SelectEvents=cms.vstring('*',)
	)
)
for prod in process.RecoTauTagAOD.outputCommands:
     process.output.outputCommands.append(prod)
process.output.outputCommands.append('keep *_selectedPatTaus_*_*')

if runSignal:
	process.output.fileName='miniAOD_TauReco_ggH.root'
	if jetCollection == 'patAK4PFJets':
		process.output.fileName='miniAOD_TauReco_ak4PFJets_ggH.root'
else:
	process.output.fileName='miniAOD_TauReco_QCD.root'
	if jetCollection == 'patAK4PFJets':
		process.output.fileName='miniAOD_TauReco_ak4PFJets_QCD.root'
process.out = cms.EndPath(process.output)

# TRYING TO MAKE THINGS MINIAOD COMPATIBLE, FROM THE START, TO THE END, 1 BY 1
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.ak4PFJetsPAT = ak4PFJets.clone(
	src=cms.InputTag("packedPFCandidates")
	)
from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import _patJets
process.patAK4PFJets = _patJets.clone(
	jetSource            = cms.InputTag("ak4PFJetsPAT"),
	addJetCorrFactors    = cms.bool(False),
	jetCorrFactorsSource = cms.VInputTag(),
	addBTagInfo          = cms.bool(False),
	addDiscriminators    = cms.bool(False),
	discriminatorSources = cms.VInputTag(),
	addAssociatedTracks  = cms.bool(False),
	addJetCharge         = cms.bool(False),
	addGenPartonMatch    = cms.bool(False),
	embedGenPartonMatch  = cms.bool(False),
	addGenJetMatch       = cms.bool(False),
	getJetMCFlavour      = cms.bool(False),
	addJetFlavourInfo    = cms.bool(False),
	)
 
# so this adds all tracks to jet in some deltaR region. we don't have tracks so don't need it :D
# process.ak4PFJetTracksAssociatorAtVertex.jets = cms.InputTag(jetCollection)

# Remove ak4PFJetTracksAssociatorAtVertex from recoTauCommonSequence
# Remove pfRecoTauTagInfoProducer from recoTauCommonSequence since it uses the jet-track association
# HOWEVER, may use https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#Isolated_Tracks
# probably needs recovery of the two modules above


process.recoTauAK4PatJets08Region = cms.EDProducer("RecoTauPatJetRegionProducer",
    deltaR = cms.double(0.8),
    maxJetAbsEta = cms.double(2.5),
    minJetPt = cms.double(14.0),
    pfCandAssocMapSrc = cms.InputTag(""),
    pfCandSrc = cms.InputTag("packedPFCandidates"),
    src = cms.InputTag(jetCollection)
)

process.recoTauPileUpVertices.src = cms.InputTag("offlineSlimmedPrimaryVertices")

process.ak4PFJetsLegacyHPSPiZeros.builders[0].qualityCuts.primaryVertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices")
process.ak4PFJetsLegacyHPSPiZeros.jetSrc = cms.InputTag(jetCollection)

for builder in process.ak4PFJetsRecoTauChargedHadrons.builders:
    builder.qualityCuts.primaryVertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices")
process.ak4PFJetsRecoTauChargedHadrons.jetSrc = cms.InputTag(jetCollection)

# FIXME - remove builder from tracks. well, because there are no tracks in miniAOD
# but see comment above!
process.ak4PFJetsRecoTauChargedHadrons.builders =  cms.VPSet(process.ak4PFJetsRecoTauChargedHadrons.builders[0], process.ak4PFJetsRecoTauChargedHadrons.builders[2])

process.recoTauCommonSequence = cms.Sequence(process.recoTauAK4PatJets08Region + process.recoTauPileUpVertices)

process.combinatoricRecoTaus = cms.EDProducer("RecoPFBaseTauProducer",
    buildNullTaus = cms.bool(False),
    builders = cms.VPSet(cms.PSet(
        decayModes = cms.VPSet(cms.PSet(
            maxPiZeros = cms.uint32(0),
            maxTracks = cms.uint32(6),
            nCharged = cms.uint32(1),
            nPiZeros = cms.uint32(0)
        ), 
            cms.PSet(
                maxPiZeros = cms.uint32(6),
                maxTracks = cms.uint32(6),
                nCharged = cms.uint32(1),
                nPiZeros = cms.uint32(1)
            ), 
            cms.PSet(
                maxPiZeros = cms.uint32(5),
                maxTracks = cms.uint32(6),
                nCharged = cms.uint32(1),
                nPiZeros = cms.uint32(2)
            ), 
            cms.PSet(
                maxPiZeros = cms.uint32(0),
                maxTracks = cms.uint32(6),
                nCharged = cms.uint32(2),
                nPiZeros = cms.uint32(0)
            ), 
            cms.PSet(
                maxPiZeros = cms.uint32(3),
                maxTracks = cms.uint32(6),
                nCharged = cms.uint32(2),
                nPiZeros = cms.uint32(1)
            ), 
            cms.PSet(
                maxPiZeros = cms.uint32(0),
                maxTracks = cms.uint32(6),
                nCharged = cms.uint32(3),
                nPiZeros = cms.uint32(0)
            ), 
            cms.PSet(
                maxPiZeros = cms.uint32(3),
                maxTracks = cms.uint32(6),
                nCharged = cms.uint32(3),
                nPiZeros = cms.uint32(1)
            )),
        isolationConeSize = cms.double(0.5),
        minAbsPhotonSumPt_insideSignalCone = cms.double(2.5),
        minAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
        minRelPhotonSumPt_insideSignalCone = cms.double(0.1),
        minRelPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
        name = cms.string('combinatoric'),
        pfCandSrc = cms.InputTag("packedPFCandidates"),
        plugin = cms.string('RecoBaseTauBuilderCombinatoricPlugin'),
        qualityCuts = cms.PSet(
            isolationQualityCuts = cms.PSet(
                maxDeltaZ = cms.double(0.2),
                maxTrackChi2 = cms.double(100.0),
                maxTransverseImpactParameter = cms.double(0.03),
                minGammaEt = cms.double(1.5),
                minTrackHits = cms.uint32(8),
                minTrackPixelHits = cms.uint32(0),
                minTrackPt = cms.double(1.0),
                minTrackVertexWeight = cms.double(-1.0)
            ),
            leadingTrkOrPFCandOption = cms.string('leadPFCand'),
            primaryVertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
            pvFindingAlgo = cms.string('closestInDeltaZ'),
            recoverLeadingTrk = cms.bool(False),
            signalQualityCuts = cms.PSet(
                maxDeltaZ = cms.double(0.4),
                maxTrackChi2 = cms.double(100.0),
                maxTransverseImpactParameter = cms.double(0.1),
                minGammaEt = cms.double(0.5),
                minNeutralHadronEt = cms.double(30.0),
                minTrackHits = cms.uint32(3),
                minTrackPixelHits = cms.uint32(0),
                minTrackPt = cms.double(0.5),
                minTrackVertexWeight = cms.double(-1.0)
            ),
            vertexTrackFiltering = cms.bool(False),
            vxAssocQualityCuts = cms.PSet(
                maxTrackChi2 = cms.double(100.0),
                maxTransverseImpactParameter = cms.double(0.1),
                minGammaEt = cms.double(0.5),
                minTrackHits = cms.uint32(3),
                minTrackPixelHits = cms.uint32(0),
                minTrackPt = cms.double(0.5),
                minTrackVertexWeight = cms.double(-1.0)
            )
        ),
        signalConeSize = cms.string('max(min(0.1, 3.0/pt()), 0.05)')
    )),
    chargedHadronSrc = cms.InputTag("ak4PFJetsRecoTauChargedHadrons"),
    jetRegionSrc = cms.InputTag("recoTauAK4PatJets08Region"),
    jetSrc = cms.InputTag(jetCollection),
    maxJetAbsEta = cms.double(2.5),
    minJetPt = cms.double(14.0),
    modifiers = cms.VPSet(cms.PSet(
        name = cms.string('sipt'),
        plugin = cms.string('RecoBaseTauImpactParameterSignificancePlugin'),
        qualityCuts = cms.PSet(
            isolationQualityCuts = cms.PSet(
                maxDeltaZ = cms.double(0.2),
                maxTrackChi2 = cms.double(100.0),
                maxTransverseImpactParameter = cms.double(0.03),
                minGammaEt = cms.double(1.5),
                minTrackHits = cms.uint32(8),
                minTrackPixelHits = cms.uint32(0),
                minTrackPt = cms.double(1.0),
                minTrackVertexWeight = cms.double(-1.0)
            ),
            leadingTrkOrPFCandOption = cms.string('leadPFCand'),
            primaryVertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
            pvFindingAlgo = cms.string('closestInDeltaZ'),
            recoverLeadingTrk = cms.bool(False),
            signalQualityCuts = cms.PSet(
                maxDeltaZ = cms.double(0.4),
                maxTrackChi2 = cms.double(100.0),
                maxTransverseImpactParameter = cms.double(0.1),
                minGammaEt = cms.double(0.5),
                minNeutralHadronEt = cms.double(30.0),
                minTrackHits = cms.uint32(3),
                minTrackPixelHits = cms.uint32(0),
                minTrackPt = cms.double(0.5),
                minTrackVertexWeight = cms.double(-1.0)
            ),
            vertexTrackFiltering = cms.bool(False),
            vxAssocQualityCuts = cms.PSet(
                maxTrackChi2 = cms.double(100.0),
                maxTransverseImpactParameter = cms.double(0.1),
                minGammaEt = cms.double(0.5),
                minTrackHits = cms.uint32(3),
                minTrackPixelHits = cms.uint32(0),
                minTrackPt = cms.double(0.5),
                minTrackVertexWeight = cms.double(-1.0)
            )
        )
    ), 
        # cms.PSet(
        #     DataType = cms.string('AOD'),
        #     EcalStripSumE_deltaEta = cms.double(0.03),
        #     EcalStripSumE_deltaPhiOverQ_maxValue = cms.double(0.5),
        #     EcalStripSumE_deltaPhiOverQ_minValue = cms.double(-0.1),
        #     EcalStripSumE_minClusEnergy = cms.double(0.1),
        #     ElecPreIDLeadTkMatch_maxDR = cms.double(0.01),
        #     ElectronPreIDProducer = cms.InputTag("elecpreid"),
        #     maximumForElectrionPreIDOutput = cms.double(-0.1),
        #     name = cms.string('elec_rej'),
        #     plugin = cms.string('RecoTauElectronRejectionPlugin')
        # ), 
        cms.PSet(
            dRaddNeutralHadron = cms.double(0.12),
            dRaddPhoton = cms.double(-1.0),
            minGammaEt = cms.double(10.0),
            minNeutralHadronEt = cms.double(50.0),
            name = cms.string('tau_en_reconstruction'),
            plugin = cms.string('PFRecoBaseTauEnergyAlgorithmPlugin'),
            verbosity = cms.int32(0)
        ), 
        cms.PSet(
            name = cms.string('tau_mass'),
            plugin = cms.string('PFRecoBaseTauMassPlugin'),
            verbosity = cms.int32(0)
        ), 
        # cms.PSet(
        #     name = cms.string('TTIworkaround'),
        #     pfTauTagInfoSrc = cms.InputTag("pfRecoTauTagInfoProducer"),
        #     plugin = cms.string('RecoTauTagInfoWorkaroundModifer')
        # )
        ),
    outputSelection = cms.string('leadPFChargedHadrCand().isNonnull()'),
    pfCandSrc = cms.InputTag("packedPFCandidates"),
    piZeroSrc = cms.InputTag("ak4PFJetsLegacyHPSPiZeros")
)

process.hpsSelectionDiscriminator = cms.EDProducer("PFRecoBaseTauDiscriminationByHPSSelection",
    PFBaseTauProducer = cms.InputTag("combinatoricRecoTaus"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and')
    ),
    decayModes = cms.VPSet(cms.PSet(
        applyBendCorrection = cms.PSet(
            eta = cms.bool(True),
            mass = cms.bool(True),
            phi = cms.bool(True)
        ),
        maxMass = cms.string('1.'),
        minMass = cms.double(-1000.0),
        nCharged = cms.uint32(1),
        nChargedPFCandsMin = cms.uint32(1),
        nPiZeros = cms.uint32(0),
        nTracksMin = cms.uint32(1)
    ), 
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(True),
                mass = cms.bool(True),
                phi = cms.bool(True)
            ),
            assumeStripMass = cms.double(0.1349),
            maxMass = cms.string('max(1.3, min(1.3*sqrt(pt/100.), 4.2))'),
            minMass = cms.double(0.3),
            nCharged = cms.uint32(1),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(1),
            nTracksMin = cms.uint32(1)
        ), 
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(True),
                mass = cms.bool(True),
                phi = cms.bool(True)
            ),
            assumeStripMass = cms.double(0.0),
            maxMass = cms.string('max(1.2, min(1.2*sqrt(pt/100.), 4.0))'),
            maxPi0Mass = cms.double(0.2),
            minMass = cms.double(0.4),
            minPi0Mass = cms.double(0.05),
            nCharged = cms.uint32(1),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(2),
            nTracksMin = cms.uint32(1)
        ), 
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(False),
                mass = cms.bool(False),
                phi = cms.bool(False)
            ),
            maxMass = cms.string('1.2'),
            minMass = cms.double(0.0),
            nCharged = cms.uint32(2),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(0),
            nTracksMin = cms.uint32(2)
        ), 
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(False),
                mass = cms.bool(False),
                phi = cms.bool(False)
            ),
            maxMass = cms.string('max(1.2, min(1.2*sqrt(pt/100.), 4.0))'),
            minMass = cms.double(0.0),
            nCharged = cms.uint32(2),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(1),
            nTracksMin = cms.uint32(2)
        ), 
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(False),
                mass = cms.bool(False),
                phi = cms.bool(False)
            ),
            maxMass = cms.string('1.5'),
            minMass = cms.double(0.8),
            nCharged = cms.uint32(3),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(0),
            nTracksMin = cms.uint32(2)
        ), 
        cms.PSet(
            applyBendCorrection = cms.PSet(
                eta = cms.bool(False),
                mass = cms.bool(False),
                phi = cms.bool(False)
            ),
            maxMass = cms.string('1.6'),
            minMass = cms.double(0.9),
            nCharged = cms.uint32(3),
            nChargedPFCandsMin = cms.uint32(1),
            nPiZeros = cms.uint32(1),
            nTracksMin = cms.uint32(2)
        )),
    matchingCone = cms.double(0.5),
    minPixelHits = cms.int32(1),
    minTauPt = cms.double(0.0),
    requireTauChargedHadronsToBeChargedPFCands = cms.bool(False)
)


process.hpsPFTauProducerSansRefs = cms.EDProducer("RecoBaseTauCleaner",
    cleaners = cms.VPSet(cms.PSet(
        name = cms.string('Charge'),
        nprongs = cms.vuint32(1, 3),
        passForCharge = cms.int32(1),
        plugin = cms.string('RecoBaseTauChargeCleanerPlugin'),
        selectionFailValue = cms.double(0)
    ), 
        cms.PSet(
            name = cms.string('HPS_Select'),
            plugin = cms.string('RecoBaseTauDiscriminantCleanerPlugin'),
            src = cms.InputTag("hpsSelectionDiscriminator")
        ), 
        cms.PSet(
            minTrackPt = cms.double(5.0),
            name = cms.string('killSoftTwoProngTaus'),
            plugin = cms.string('RecoBaseTauSoftTwoProngTausCleanerPlugin')
        ), 
        cms.PSet(
            name = cms.string('ChargedHadronMultiplicity'),
            plugin = cms.string('RecoBaseTauChargedHadronMultiplicityCleanerPlugin')
        ), 
        cms.PSet(
            name = cms.string('Pt'),
            plugin = cms.string('RecoBaseTauStringCleanerPlugin'),
            selection = cms.string('leadPFCand().isNonnull()'),
            selectionFailValue = cms.double(1000.0),
            selectionPassFunction = cms.string('-pt()'),
            tolerance = cms.double(0.01)
        ), 
        cms.PSet(
            name = cms.string('StripMultiplicity'),
            plugin = cms.string('RecoBaseTauStringCleanerPlugin'),
            selection = cms.string('leadPFCand().isNonnull()'),
            selectionFailValue = cms.double(1000.0),
            selectionPassFunction = cms.string('-signalPiZeroCandidates().size()')
        ), 
        cms.PSet(
            name = cms.string('CombinedIsolation'),
            plugin = cms.string('RecoBaseTauStringCleanerPlugin'),
            selection = cms.string('leadPFCand().isNonnull()'),
            selectionFailValue = cms.double(1000.0),
            selectionPassFunction = cms.string('isolationPFChargedHadrCandsPtSum() + isolationPFGammaCandsEtSum()')
        )),
    src = cms.InputTag("combinatoricRecoTaus")
)

process.hpsPFTauProducer = cms.EDProducer("RecoBaseTauPiZeroUnembedder",
    src = cms.InputTag("hpsPFTauProducerSansRefs")
)


def convertModuleToBaseTau(process, name):
    module = getattr(process, name)
    module.__dict__['_TypedParameterizable__type'] = module.type_().replace('RecoTau', 'RecoBaseTau')
    if hasattr(module, 'PFTauProducer'):
        module.PFBaseTauProducer = module.PFTauProducer
        # del module.PFTauProducer
    if hasattr(module, 'particleFlowSrc'):
        module.particleFlowSrc = cms.InputTag("packedPFCandidates", "", "")
    if hasattr(module, 'vertexSrc'):
        module.vertexSrc = cms.InputTag('offlineSlimmedPrimaryVertices')
    if hasattr(module, 'qualityCuts') and hasattr(module.qualityCuts, 'primaryVertexSrc'):
        module.qualityCuts.primaryVertexSrc = cms.InputTag('offlineSlimmedPrimaryVertices')

convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByDecayModeFindingNewDMs')
convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByDecayModeFindingOldDMs')
convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByDecayModeFinding')
convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByLooseChargedIsolation')
convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByLooseIsolation')
convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits')
convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits')
convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits')
convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits')

convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03')
convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03')
convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03')

convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByLoosePileupWeightedIsolation3Hits')
convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByMediumPileupWeightedIsolation3Hits')
convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByTightPileupWeightedIsolation3Hits')
convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByRawPileupWeightedIsolation3Hits')

convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone')

convertModuleToBaseTau(process, 'hpsPFTauChargedIsoPtSum')
convertModuleToBaseTau(process, 'hpsPFTauNeutralIsoPtSum')
convertModuleToBaseTau(process, 'hpsPFTauPUcorrPtSum')
convertModuleToBaseTau(process, 'hpsPFTauNeutralIsoPtSumWeight')
convertModuleToBaseTau(process, 'hpsPFTauFootprintCorrection')
convertModuleToBaseTau(process, 'hpsPFTauPhotonPtSumOutsideSignalCone')

process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByTightElectronRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMediumElectronRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByLooseElectronRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVA6rawElectronRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVA6VLooseElectronRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVA6LooseElectronRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVA6MediumElectronRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVA6TightElectronRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByMVA6VTightElectronRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByDeadECALElectronRejection)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByLooseMuonRejection3)
process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauDiscriminationByTightMuonRejection3)

# process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauVertexAndImpactParametersSeq)
# process.produceAndDiscriminateHPSPFTaus.remove(process.hpsPFTauMVAIsolation2Seq)
process.hpsPFTauVertexAndImpactParametersSeq = cms.Sequence(process.hpsPFTauTransverseImpactParameters)

process.hpsPFTauTransverseImpactParameters = cms.EDProducer("PFBaseTauTransverseImpactParameters",
    PFTauTag = cms.InputTag("hpsPFTauProducer"),
    useFullCalculation = cms.bool(True),
    leadingTrkOrPFCandOption = cms.string('leadPFCand'),
    primaryVertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    pvFindingAlgo = cms.string('closestInDeltaZ'),
    recoverLeadingTrk = cms.bool(False),
    vxAssocQualityCuts = cms.PSet(
        maxTrackChi2 = cms.double(100.0),
        maxTransverseImpactParameter = cms.double(0.1),
        minGammaEt = cms.double(0.5),
        minTrackHits = cms.uint32(3),
        minTrackPixelHits = cms.uint32(0),
        minTrackPt = cms.double(0.5),
        minTrackVertexWeight = cms.double(-1.0)
    ),
    vertexTrackFiltering = cms.bool(False),

)

convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw')
for wp in ['VLoose', 'Loose', 'Medium', 'Tight', 'VTight', 'VVTight']:
    convertModuleToBaseTau(process, 'hpsPFTauDiscriminationBy{}IsolationMVArun2v1DBoldDMwLT'.format(wp))

convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw')
for wp in ['VLoose', 'Loose', 'Medium', 'Tight', 'VTight', 'VVTight']:
    convertModuleToBaseTau(process, 'hpsPFTauDiscriminationBy{}IsolationMVArun2v1DBnewDMwLT'.format(wp))

convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw')
for wp in ['VLoose', 'Loose', 'Medium', 'Tight', 'VTight', 'VVTight']:
    convertModuleToBaseTau(process, 'hpsPFTauDiscriminationBy{}IsolationMVArun2v1PWoldDMwLT'.format(wp))

convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw')
for wp in ['VLoose', 'Loose', 'Medium', 'Tight', 'VTight', 'VVTight']:
    convertModuleToBaseTau(process, 'hpsPFTauDiscriminationBy{}IsolationMVArun2v1PWnewDMwLT'.format(wp))

convertModuleToBaseTau(process, 'hpsPFTauChargedIsoPtSumdR03')
convertModuleToBaseTau(process, 'hpsPFTauNeutralIsoPtSumdR03')
convertModuleToBaseTau(process, 'hpsPFTauPUcorrPtSumdR03')
convertModuleToBaseTau(process, 'hpsPFTauNeutralIsoPtSumWeightdR03')
convertModuleToBaseTau(process, 'hpsPFTauFootprintCorrectiondR03')
convertModuleToBaseTau(process, 'hpsPFTauPhotonPtSumOutsideSignalConedR03')


convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw')
for wp in ['VLoose', 'Loose', 'Medium', 'Tight', 'VTight', 'VVTight']:
    convertModuleToBaseTau(process, 'hpsPFTauDiscriminationBy{}IsolationMVArun2v1DBdR03oldDMwLT'.format(wp))

convertModuleToBaseTau(process, 'hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw')
for wp in ['VLoose', 'Loose', 'Medium', 'Tight', 'VTight', 'VVTight']:
    convertModuleToBaseTau(process, 'hpsPFTauDiscriminationBy{}IsolationMVArun2v1PWdR03oldDMwLT'.format(wp))

# OK NOW COMES PATTY PAT

# FIXME - check both if this is the OK collection...
process.tauGenJets.GenParticles = cms.InputTag("prunedGenParticles")
process.tauMatch.genParticles = cms.InputTag("prunedGenParticles")
process.tauMatch.matched = cms.InputTag("prunedGenParticles")


process.patTaus.__dict__['_TypedParameterizable__type'] = 'PATTauBaseProducer'
convertModuleToBaseTau(process, 'patTaus')

process.patTaus.tauIDSources = tauIDSources = cms.PSet(
        # againstElectronLooseMVA6 = cms.InputTag("hpsPFTauDiscriminationByMVA6LooseElectronRejection"),
        # againstElectronMVA6Raw = cms.InputTag("hpsPFTauDiscriminationByMVA6rawElectronRejection"),
        # againstElectronMVA6category = cms.InputTag("hpsPFTauDiscriminationByMVA6rawElectronRejection","category"),
        # againstElectronMediumMVA6 = cms.InputTag("hpsPFTauDiscriminationByMVA6MediumElectronRejection"),
        # againstElectronTightMVA6 = cms.InputTag("hpsPFTauDiscriminationByMVA6TightElectronRejection"),
        # againstElectronVLooseMVA6 = cms.InputTag("hpsPFTauDiscriminationByMVA6VLooseElectronRejection"),
        # againstElectronVTightMVA6 = cms.InputTag("hpsPFTauDiscriminationByMVA6VTightElectronRejection"),
        # againstMuonLoose3 = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection3"),
        # againstMuonTight3 = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection3"),
        byCombinedIsolationDeltaBetaCorrRaw3Hits = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits"),
        byIsolationMVArun2v1DBdR03oldDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw"),
        byIsolationMVArun2v1DBnewDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw"),
        byIsolationMVArun2v1DBoldDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw"),
        byIsolationMVArun2v1PWdR03oldDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw"),
        byIsolationMVArun2v1PWnewDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw"),
        byIsolationMVArun2v1PWoldDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw"),
        byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits"),
        byLooseIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT"),
        byLooseIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT"),
        byLooseIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT"),
        byLooseIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT"),
        byLooseIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT"),
        byLooseIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT"),
        byMediumCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits"),
        byMediumIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT"),
        byMediumIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT"),
        byMediumIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT"),
        byMediumIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT"),
        byMediumIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT"),
        byMediumIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT"),
        byPhotonPtSumOutsideSignalCone = cms.InputTag("hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone"),
        byTightCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits"),
        byTightIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT"),
        byTightIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT"),
        byTightIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT"),
        byTightIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT"),
        byTightIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT"),
        byTightIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT"),
        byVLooseIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT"),
        byVLooseIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT"),
        byVLooseIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT"),
        byVLooseIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT"),
        byVLooseIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT"),
        byVLooseIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT"),
        byVTightIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT"),
        byVTightIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT"),
        byVTightIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT"),
        byVTightIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT"),
        byVTightIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT"),
        byVTightIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT"),
        byVVTightIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT"),
        byVVTightIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT"),
        byVVTightIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT"),
        byVVTightIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT"),
        byVVTightIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT"),
        byVVTightIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT"),
        chargedIsoPtSum = cms.InputTag("hpsPFTauChargedIsoPtSum"),
        chargedIsoPtSumdR03 = cms.InputTag("hpsPFTauChargedIsoPtSumdR03"),
        decayModeFinding = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
        decayModeFindingNewDMs = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
        footprintCorrection = cms.InputTag("hpsPFTauFootprintCorrection"),
        footprintCorrectiondR03 = cms.InputTag("hpsPFTauFootprintCorrectiondR03"),
        neutralIsoPtSum = cms.InputTag("hpsPFTauNeutralIsoPtSum"),
        neutralIsoPtSumWeight = cms.InputTag("hpsPFTauNeutralIsoPtSumWeight"),
        neutralIsoPtSumWeightdR03 = cms.InputTag("hpsPFTauNeutralIsoPtSumWeightdR03"),
        neutralIsoPtSumdR03 = cms.InputTag("hpsPFTauNeutralIsoPtSumdR03"),
        photonPtSumOutsideSignalCone = cms.InputTag("hpsPFTauPhotonPtSumOutsideSignalCone"),
        photonPtSumOutsideSignalConedR03 = cms.InputTag("hpsPFTauPhotonPtSumOutsideSignalConedR03"),
        puCorrPtSum = cms.InputTag("hpsPFTauPUcorrPtSum")
    )

process.makePatTausTask = cms.Task(
    # reco pre-production
    process.patHPSPFTauDiscriminationTask,
    # patPFCandidateIsoDepositSelectionTask,
    # process.patPFTauIsolationTask,
    #patTauJetCorrections *
    # pat specifics
    process.tauMatch,
    process.tauGenJets,
    process.tauGenJetsSelectorAllHadrons,
    process.tauGenJetMatch,
    # object production
    process.patTaus,    
)
if jetCollection == 'patAK4PFJets':
	process.PATTauSequence.replace(process.recoTauAK4PatJets08Region,
				       process.ak4PFJetsPAT+
				       process.patAK4PFJets+
				       process.recoTauAK4PatJets08Region)

process.load('FWCore.MessageService.MessageLogger_cfi')
if process.maxEvents.input.value()>10:
     process.MessageLogger.cerr.FwkReport.reportEvery = process.maxEvents.input.value()//10
if process.maxEvents.input.value()>10000 or process.maxEvents.input.value()<0:
     process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
)
process.options.numberOfThreads=cms.untracked.uint32(4)
#process.options.numberOfThreads=cms.untracked.uint32(1)
process.options.numberOfStreams=cms.untracked.uint32(0)
print '\t No. of threads:', process.options.numberOfThreads.value(),', no. of streams:',process.options.numberOfStreams.value()

process.options = cms.untracked.PSet(
     process.options,
     wantSummary = cms.untracked.bool(True)
)
