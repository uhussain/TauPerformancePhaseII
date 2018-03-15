// -*- C++ -*-
//
// Package:    RecoTauTag/phase2Taus
// Class:      phase2Taus
// 
/**\class phase2Taus phase2Taus.cc RecoTauTag/phase2Taus/plugins/phase2Taus.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Isabel Ojalvo
//         Created:  Tue, 15 Nov 2016 16:00:32 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "iostream"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class phase2Taus : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit phase2Taus(const edm::ParameterSet&);
      ~phase2Taus();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
  edm::EDGetTokenT<std::vector<reco::Vertex> > vtxToken_;
  edm::EDGetTokenT<pat::TauCollection> tauToken_; 
  //edm::EDGetTokenT<pat::TauCollection> tauOrgToken_; //Fix for puppi
  std::string tauID_;
  edm::EDGetTokenT<std::vector <reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<std::vector<pat::Jet> > jetsAK4Token_;
  edm::EDGetTokenT<std::vector<reco::GenJet> > genJetsToken_;

  TTree* tree;
  //TTree* OrgTaustree; //Fix
  TTree* jetTree;
  //TTree* OrgTausjetTree;//Fix
  int     run_;
  int  event_;
  int     lumis_;
  double tauPt_;
  double tauEta_;
  double tauMass_;
  double genTauPt_;
  double genTauEta_;

  //tau ID and isolation
  bool   taupfTausDiscriminationByDecayModeFinding_;
  bool   taupfTausDiscriminationByDecayModeFindingNewDMs_;
  bool   tauByLooseCombinedIsolationDeltaBetaCorr3Hits_;
  bool   tauByMediumCombinedIsolationDeltaBetaCorr3Hits_;
  bool   tauByTightCombinedIsolationDeltaBetaCorr3Hits_;
  float   tauCombinedIsolationDeltaBetaCorrRaw3Hits_;

  float        tauByIsolationMVArun2v1DBnewDMwLTraw_;
  float        tauByIsolationMVArun2v1DBoldDMwLTraw_;
  float        tauByIsolationMVArun2v1PWnewDMwLTraw_;
  float        tauByIsolationMVArun2v1PWoldDMwLTraw_;

  //Tau Ingredients
  float tauChargedIsoPtSum_;
  float tauNeutralIsoPtSum_;
  float tauPuCorrPtSum_;
  float taufootprintCorrection_;
  float tauphotonPtSumOutsideSignalCone_;


  double jetPt_;
  double jetEta_;
  double vtxX_, vtxY_, vtxZ_;
  int nvtx_;
  int vtxIndex_;
  double dz_tt_;
  int dm_;
  int goodReco_;
  int genTauMatch_;
  int genJetMatch_;
  int jetTauMatch_;
  
  reco::Candidate::LorentzVector GetVisibleP4(std::vector<const reco::GenParticle*>& daughters);
  void findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters);
  bool isNeutrino(const reco::Candidate* daughter);
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
phase2Taus::phase2Taus(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertices"))),
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
//  tauOrgToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tausOrg"))),//Fix
  prunedGenToken_(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  jetsAK4Token_(consumes<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
  genJetsToken_(consumes<std::vector<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("genJets")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");

   tauID_    = iConfig.getParameter<std::string>("tauID");
   edm::Service<TFileService> fs;

   tree = fs->make<TTree>("ModFixedStripTaus", "ModFixedStripTaus");
   tree->Branch("run",     &run_);
   tree->Branch("event",   &event_);
   tree->Branch("lumis",   &lumis_);
   tree->Branch("tauPt", &tauPt_,"tauPt_/D");
   tree->Branch("tauEta", &tauEta_,"tauEta_/D");
   tree->Branch("genTauPt", &genTauPt_,"genTauPt_/D");
   tree->Branch("genTauEta", &genTauEta_,"genTauEta_/D");
   tree->Branch("genTauMatch", &genTauMatch_,"genTauMatch_/I");
   tree->Branch("nvtx",&nvtx_,"nvtx_/I");
   tree->Branch("vtxX",         &vtxX_,        "vtxX_/D"         );
   tree->Branch("vtxY",         &vtxY_,        "vtxY_/D"         );
   tree->Branch("vtxZ",         &vtxZ_,        "vtxZ_/D"         );
   tree->Branch("vtxIndex",     &vtxIndex_,    "vtxIndex_/I"     );
   tree->Branch("dz_tt", &dz_tt_, "dz_tt_/D"); //////Dem
   tree->Branch("dm",&dm_,"dm_/I");
   tree->Branch("goodReco",&goodReco_,"goodReco_/I");
   tree->Branch("tauMass",&tauMass_,"tauMass_/D");
   tree->Branch("taupfTausDiscriminationByDecayModeFinding", &taupfTausDiscriminationByDecayModeFinding_);
   tree->Branch("taupfTausDiscriminationByDecayModeFindingNewDMs", &taupfTausDiscriminationByDecayModeFindingNewDMs_);
   tree->Branch("tauByIsolationMVArun2v1DBnewDMwLTraw", &tauByIsolationMVArun2v1DBnewDMwLTraw_);
   tree->Branch("tauByIsolationMVArun2v1DBoldDMwLTraw", &tauByIsolationMVArun2v1DBoldDMwLTraw_);
   tree->Branch("tauByIsolationMVArun2v1PWnewDMwLTraw", &tauByIsolationMVArun2v1PWnewDMwLTraw_);
   tree->Branch("tauByIsolationMVArun2v1PWoldDMwLTraw", &tauByIsolationMVArun2v1PWoldDMwLTraw_);
   tree->Branch("tauByIsolationMVArun2v1DBnewDMwLTraw", &tauByIsolationMVArun2v1DBnewDMwLTraw_);
   tree->Branch("tauByIsolationMVArun2v1DBoldDMwLTraw", &tauByIsolationMVArun2v1DBoldDMwLTraw_);
   tree->Branch("tauByIsolationMVArun2v1PWnewDMwLTraw", &tauByIsolationMVArun2v1PWnewDMwLTraw_);
   tree->Branch("tauByIsolationMVArun2v1PWoldDMwLTraw", &tauByIsolationMVArun2v1PWoldDMwLTraw_);
   tree->Branch("tauChargedIsoPtSum"  ,&tauChargedIsoPtSum_);
   tree->Branch("tauNeutralIsoPtSum"  ,&tauNeutralIsoPtSum_);
   tree->Branch("tauByLooseCombinedIsolationDeltaBetaCorr3Hits", &tauByLooseCombinedIsolationDeltaBetaCorr3Hits_);
   tree->Branch("tauByMediumCombinedIsolationDeltaBetaCorr3Hits", &tauByMediumCombinedIsolationDeltaBetaCorr3Hits_);
   tree->Branch("tauByTightCombinedIsolationDeltaBetaCorr3Hits", &tauByTightCombinedIsolationDeltaBetaCorr3Hits_);
   tree->Branch("tauCombinedIsolationDeltaBetaCorrRaw3Hits", &tauCombinedIsolationDeltaBetaCorrRaw3Hits_);
   tree->Branch("tauPuCorrPtSum"  ,&tauPuCorrPtSum_);
   tree->Branch("taufootprintCorrection"  ,&taufootprintCorrection_);
   tree->Branch("tauphotonPtSumOutsideSignalCone"  ,&tauphotonPtSumOutsideSignalCone_);
   /*
   OrgTaustree = fs->make<TTree>("OrgPFTaus", "OrgPFTaus");
   OrgTaustree->Branch("run",     &run_);
   OrgTaustree->Branch("event",   &event_);
   OrgTaustree->Branch("lumis",   &lumis_);
   OrgTaustree->Branch("tauPt", &tauPt_,"tauPt_/D");
   OrgTaustree->Branch("tauEta", &tauEta_,"tauEta_/D");
   OrgTaustree->Branch("genTauPt", &genTauPt_,"genTauPt_/D");
   OrgTaustree->Branch("genTauEta", &genTauEta_,"genTauEta_/D");
   OrgTaustree->Branch("genTauMatch", &genTauMatch_,"genTauMatch_/I");
   OrgTaustree->Branch("nvtx",&nvtx_,"nvtx_/I");
   OrgTaustree->Branch("vtxX",         &vtxX_,        "vtxX_/D"         );
   OrgTaustree->Branch("vtxY",         &vtxY_,        "vtxY_/D"         );
   OrgTaustree->Branch("vtxZ",         &vtxZ_,        "vtxZ_/D"         );
   OrgTaustree->Branch("vtxIndex",     &vtxIndex_,    "vtxIndex_/I"     );
   OrgTaustree->Branch("dm",&dm_,"dm_/I");
   OrgTaustree->Branch("goodReco",&goodReco_,"goodReco_/I");
   OrgTaustree->Branch("tauMass",&tauMass_,"tauMass_/D");
   OrgTaustree->Branch("taupfTausDiscriminationByDecayModeFinding", &taupfTausDiscriminationByDecayModeFinding_);
   OrgTaustree->Branch("taupfTausDiscriminationByDecayModeFindingNewDMs", &taupfTausDiscriminationByDecayModeFindingNewDMs_);
   OrgTaustree->Branch("tauByIsolationMVArun2v1DBnewDMwLTraw", &tauByIsolationMVArun2v1DBnewDMwLTraw_);
   OrgTaustree->Branch("tauByIsolationMVArun2v1DBoldDMwLTraw", &tauByIsolationMVArun2v1DBoldDMwLTraw_);
   OrgTaustree->Branch("tauByIsolationMVArun2v1PWnewDMwLTraw", &tauByIsolationMVArun2v1PWnewDMwLTraw_);
   OrgTaustree->Branch("tauByIsolationMVArun2v1PWoldDMwLTraw", &tauByIsolationMVArun2v1PWoldDMwLTraw_);
   OrgTaustree->Branch("tauByIsolationMVArun2v1DBnewDMwLTraw", &tauByIsolationMVArun2v1DBnewDMwLTraw_);
   OrgTaustree->Branch("tauByIsolationMVArun2v1DBoldDMwLTraw", &tauByIsolationMVArun2v1DBoldDMwLTraw_);
   OrgTaustree->Branch("tauByIsolationMVArun2v1PWnewDMwLTraw", &tauByIsolationMVArun2v1PWnewDMwLTraw_);
   OrgTaustree->Branch("tauByIsolationMVArun2v1PWoldDMwLTraw", &tauByIsolationMVArun2v1PWoldDMwLTraw_);
   OrgTaustree->Branch("tauChargedIsoPtSum"  ,&tauChargedIsoPtSum_);
   OrgTaustree->Branch("tauNeutralIsoPtSum"  ,&tauNeutralIsoPtSum_);
   OrgTaustree->Branch("tauByLooseCombinedIsolationDeltaBetaCorr3Hits", &tauByLooseCombinedIsolationDeltaBetaCorr3Hits_);
   OrgTaustree->Branch("tauByMediumCombinedIsolationDeltaBetaCorr3Hits", &tauByMediumCombinedIsolationDeltaBetaCorr3Hits_);
   OrgTaustree->Branch("tauByTightCombinedIsolationDeltaBetaCorr3Hits", &tauByTightCombinedIsolationDeltaBetaCorr3Hits_);
   OrgTaustree->Branch("tauCombinedIsolationDeltaBetaCorrRaw3Hits", &tauCombinedIsolationDeltaBetaCorrRaw3Hits_);
   OrgTaustree->Branch("tauPuCorrPtSum"  ,&tauPuCorrPtSum_);
   OrgTaustree->Branch("taufootprintCorrection"  ,&taufootprintCorrection_);
   OrgTaustree->Branch("tauphotonPtSumOutsideSignalCone"  ,&tauphotonPtSumOutsideSignalCone_);
*/ //Fix for Puppi
   jetTree = fs->make<TTree>(      "jetModFixedStripTaus",   "jetModFixedStripTaus"       );
   jetTree->Branch("run",     &run_);
   jetTree->Branch("event",   &event_);
   jetTree->Branch("lumis",   &lumis_);  
   jetTree->Branch("tauPt",        &tauPt_,       "tauPt_/D"        );
   jetTree->Branch("tauEta",       &tauEta_,      "tauEta_/D"       );
   jetTree->Branch("jetPt",        &jetPt_,       "jetPt_/D"        );
   jetTree->Branch("jetEta",       &jetEta_,      "jetEta_/D"       );
   jetTree->Branch("genJetMatch",  &genJetMatch_, "genJetMatch_/I"  );
   jetTree->Branch("jetTauMatch",  &jetTauMatch_, "jetTauMatch_/I"  );
   jetTree->Branch("nvtx",         &nvtx_,        "nvtx_/I"         );
   jetTree->Branch("vtxX",         &vtxX_,        "vtxX_/D"         );
   jetTree->Branch("vtxY",         &vtxY_,        "vtxY_/D"         );
   jetTree->Branch("vtxZ",         &vtxZ_,        "vtxZ_/D"         ); 
   jetTree->Branch("dz_tt", &dz_tt_, "dz_tt_/D");
   //jetTree->Branch("vtxIndex",     &vtxIndex_,    "vtxIndex_/I"     );
   jetTree->Branch("dm",          &dm_,         "dm_/I"          );
   jetTree->Branch("tauMass",      &tauMass_,     "tauMass_/D"      );
   jetTree->Branch("taupfTausDiscriminationByDecayModeFinding", &taupfTausDiscriminationByDecayModeFinding_);
   jetTree->Branch("taupfTausDiscriminationByDecayModeFindingNewDMs", &taupfTausDiscriminationByDecayModeFindingNewDMs_);
   jetTree->Branch("tauByIsolationMVArun2v1DBnewDMwLTraw", &tauByIsolationMVArun2v1DBnewDMwLTraw_);
   jetTree->Branch("tauByIsolationMVArun2v1DBoldDMwLTraw", &tauByIsolationMVArun2v1DBoldDMwLTraw_);
   jetTree->Branch("tauByIsolationMVArun2v1PWnewDMwLTraw", &tauByIsolationMVArun2v1PWnewDMwLTraw_);
   jetTree->Branch("tauByIsolationMVArun2v1PWoldDMwLTraw", &tauByIsolationMVArun2v1PWoldDMwLTraw_);
   jetTree->Branch("tauByIsolationMVArun2v1DBnewDMwLTraw", &tauByIsolationMVArun2v1DBnewDMwLTraw_);
   jetTree->Branch("tauByIsolationMVArun2v1DBoldDMwLTraw", &tauByIsolationMVArun2v1DBoldDMwLTraw_);
   jetTree->Branch("tauByIsolationMVArun2v1PWnewDMwLTraw", &tauByIsolationMVArun2v1PWnewDMwLTraw_);
   jetTree->Branch("tauByIsolationMVArun2v1PWoldDMwLTraw", &tauByIsolationMVArun2v1PWoldDMwLTraw_);
   jetTree->Branch("tauChargedIsoPtSum"  ,&tauChargedIsoPtSum_);
   jetTree->Branch("tauNeutralIsoPtSum"  ,&tauNeutralIsoPtSum_);
   jetTree->Branch("tauByLooseCombinedIsolationDeltaBetaCorr3Hits", &tauByLooseCombinedIsolationDeltaBetaCorr3Hits_);
   jetTree->Branch("tauByMediumCombinedIsolationDeltaBetaCorr3Hits", &tauByMediumCombinedIsolationDeltaBetaCorr3Hits_);
   jetTree->Branch("tauByTightCombinedIsolationDeltaBetaCorr3Hits", &tauByTightCombinedIsolationDeltaBetaCorr3Hits_);
   jetTree->Branch("tauCombinedIsolationDeltaBetaCorrRaw3Hits", &tauCombinedIsolationDeltaBetaCorrRaw3Hits_);
   jetTree->Branch("tauPuCorrPtSum"  ,&tauPuCorrPtSum_);
   jetTree->Branch("taufootprintCorrection"  ,&taufootprintCorrection_);
   jetTree->Branch("tauphotonPtSumOutsideSignalCone"  ,&tauphotonPtSumOutsideSignalCone_);
   /*
   OrgTausjetTree = fs->make<TTree>(      "jetOrgPFTaus",   "jetOrgPFTaus"       );
   OrgTausjetTree->Branch("run",     &run_);
   OrgTausjetTree->Branch("event",   &event_);
   OrgTausjetTree->Branch("lumis",   &lumis_);  
   OrgTausjetTree->Branch("tauPt",        &tauPt_,       "tauPt_/D"        );
   OrgTausjetTree->Branch("tauEta",       &tauEta_,      "tauEta_/D"       );
   OrgTausjetTree->Branch("jetPt",        &jetPt_,       "jetPt_/D"        );
   OrgTausjetTree->Branch("jetEta",       &jetEta_,      "jetEta_/D"       );
   OrgTausjetTree->Branch("genJetMatch",  &genJetMatch_, "genJetMatch_/I"  );
   OrgTausjetTree->Branch("jetTauMatch",  &jetTauMatch_, "jetTauMatch_/I"  );
   OrgTausjetTree->Branch("nvtx",         &nvtx_,        "nvtx_/I"         );
   OrgTausjetTree->Branch("vtxX",         &vtxX_,        "vtxX_/D"         );
   OrgTausjetTree->Branch("vtxY",         &vtxY_,        "vtxY_/D"         );
   OrgTausjetTree->Branch("vtxZ",         &vtxZ_,        "vtxZ_/D"         ); 
   //OrgTausjetTree->Branch("vtxIndex",     &vtxIndex_,    "vtxIndex_/I"     );
   OrgTausjetTree->Branch("dm",          &dm_,         "dm_/I"          );
   OrgTausjetTree->Branch("tauMass",      &tauMass_,     "tauMass_/D"      );
   OrgTausjetTree->Branch("taupfTausDiscriminationByDecayModeFinding", &taupfTausDiscriminationByDecayModeFinding_);
   OrgTausjetTree->Branch("taupfTausDiscriminationByDecayModeFindingNewDMs", &taupfTausDiscriminationByDecayModeFindingNewDMs_);
   OrgTausjetTree->Branch("tauByIsolationMVArun2v1DBnewDMwLTraw", &tauByIsolationMVArun2v1DBnewDMwLTraw_);
   OrgTausjetTree->Branch("tauByIsolationMVArun2v1DBoldDMwLTraw", &tauByIsolationMVArun2v1DBoldDMwLTraw_);
   OrgTausjetTree->Branch("tauByIsolationMVArun2v1PWnewDMwLTraw", &tauByIsolationMVArun2v1PWnewDMwLTraw_);
   OrgTausjetTree->Branch("tauByIsolationMVArun2v1PWoldDMwLTraw", &tauByIsolationMVArun2v1PWoldDMwLTraw_);
   OrgTausjetTree->Branch("tauByIsolationMVArun2v1DBnewDMwLTraw", &tauByIsolationMVArun2v1DBnewDMwLTraw_);
   OrgTausjetTree->Branch("tauByIsolationMVArun2v1DBoldDMwLTraw", &tauByIsolationMVArun2v1DBoldDMwLTraw_);
   OrgTausjetTree->Branch("tauByIsolationMVArun2v1PWnewDMwLTraw", &tauByIsolationMVArun2v1PWnewDMwLTraw_);
   OrgTausjetTree->Branch("tauByIsolationMVArun2v1PWoldDMwLTraw", &tauByIsolationMVArun2v1PWoldDMwLTraw_);
   OrgTausjetTree->Branch("tauChargedIsoPtSum"  ,&tauChargedIsoPtSum_);
   OrgTausjetTree->Branch("tauNeutralIsoPtSum"  ,&tauNeutralIsoPtSum_);
   OrgTausjetTree->Branch("tauByLooseCombinedIsolationDeltaBetaCorr3Hits", &tauByLooseCombinedIsolationDeltaBetaCorr3Hits_);
   OrgTausjetTree->Branch("tauByMediumCombinedIsolationDeltaBetaCorr3Hits", &tauByMediumCombinedIsolationDeltaBetaCorr3Hits_);
   OrgTausjetTree->Branch("tauByTightCombinedIsolationDeltaBetaCorr3Hits", &tauByTightCombinedIsolationDeltaBetaCorr3Hits_);
   OrgTausjetTree->Branch("tauCombinedIsolationDeltaBetaCorrRaw3Hits", &tauCombinedIsolationDeltaBetaCorrRaw3Hits_);
   OrgTausjetTree->Branch("tauPuCorrPtSum"  ,&tauPuCorrPtSum_);
   OrgTausjetTree->Branch("taufootprintCorrection"  ,&taufootprintCorrection_);
   OrgTausjetTree->Branch("tauphotonPtSumOutsideSignalCone"  ,&tauphotonPtSumOutsideSignalCone_);
*/ //Fix for Puppi
}


phase2Taus::~phase2Taus()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
phase2Taus::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   Handle<std::vector<reco::Vertex> > vertices;
   iEvent.getByToken(vtxToken_, vertices);
   nvtx_=vertices->size();
   
   Handle<pat::TauCollection> taus;
   iEvent.getByToken(tauToken_, taus);
   
//   Handle<pat::TauCollection> tausOrg; //Fix for Puppi
//   iEvent.getByToken(tauOrgToken_, tausOrg); //Fix for Puppi

   Handle<reco::GenJetCollection> genJets;
   iEvent.getByToken(genJetsToken_,genJets);
   
   Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByToken(prunedGenToken_, genParticles);
   
   Handle<std::vector<pat::Jet> > jetHandle;
   iEvent.getByToken(jetsAK4Token_, jetHandle);

   if (!jetHandle.isValid()) {
        edm::LogWarning("ggNtuplizer") << "no pat::Jets (AK4) in event";
        return;
    }
    
   run_    = iEvent.id().run();
   event_  = iEvent.id().event();
   lumis_  = iEvent.luminosityBlock();
   
   std::vector<const reco::GenParticle*> GenTaus;
   std::vector<const reco::GenParticle*> GenEles;
   std::vector<const reco::GenParticle*> GenMus;
   for(std::vector<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++){
      if(TMath::Abs(genParticle->pdgId()) == 15 && genParticle->isLastCopy() && genParticle->statusFlags().fromHardProcess()) GenTaus.push_back(&(*genParticle)); 
      if(TMath::Abs(genParticle->pdgId()) == 11 && genParticle->isLastCopy() && genParticle->statusFlags().fromHardProcess()) GenEles.push_back(&(*genParticle));
      if(TMath::Abs(genParticle->pdgId()) == 13 && genParticle->isLastCopy() && genParticle->statusFlags().fromHardProcess()) GenMus.push_back(&(*genParticle));
     //if(TMath::Abs(genParticle->pdgId()) == 11) GenEles.push_back(&(*genParticle));
     //if(TMath::Abs(genParticle->pdgId()) == 13) GenMus.push_back(&(*genParticle));
   }

   std::vector<pat::Jet> Jets;
   for (unsigned int iJet = 0; iJet < jetHandle->size() ; ++iJet){
     pat::JetRef jetCand(jetHandle, iJet);
     if(jetCand->pt() < 18 )continue;
     bool isATau=false;
     for(auto genTau : GenTaus){
       std::vector<const reco::GenParticle*> genTauDaughters;
       findDaughters(genTau, genTauDaughters);
       reco::Candidate::LorentzVector genTauVis = GetVisibleP4(genTauDaughters);
       genTauPt_  = (float) genTauVis.pt();
       genTauEta_ = (float) genTauVis.eta();
       if (reco::deltaR(jetCand->eta(),jetCand->phi(),genTauVis.eta(),genTauVis.phi()) < 0.5)
	        isATau=true;
       
     }
     bool isAEle=false;
     for(auto genEle : GenEles){
       if (reco::deltaR(jetCand->eta(),jetCand->phi(),genEle->eta(),genEle->phi()) < 0.5)
	        isAEle=true;
     }
     bool isAMu=false;
     for(auto genMu : GenMus){
       if (reco::deltaR(jetCand->eta(),jetCand->phi(),genMu->eta(),genMu->phi()) < 0.5)
	        isAMu=true;
     }
     if(!isATau && !isAEle && !isAMu)
       Jets.push_back(*jetCand);
   }   
   std::cout<<run_<<":"<<event_<<":"<<lumis_<<std::endl;
   std::cout<<"RecoTausSize: " <<taus->size()<<std::endl;
   std::cout<<"GenTausSize: "<<GenTaus.size()<<std::endl;
   for(auto genTau : GenTaus){
     tauPt_=-10;
     tauEta_=-10;
     tauMass_=-10;
     
     taupfTausDiscriminationByDecayModeFinding_=0;
     taupfTausDiscriminationByDecayModeFindingNewDMs_=0; 
     tauByLooseCombinedIsolationDeltaBetaCorr3Hits_=0;
     tauByMediumCombinedIsolationDeltaBetaCorr3Hits_=0;
     tauByTightCombinedIsolationDeltaBetaCorr3Hits_=0;
       
     tauCombinedIsolationDeltaBetaCorrRaw3Hits_=-10;  
     tauByIsolationMVArun2v1DBnewDMwLTraw_=-10;
     tauByIsolationMVArun2v1DBoldDMwLTraw_=-10;
     tauByIsolationMVArun2v1PWnewDMwLTraw_=-10;
     tauByIsolationMVArun2v1PWoldDMwLTraw_=-10;
     
     tauChargedIsoPtSum_=-10;
     tauNeutralIsoPtSum_=-10;
     tauPuCorrPtSum_=-10;

     taufootprintCorrection_=-10;
     tauphotonPtSumOutsideSignalCone_=-10;

     genTauPt_=-10;
     genTauEta_=-10;
     
     nvtx_=-10;
     dm_=-10;
     dz_tt_=-10;
     goodReco_=0;
     genTauMatch_=0;

     std::vector<const reco::GenParticle*> genTauDaughters;
     findDaughters(genTau, genTauDaughters);
     reco::Candidate::LorentzVector genTauVis = GetVisibleP4(genTauDaughters);
     genTauPt_  = (float) genTauVis.pt();
     genTauEta_ = (float) genTauVis.eta();

     //std::cout<<" pt "<<genTauPt_<<" eta "<<genTauEta_<<std::endl;
     for(const pat::Tau &tau : *taus){
       for(auto cand : tau.isolationChargedHadrCands()){  ///dem
         if (abs(cand->charge())<0)continue;   ///dem >or < Fix
         if (cand->pt()<=0.5 or cand->dxy(vertices[tau_vertex_idxpf].position())>=0.1)continue; /////FIX ME dem
         if (cand->hasTrackDetails()){//dem
            tt = cand->pseudoTrack();//dem
            if (tt.normalizedChi2()>=100. or cand.numberOfHits()<3)continue;//dem
	 }//dem
         dz_tt =cand.dz(vertices[tau_vertex_idxpf].position()); //dem
       } //dem
       if (reco::deltaR(tau.eta(),tau.phi(),genTauVis.eta(),genTauVis.phi()) < 0.5 && tau.tauID("decayModeFinding")>-1 && tau.tauID("chargedIsoPtSum")>-1){ 
        //std::cout<<"RecoTauMatched"<<std::endl;
        genTauMatch_ = 1;
	      tauPt_  =tau.pt();
	      tauEta_ =tau.eta();
        tauMass_=tau.mass();
	      dm_ =tau.decayMode();


        taupfTausDiscriminationByDecayModeFinding_=tau.tauID("decayModeFinding");
        taupfTausDiscriminationByDecayModeFindingNewDMs_=tau.tauID("decayModeFindingNewDMs");   
        tauByLooseCombinedIsolationDeltaBetaCorr3Hits_=tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
        tauByMediumCombinedIsolationDeltaBetaCorr3Hits_=tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
        tauByTightCombinedIsolationDeltaBetaCorr3Hits_=tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
        tauCombinedIsolationDeltaBetaCorrRaw3Hits_=tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        tauByIsolationMVArun2v1DBnewDMwLTraw_=tau.tauID("byIsolationMVArun2v1DBnewDMwLTraw");
        tauByIsolationMVArun2v1DBoldDMwLTraw_=tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
        tauByIsolationMVArun2v1PWnewDMwLTraw_=tau.tauID("byIsolationMVArun2v1PWnewDMwLTraw");
        tauByIsolationMVArun2v1PWoldDMwLTraw_=tau.tauID("byIsolationMVArun2v1PWoldDMwLTraw");

        tauChargedIsoPtSum_=tau.tauID("chargedIsoPtSum");
        tauNeutralIsoPtSum_=tau.tauID("neutralIsoPtSum");
        tauPuCorrPtSum_=tau.tauID("puCorrPtSum");
        taufootprintCorrection_=tau.tauID("footprintCorrection");
        tauphotonPtSumOutsideSignalCone_=tau.tauID("photonPtSumOutsideSignalCone");

        goodReco_ = (bool) tau.tauID(tauID_) >0.5;

        //get the matched vertex
        //int vtx_index = -1;
        //float max_weight = 0.f;
        //for( unsigned i = 0; i < vertices->size(); ++i ) {
        //  const auto& vtx = (*vertices)[i];     
        //  //std::cout<<"fine here"<<std::endl;
        //  const float weight = vtx.trackWeight(tau.leadChargedHadrCand()->trackRef());// check me -> Get track ref for charged hadron candidate 
        // // std::cout<<"fails here"<<std::endl;
        //  if( weight > max_weight ) {
        //    max_weight = weight;
        //    vtx_index = i;
        //   }
        //      }

       // if (vertices->size()>0) {
       //       pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
       //       float max_weight = 0.f;
       //       for( unsigned i = 0; i < vertices->size(); ++i ) {
       //         const auto& vtx = (*vertices)[i];     
       //         //std::cout<<"fine here"<<std::endl;
       //         reco::TrackBaseRef const& trk = dynamic_cast<reco::TrackBaseRef const&>(packedLeadTauCand->bestTrack());
       //         //const auto& trk = *(packedLeadTauCand->bestTrack());
       //         const float weight = vtx.trackWeight(trk);// check me -> Get track ref for charged hadron candidate 
       //        // std::cout<<"fails here"<<std::endl;
       //         if( weight > max_weight ) {
       //           max_weight = weight;
       //           vtx_index = i;
       //          }
       //             }
        /*
        pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
        const auto& TauVtx = packedLeadTauCand->vertexRef();


        for( unsigned i = 0; i < vertices->size(); ++i ) {
          const auto& vtx = (*vertices)[i]; 
          if(vtx.x()==TauVtx->x() && vtx.y()==TauVtx->y() && vtx.z()==TauVtx->z()){
            vtx_index=i;
          }
        }
        //const auto& TauVtx= tau.leadChargedHadrCand()->vertexRef();
        std::cout<<"vtx_index: "<<vtx_index<<std::endl;
	      //now do vtx variable filling
	      vtxIndex_ = vtx_index;
	      const reco::Vertex& vtx = (vtx_index == -1 ? (*vertices)[0] : (*vertices)[vtx_index]);
	      vtxX_ = vtx.x();
	      vtxY_ = vtx.y();
	      vtxZ_ = vtx.z();
	      break;
*/ //Fix for Puppi
            }
          }

    tree->Fill(); 
   }
/*
   //duplicating the code above for OrigTaus (built with Orig PFTaus from AOD)
   for(auto genTau : GenTaus){
     tauPt_=-10;
     tauEta_=-10;
     tauMass_=-10;
     
     taupfTausDiscriminationByDecayModeFinding_=0;
     taupfTausDiscriminationByDecayModeFindingNewDMs_=0; 
     tauByLooseCombinedIsolationDeltaBetaCorr3Hits_=0;
     tauByMediumCombinedIsolationDeltaBetaCorr3Hits_=0;
     tauByTightCombinedIsolationDeltaBetaCorr3Hits_=0;
       
     tauCombinedIsolationDeltaBetaCorrRaw3Hits_=-10;  
     tauByIsolationMVArun2v1DBnewDMwLTraw_=-10;
     tauByIsolationMVArun2v1DBoldDMwLTraw_=-10;
     tauByIsolationMVArun2v1PWnewDMwLTraw_=-10;
     tauByIsolationMVArun2v1PWoldDMwLTraw_=-10;
     
     tauChargedIsoPtSum_=-10;
     tauNeutralIsoPtSum_=-10;
     tauPuCorrPtSum_=-10;

     taufootprintCorrection_=-10;
     tauphotonPtSumOutsideSignalCone_=-10;
     genTauPt_=-10;
     genTauEta_=-10;
     
     nvtx_=-10;
     dm_=-10;
     goodReco_=0;
     genTauMatch_=0;

     std::vector<const reco::GenParticle*> genTauDaughters;
     findDaughters(genTau, genTauDaughters);
     reco::Candidate::LorentzVector genTauVis = GetVisibleP4(genTauDaughters);
     genTauPt_  = (float) genTauVis.pt();
     genTauEta_ = (float) genTauVis.eta();

     //std::cout<<" pt "<<genTauPt_<<" eta "<<genTauEta_<<std::endl;
     for(const pat::Tau &tau : *tausOrg){
       if (reco::deltaR(tau.eta(),tau.phi(),genTauVis.eta(),genTauVis.phi()) < 0.5 && tau.tauID("decayModeFinding")>-1 && tau.tauID("chargedIsoPtSum")>-1){ 
        std::cout<<"RecoTauMatched"<<std::endl;
        genTauMatch_ = 1;
	      tauPt_  =tau.pt();
	      tauEta_ =tau.eta();
        tauMass_=tau.mass();
	      dm_ =tau.decayMode();

        taupfTausDiscriminationByDecayModeFinding_=tau.tauID("decayModeFinding");
        taupfTausDiscriminationByDecayModeFindingNewDMs_=tau.tauID("decayModeFindingNewDMs");   
        tauByLooseCombinedIsolationDeltaBetaCorr3Hits_=tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
        tauByMediumCombinedIsolationDeltaBetaCorr3Hits_=tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
        tauByTightCombinedIsolationDeltaBetaCorr3Hits_=tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
        tauCombinedIsolationDeltaBetaCorrRaw3Hits_=tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        tauByIsolationMVArun2v1DBnewDMwLTraw_=tau.tauID("byIsolationMVArun2v1DBnewDMwLTraw");
        tauByIsolationMVArun2v1DBoldDMwLTraw_=tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
        tauByIsolationMVArun2v1PWnewDMwLTraw_=tau.tauID("byIsolationMVArun2v1PWnewDMwLTraw");
        tauByIsolationMVArun2v1PWoldDMwLTraw_=tau.tauID("byIsolationMVArun2v1PWoldDMwLTraw");

        tauChargedIsoPtSum_=tau.tauID("chargedIsoPtSum");
        tauNeutralIsoPtSum_=tau.tauID("neutralIsoPtSum");
        tauPuCorrPtSum_=tau.tauID("puCorrPtSum");
        taufootprintCorrection_=tau.tauID("footprintCorrection");
        tauphotonPtSumOutsideSignalCone_=tau.tauID("photonPtSumOutsideSignalCone");
	      goodReco_ = (bool) tau.tauID(tauID_) >0.5;

        //get the matched vertex
        int vtx_index = -1;
        //float max_weight = 0.f;
        //for( unsigned i = 0; i < vertices->size(); ++i ) {
        //  const auto& vtx = (*vertices)[i];     
        //  //std::cout<<"fine here"<<std::endl;
        //  const float weight = vtx.trackWeight(tau.leadChargedHadrCand()->trackRef());// check me -> Get track ref for charged hadron candidate 
        // // std::cout<<"fails here"<<std::endl;
        //  if( weight > max_weight ) {
        //    max_weight = weight;
        //    vtx_index = i;
        //   }
        //      }

       // if (vertices->size()>0) {
       //       pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
       //       float max_weight = 0.f;
       //       for( unsigned i = 0; i < vertices->size(); ++i ) {
       //         const auto& vtx = (*vertices)[i];     
       //         //std::cout<<"fine here"<<std::endl;
       //         reco::TrackBaseRef const& trk = dynamic_cast<reco::TrackBaseRef const&>(packedLeadTauCand->bestTrack());
       //         //const auto& trk = *(packedLeadTauCand->bestTrack());
       //         const float weight = vtx.trackWeight(trk);// check me -> Get track ref for charged hadron candidate 
       //        // std::cout<<"fails here"<<std::endl;
       //         if( weight > max_weight ) {
       //           max_weight = weight;
       //           vtx_index = i;
       //          }
       //             }
        
        pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
        const auto& TauVtx = packedLeadTauCand->vertexRef();

        for( unsigned i = 0; i < vertices->size(); ++i ) {
          const auto& vtx = (*vertices)[i]; 
          if(vtx.x()==TauVtx->x() && vtx.y()==TauVtx->y() && vtx.z()==TauVtx->z()){
            vtx_index=i;
          }
        }
        //const auto& TauVtx=tau.leadChargedHadrCand()->vertexRef();
        std::cout<<"vtx_index: "<<vtx_index<<std::endl;
	      //now do vtx variable filling
	      vtxIndex_ = vtx_index;
	      const reco::Vertex& vtx = (vtx_index == -1 ? (*vertices)[0] : (*vertices)[vtx_index]);
	      vtxX_ = vtx.x();
	      vtxY_ = vtx.y();
	      vtxZ_ = vtx.z();
	      break;
            }
          }
    OrgTaustree->Fill(); 
   }
   //jetTree for ModFixedStripTaus
   for(auto jet : Jets){
     tauPt_=-10;
     tauEta_=-10;
     tauMass_=-10;

     taupfTausDiscriminationByDecayModeFinding_=0;
     taupfTausDiscriminationByDecayModeFindingNewDMs_=0; 
     tauByLooseCombinedIsolationDeltaBetaCorr3Hits_=0;
     tauByMediumCombinedIsolationDeltaBetaCorr3Hits_=0;
     tauByTightCombinedIsolationDeltaBetaCorr3Hits_=0;
       
     tauCombinedIsolationDeltaBetaCorrRaw3Hits_=-10;  
     tauByIsolationMVArun2v1DBnewDMwLTraw_=-10;
     tauByIsolationMVArun2v1DBoldDMwLTraw_=-10;
     tauByIsolationMVArun2v1PWnewDMwLTraw_=-10;
     tauByIsolationMVArun2v1PWoldDMwLTraw_=-10;
     
     tauChargedIsoPtSum_=-10;
     tauNeutralIsoPtSum_=-10;
     tauPuCorrPtSum_=-10;

     taufootprintCorrection_=-10;
     tauphotonPtSumOutsideSignalCone_=-10;
     jetPt_=jet.pt();
     jetEta_=jet.eta();
     
     //nvtx_=-10;
     dm_=-10;
     jetTauMatch_=0;
     genJetMatch_ = 0;

     //for (unsigned int iGenJet = 0; iGenJet < genJets->size() ; ++iGenJet){
       //reco::GenJetRef genJet(genJets, iGenJet);
      // if (reco::deltaR(genJet->eta(),genJet->phi(),jet.eta(),jet.phi()) < 0.4)
        // genJetMatch_ = 1;
    //}

     for(const pat::Tau &tau : *taus){
       if (reco::deltaR(tau.eta(),tau.phi(),jet.eta(),jet.phi()) < 0.3 && tau.tauID("decayModeFinding")>-1 && tau.tauID("chargedIsoPtSum")>-1){
	      jetTauMatch_ = 1;
	      tauPt_  =tau.pt();
	      tauEta_ =tau.eta();
        tauMass_=tau.mass();
	      dm_ =tau.decayMode();

        taupfTausDiscriminationByDecayModeFinding_=tau.tauID("decayModeFinding");
        taupfTausDiscriminationByDecayModeFindingNewDMs_=tau.tauID("decayModeFindingNewDMs");   
        tauByLooseCombinedIsolationDeltaBetaCorr3Hits_=tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
        tauByMediumCombinedIsolationDeltaBetaCorr3Hits_=tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
        tauByTightCombinedIsolationDeltaBetaCorr3Hits_=tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
        tauCombinedIsolationDeltaBetaCorrRaw3Hits_=tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        tauByIsolationMVArun2v1DBnewDMwLTraw_=tau.tauID("byIsolationMVArun2v1DBnewDMwLTraw");
        tauByIsolationMVArun2v1DBoldDMwLTraw_=tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
        tauByIsolationMVArun2v1PWnewDMwLTraw_=tau.tauID("byIsolationMVArun2v1PWnewDMwLTraw");
        tauByIsolationMVArun2v1PWoldDMwLTraw_=tau.tauID("byIsolationMVArun2v1PWoldDMwLTraw");

        tauChargedIsoPtSum_=tau.tauID("chargedIsoPtSum");
        tauNeutralIsoPtSum_=tau.tauID("neutralIsoPtSum");
        tauPuCorrPtSum_=tau.tauID("puCorrPtSum");
        taufootprintCorrection_=tau.tauID("footprintCorrection");
        tauphotonPtSumOutsideSignalCone_=tau.tauID("photonPtSumOutsideSignalCone");
        //get the matched vertex
        int vtx_index = -1;
        //float max_weight = 0.f;
        //for( unsigned i = 0; i < vertices->size(); ++i ) {
        //  const auto& vtx = (*vertices)[i];     
        //  const float weight = vtx.trackWeight(tau.leadPFChargedHadrCand()->trackRef());// check me -> Get track ref for charged hadron candidate
        //  if( weight > max_weight ) {
        //    max_weight = weight;
        //    vtx_index = i;
        //   }
        //      }

	      //now do vtx variable filling
	      //vtxIndex_ = vtx_index;
	      const reco::Vertex& vtx = (vtx_index == -1 ? (*vertices)[0] : (*vertices)[vtx_index]);
	      vtxX_ = vtx.x();
	      vtxY_ = vtx.y();
	      vtxZ_ = vtx.z();
	      
        break;
            }
          }
    jetTree->Fill(); 

    }
*/
   //jetTree for OrigTaus
   for(auto jet : Jets){
     tauPt_=-10;
     tauEta_=-10;
     tauMass_=-10;

     taupfTausDiscriminationByDecayModeFinding_=0;
     taupfTausDiscriminationByDecayModeFindingNewDMs_=0; 
     tauByLooseCombinedIsolationDeltaBetaCorr3Hits_=0;
     tauByMediumCombinedIsolationDeltaBetaCorr3Hits_=0;
     tauByTightCombinedIsolationDeltaBetaCorr3Hits_=0;
       
     tauCombinedIsolationDeltaBetaCorrRaw3Hits_=-10;  
     tauByIsolationMVArun2v1DBnewDMwLTraw_=-10;
     tauByIsolationMVArun2v1DBoldDMwLTraw_=-10;
     tauByIsolationMVArun2v1PWnewDMwLTraw_=-10;
     tauByIsolationMVArun2v1PWoldDMwLTraw_=-10;
     
     tauChargedIsoPtSum_=-10;
     tauNeutralIsoPtSum_=-10;
     tauPuCorrPtSum_=-10;

     taufootprintCorrection_=-10;
     tauphotonPtSumOutsideSignalCone_=-10;
     jetPt_=jet.pt();
     jetEta_=jet.eta();
     
     //nvtx_=-10;
     dm_=-10;
     dz_tt_=-10;
     jetTauMatch_=0;

     genJetMatch_ = 0;

     for (unsigned int iGenJet = 0; iGenJet < genJets->size() ; ++iGenJet){
       reco::GenJetRef genJet(genJets, iGenJet);
       if(genJet->pt() < 18 )continue;
       if (reco::deltaR(genJet->eta(),genJet->phi(),jet.eta(),jet.phi()) < 0.1){
          //&& genJet->pt() > 20 && abs(genJet->pt()-jet.pt())/abs(jet.pt())<=2 ???
         genJetMatch_ = 1;
         break;
        }
    }

     for(const pat::Tau &tau : *taus){   // Fix it was *tausOrg??
       for(auto cand : tau.isolationChargedHadrCands()){  ///dem
         if (abs(cand.charge())<0)continue;   ///dem >or < Fix
         if (cand.pt()<=0.5 or cand.dxy(vertices[tau_vertex_idxpf].position())>=0.1)continue; /////FIX ME dem
         if (cand.hasTrackDetails()){//dem
            tt = cand.pseudoTrack();//dem
            if (tt.normalizedChi2()>=100. or cand.numberOfHits()<3)continue;//dem
         }//dem
         dz_tt =cand.dz(vertices[tau_vertex_idxpf].position()); //dem
       } //dem

       if (reco::deltaR(tau.eta(),tau.phi(),jet.eta(),jet.phi()) < 0.3 && tau.tauID("decayModeFinding")>-1 && tau.tauID("chargedIsoPtSum")>-1){
	      jetTauMatch_ = 1;
	      tauPt_  =tau.pt();
	      tauEta_ =tau.eta();
        tauMass_=tau.mass();
	      dm_ =tau.decayMode();

        taupfTausDiscriminationByDecayModeFinding_=tau.tauID("decayModeFinding");
        taupfTausDiscriminationByDecayModeFindingNewDMs_=tau.tauID("decayModeFindingNewDMs");   
        tauByLooseCombinedIsolationDeltaBetaCorr3Hits_=tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
        tauByMediumCombinedIsolationDeltaBetaCorr3Hits_=tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
        tauByTightCombinedIsolationDeltaBetaCorr3Hits_=tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
        tauCombinedIsolationDeltaBetaCorrRaw3Hits_=tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        tauByIsolationMVArun2v1DBnewDMwLTraw_=tau.tauID("byIsolationMVArun2v1DBnewDMwLTraw");
        tauByIsolationMVArun2v1DBoldDMwLTraw_=tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
        tauByIsolationMVArun2v1PWnewDMwLTraw_=tau.tauID("byIsolationMVArun2v1PWnewDMwLTraw");
        tauByIsolationMVArun2v1PWoldDMwLTraw_=tau.tauID("byIsolationMVArun2v1PWoldDMwLTraw");

        tauChargedIsoPtSum_=tau.tauID("chargedIsoPtSum");
        tauNeutralIsoPtSum_=tau.tauID("neutralIsoPtSum");
        tauPuCorrPtSum_=tau.tauID("puCorrPtSum");
        taufootprintCorrection_=tau.tauID("footprintCorrection");
        tauphotonPtSumOutsideSignalCone_=tau.tauID("photonPtSumOutsideSignalCone");
        //get the matched vertex
       // int vtx_index = -1;// fix
        //float max_weight = 0.f;
        //for( unsigned i = 0; i < vertices->size(); ++i ) {
        //  const auto& vtx = (*vertices)[i];     
        //  const float weight = vtx.trackWeight(tau.leadPFChargedHadrCand()->trackRef());// check me -> Get track ref for charged hadron candidate
        //  if( weight > max_weight ) {
        //    max_weight = weight;
        //    vtx_index = i;
        //   }
        //      }

	      //now do vtx variable filling
	      //vtxIndex_ = vtx_index;
//	      const reco::Vertex& vtx = (vtx_index == -1 ? (*vertices)[0] : (*vertices)[vtx_index]);
//	      vtxX_ = vtx.x();
//	      vtxY_ = vtx.y();
//	      vtxZ_ = vtx.z();
	      
  //      break;
            }
          }
    jetTree->Fill(); //Fix it was jetOrgjetTree

    }

}

// ------------ method called once each job just before starting event loop  ------------
void 
phase2Taus::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
phase2Taus::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
phase2Taus::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
}


//Gets visible 4-momentum of a particle from list of daughters
reco::Candidate::LorentzVector phase2Taus::GetVisibleP4(std::vector<const reco::GenParticle*>& daughters){
  reco::Candidate::LorentzVector p4_vis(0,0,0,0);
  for(size_t i = 0; i < daughters.size(); ++i){
    if (!isNeutrino(daughters[i]) && daughters[i]->status() == 1){  //list should include intermediate daughters, so check for status = 1
      p4_vis += daughters[i]->p4();
    }
  }
  return p4_vis;
}


//Creates a vector of all (including intermediate) daughters for a given mother particle
void phase2Taus::findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters)
{
  unsigned numDaughters = mother->numberOfDaughters();
  if (numDaughters == 0) std::cout << " none ";
  for (unsigned iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
    const reco::GenParticle* daughter = mother->daughterRef(iDaughter).get();
    if (daughter->status() == 1){  //status = 1 is a final state daughter
      daughters.push_back(daughter); 
    }
    if (daughter->status() == 2){  //status = 2 is an intermediate daughter; will decay further
      daughters.push_back(daughter); 
      findDaughters(daughter, daughters);
    }
  }
}

bool phase2Taus::isNeutrino(const reco::Candidate* daughter)
{
  return (TMath::Abs(daughter->pdgId()) == 12 || TMath::Abs(daughter->pdgId()) == 14 || TMath::Abs(daughter->pdgId()) == 16 || TMath::Abs(daughter->pdgId()) == 18);
}


//define this as a plug-in
DEFINE_FWK_MODULE(phase2Taus);
