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
// Author: Usama Hussain
// Code adapted originally from Isobel Ojalvo
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
  //edm::EDGetTokenT<pat::TauCollection> tauOrgToken_;
  std::string tauID_;
  edm::EDGetTokenT<std::vector <reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<std::vector<pat::Jet> > jetsAK4Token_;
  edm::EDGetTokenT<std::vector<reco::GenJet> > genJetsToken_;

  TTree* tree;
  //TTree* OrgTaustree;
  TTree* jetTree;
  //TTree* OrgTausjetTree;
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
//   std::vector<reco::Jet> Jets;
//   for (unsigned int iJet = 0; iJet < genJets->size() ; ++iJet){
     pat::JetRef jetCand(jetHandle, iJet);
     //const reco::Jet *jetCand = &(*genJets)[iJet];
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
     goodReco_=0;
     genTauMatch_=0;

     std::vector<const reco::GenParticle*> genTauDaughters;
     findDaughters(genTau, genTauDaughters);
     reco::Candidate::LorentzVector genTauVis = GetVisibleP4(genTauDaughters);
     genTauPt_  = (float) genTauVis.pt();
     genTauEta_ = (float) genTauVis.eta();

     //std::cout<<" pt "<<genTauPt_<<" eta "<<genTauEta_<<std::endl;
     for(const pat::Tau &tau : *taus){
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

            }
          }
    tree->Fill(); 
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

     for (unsigned int iGenJet = 0; iGenJet < genJets->size() ; ++iGenJet){
       reco::GenJetRef genJet(genJets, iGenJet);
       if(genJet->pt() < 18 )continue;
       //std::cout<<"genJetPt: "<<genJet->pt()<<std::endl;
       if (reco::deltaR(genJet->eta(),genJet->phi(),jet.eta(),jet.phi()) < 0.1){
         genJetMatch_ = 1;
         //std::cout<<"genJetMatched"<<std::endl;
         break;
       }
    }

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
            }
          }
    jetTree->Fill(); 

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
