#include "tdrstyle.C"
#include "PlotStyles.C"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include <math.h>
#include <TEfficiency.h>
#include <TMath.h>
#include <TFormula.h>
#include <iostream>
#include <string>
#include <iostream>
#include <cmath>
#include "TLegend.h"
#include "TPad.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TF1.h"
#include "THStack.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TTree.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"

void applyPadStyle(TPad* pad1){
  pad1->SetFillColor(0);
  pad1->Draw();  pad1->cd();  pad1->SetLeftMargin(0.15);  pad1->SetBottomMargin(0.13); pad1->SetRightMargin(0.05);
  //pad1->SetGrid();
  pad1->SetGrid(10,10);
}

void plotTwoTauClopperEfficiency(){
  TString fileName = " /data/uhussain/TwoTausEff_July26_hadd/Ztt_pre4_miniAOD_TwoTaus_pu0.root";  
  TString fileName1 = "/data/uhussain/TwoTausEff_July26_hadd/Ztt_pre4_miniAOD_TwoTaus_pu140.root"; 
  //TString fileName2 = "/data/uhussain/TauTiming_July_hadd/Ztt_pre4_pu200.root"; 
  TString treePath = "cutBased/OrgPFTaus"; 
  TString treePath2 = "cutBased/ModFixedStripTaus";
  //int bins = 10;
  //double low  = 0;
  //double high = 200;

  float binarray[]={0,10,20,30,40,50,60,70,80,90,100,150,200};
  int bins = sizeof(binarray)/sizeof(float) -1 ;

  TString z1("abs(vtxZ) < 0.5"), 
    z2("abs(vtxZ) < 4.0 && abs(vtxZ) > 1.5"), 
    z3("abs(vtxZ) < 6.0 && abs(vtxZ) > 4.0"), 
    z4("abs(vtxZ) < 9.0 && abs(vtxZ) > 7.75");

  TString isoCut = "2";
  //Plotting Variables
  TString variable = "genTauPt";
  TString GenCut= "genTauPt> 20 && abs(genTauEta) < 2.3 && (dm!=5&&dm!=6 && dm > -1)";
  
  //TString GenCut1= "genTauPt > 22 && abs(genTauEta)> 2.3 && abs(genTauEta) <4.0 && (dmf!=5&&dmf!=6 && dmf > -1) && (dmf == 10) &&"+z3;
 
  //TString RecoCut= "tauPt > 20 && abs(tauEta)<2.3 && taupfTausDiscriminationByDecayModeFinding==1 &&" + GenCut;
  
  TString RecoCut= "tauPt > 20 && abs(tauEta)<2.3 && taupfTausDiscriminationByDecayModeFinding==1 && tauCombinedIsolationDeltaBetaCorrRaw3Hits<1.5 &&" + GenCut;
  //TString RecoCut1= "tauPt>22&&abs(tauEta)>2.3&& abs(genTauEta)<4.0&&" + GenCut1;
  //TString RecoCut1= "tauPt>22&&PFCharged<"+isoCut+"&&" + GenCut;
  //TString RecoCut2= "tauPt>22&&PFCharged<"+isoCut+"&&" + GenCut;


  //Style
  TString xaxis = "Gen #tau P_{T}^{vis}";
  int markerstyle = 20;
  Color_t color = TColor::GetColor("#283593");//dark blue color1
  Color_t color1 = TColor::GetColor("#F44336");//red color4
  Color_t color2 = TColor::GetColor("#0288D1"); //green blue color2 
  Color_t color3 = TColor::GetColor("#FF00FF"); //magenta (Signal before fix)
  TString outFileName = "plot-ZTT-PT-pu0-140_olddmf_dBcorrIso_";
  TString legLabel = "Z #rightarrow #tau #tau oldDMF && combinedIso (dBcorr) < 1.5 GeV";
  setTDRStyle();

  TH1F *basehist = new TH1F("basehist","",100,0,2.5);

  basehist->SetStats(false);

  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  TFile *tauFile    = new TFile(fileName);
  //PU0
  if(!tauFile->IsOpen()||tauFile==0){
    std::cout<<"ERROR FILE "<< fileName<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }
  //PU140
  TFile *tauFile1    = new TFile(fileName1);

  if(!tauFile1->IsOpen()||tauFile1==0){
    std::cout<<"ERROR FILE "<< fileName1<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }

  //TFile *tauFile2    = new TFile(fileName1);

  //if(!tauFile2->IsOpen()||tauFile2==0){
  //  std::cout<<"ERROR FILE "<< fileName1<<" NOT FOUND; EXITING"<<std::endl;
  //  exit(0);
  //}
  TCanvas *Tcan= new TCanvas("Tcan","",100,20,600,600); Tcan->cd();  Tcan->SetFillColor(0);
  TPad* pad1 = new TPad("pad1","The pad",0,0,0.98,0.98);

  applyPadStyle(pad1);
  gStyle->SetOptFit(0);
  gStyle->SetEndErrorSize(0);
  
  //PU0 Before Fix
  TTree* tauTree = (TTree*)tauFile->Get(treePath);
  if(tauTree == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }
  //PU After Fix
  TTree* tauTree1 = (TTree*)tauFile->Get(treePath2);
  if(tauTree1 == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree1<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }
  //PU140 Before Fix
  TTree* tauTree2 = (TTree*)tauFile1->Get(treePath);
  if(tauTree2 == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree2<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }
  //PU140 After Fix
  TTree* tauTree3 = (TTree*)tauFile1->Get(treePath2);
  if(tauTree3 == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree3<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }
  /// first
  TH1F* Denom;
  Denom = new TH1F("Denom","Denom",bins,binarray);
  Denom->Sumw2();
  tauTree->Draw(variable+">>+Denom",GenCut);

  TH1F* Num;
  Num = new TH1F("Num","Num",bins,binarray);
  tauTree->Draw(variable+">>+Num",RecoCut);

  TGraphAsymmErrors *eff = new TGraphAsymmErrors();
  eff->Divide(Num,Denom,"cp");

  eff->SetMarkerStyle(markerstyle);
  eff->SetMarkerColor(color);
  eff->SetFillStyle(1001);
  eff->SetLineWidth((short)1.5);
  //Num->Divide(Denom);

  /// second
  TH1F* Denom1;
  Denom1 = new TH1F("Denom1","Denom1",bins,binarray);
  Denom1->Sumw2();
  tauTree1->Draw(variable+">>+Denom1",GenCut);

  TH1F* Num1;
  Num1 = new TH1F("Num1","Num1",bins,binarray);
  tauTree1->Draw(variable+">>+Num1",RecoCut);

  TGraphAsymmErrors *eff_1 = new TGraphAsymmErrors();
 
  //Num1->Divide(Denom1);
  //for( unsigned ibin=1; ibin < Num1->GetNbinsX(); ++ibin )
  // std::cout << ibin << " pass/total: " <<Num1->GetBinContent(ibin) << std::endl;
  
  eff_1->Divide(Num1,Denom1,"cp");
  
  //Num1->Divide(Denom1);
  eff_1->SetMarkerStyle(markerstyle);
  eff_1->SetMarkerColor(color1);

  eff_1->SetFillStyle(1001);
  eff_1->SetLineWidth((short)1.5);
  ////
  /// third
  TH1F* Denom2;
  Denom2 = new TH1F("Denom2","Denom2",bins,binarray);
  Denom2->Sumw2();
  tauTree2->Draw(variable+">>+Denom2",GenCut);

  TH1F* Num2;
  Num2 = new TH1F("Num2","Num2",bins,binarray);
  tauTree2->Draw(variable+">>+Num2",RecoCut);

  TGraphAsymmErrors *eff_2 = new TGraphAsymmErrors();

  eff_2->Divide(Num2,Denom2,"cp");
  //Num2->Divide(Denom2);

  eff_2->SetMarkerStyle(markerstyle);
  eff_2->SetMarkerColor(color2);

  eff_2->SetFillStyle(1001);
  eff_2->SetLineWidth((short)1.5);
  ////
  ///fourth

  TH1F* Denom3;
  Denom3 = new TH1F("Denom3","Denom3",bins,binarray);
  Denom3->Sumw2();
  tauTree3->Draw(variable+">>+Denom3",GenCut);

  TH1F* Num3;
  Num3 = new TH1F("Num3","Num3",bins,binarray);
  tauTree3->Draw(variable+">>+Num3",RecoCut);

  TGraphAsymmErrors *eff_3 = new TGraphAsymmErrors();

  eff_3->Divide(Num3,Denom3,"cp");
  //Num3->Divide(Denom3);

  eff_3->SetMarkerStyle(markerstyle);
  eff_3->SetMarkerColor(color3);

  eff_3->SetFillStyle(1001);
  eff_3->SetLineWidth((short)1.5);
  
  gStyle->SetErrorX(0.5);
  
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("Phase 2 Preliminary Studies");

  mg->Add(eff);
  mg->Add(eff_1);
  mg->Add(eff_2);
  mg->Add(eff_3);
  //Denom->Draw("P same");

  mg->Draw("AP");
  
  mg->GetXaxis()->SetTitle(xaxis);
  mg->GetYaxis()->SetTitle("Efficiency");
  mg->GetYaxis()->SetTitleOffset(1.1);
  mg->SetMaximum(1.5);
  mg->SetMinimum(0);
  

  TLegend *leg = new TLegend(.30,.744,.67,.925,legLabel,"nbNDC");
  //setLegendStyles(leg,legLabel, 1);
  leg->SetBorderSize(0);
  leg->SetShadowColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(0.025);  
  leg->AddEntry(eff,"Ztt 900_pre4 (PU0)-BeforeFix","PL");
  leg->AddEntry(eff_1,"Ztt 900_pre4 (PU0)-AfterFix","PL");
  leg->AddEntry(eff_2,"Ztt 900_pre4 (PU140)-BeforeFix","PL");
  leg->AddEntry(eff_3,"Ztt 900_pre4 (PU140)-AfterFix","PL"); 
  //leg->AddEntry(Num,"0PU numerator","PL");
  //leg->AddEntry(Denom,"0PU denom","PL");
  leg->Draw();


  Tcan->cd();
  Tcan->SaveAs("FakePlots/"+outFileName+".pdf");
}
