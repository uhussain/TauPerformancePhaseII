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
#include "TH1F.h"
#include "TH2.h"
#include "THStack.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TTree.h"
#include "TPaveText.h"
#include "tdrstyle.C"
void applyPadStyle(TPad* pad1){
    pad1->SetFillColor(0);
    pad1->Draw();  pad1->cd();  pad1->SetLeftMargin(0.15);  pad1->SetBottomMargin(0.13); pad1->SetRightMargin(0.05);
    pad1->SetLogy();
    pad1->SetGrid(10,10);
      }
void plotFakeTauDistributions(){ 
  
  TString fileName = "/data/uhussain/TwoTausEff_July26_hadd/QCD_pre4_miniAOD_TwoTaus_pu0.root";  
  TString fileName1 = "/data/uhussain/TwoTausEff_July26_hadd/QCD_pre4_miniAOD_TwoTaus_pu140.root"; 
  
  //TString fileName2 = "/data/uhussain/TwoTausEff_July18_hadd/GluHtoTT_911_miniAOD_TwoTaus_pu0.root"; 
  TString treePath = "cutBased/jetOrgPFTaus"; 
  TString treePath2 = "cutBased/jetModFixedStripTaus";
  int bins = 50;
  double low  = 0;
  double high = 300;

  //Plotting variables
  //String variable = "(genTauPt-genTauPt)/(genTauPt)";
  //TString variable = "dm";
   TString variable = "tauPuCorrPtSum";
   //TString variable = "tauCombinedIsolationDeltaBetaCorrRaw3Hits";
  //TString variable2 = "(genTauPt-pfChHadrPt)/(genTauPt)";

  TString basic = "abs(tauEta)>0 && abs(tauEta)<2.3";
  //TString t3 = "vtxIndex!=-1 && isPV_zt3==1 && dmf==10";
  TString t3 = "genTauPt>30&&"+basic;
  //TString t4 = "dmf==1&&" + basic;
  //TString t5 = "dmf==10&&" + basic;
  //TString t6 = "tauIsPVt6==1";

  //Style
  TString xaxis = variable;

  Color_t color = TColor::GetColor("#283593");//dark blue color1
  Color_t color1 = TColor::GetColor("#F44336");//red color4 (Signal After fix)
  Color_t color2 = TColor::GetColor("#0288D1"); //green blue color2
  //Color_t colort3 = TColor::GetColor("#00FFFF");//Cyan
  //Color_t colort4 = TColor::GetColor("#F44336");//red color4
  //Color_t colort5 = TColor::GetColor("#0288D1"); //green blue color2
  Color_t color3 = TColor::GetColor("#FF00FF"); //magenta (Signal before fix)
  
  TString outFileName = "PuCorrPtSum-eta2.3-puall_Fakes";
  setTDRStyle();

  TFile *tauFile    = new TFile(fileName);

  if(!tauFile->IsOpen()||tauFile==0){
    std::cout<<"ERROR FILE "<< fileName<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
           }

  TFile *tauFile1    = new TFile(fileName1);

  if(!tauFile1->IsOpen()||tauFile1==0){
    std::cout<<"ERROR FILE "<< fileName1<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
           }
  //Signal file
  //TFile *tauFile2    = new TFile(fileName2);

  //if(!tauFile2->IsOpen()||tauFile2==0){
  //  std::cout<<"ERROR FILE "<< fileName2<<" NOT FOUND; EXITING"<<std::endl;
  //  exit(0);
  //         }
  TCanvas *Tcan= new TCanvas("Tcan","",100,20,600,600); Tcan->cd();  Tcan->SetFillColor(0);
  TPad* pad1 = new TPad("pad1","The pad",0,0,0.98,0.98);
  
  
  applyPadStyle(pad1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  TTree* tauTree = (TTree*)tauFile->Get(treePath);
  if(tauTree == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
         }
  //PU0 Signal after fix
  TTree* tauTree1 = (TTree*)tauFile->Get(treePath2);
  if(tauTree1 == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree1<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
         }
  TTree* tauTree2 = (TTree*)tauFile1->Get(treePath);
  if(tauTree2 == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree2<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
         }
  //PU200 Signal after fix
  TTree* tauTree3 = (TTree*)tauFile1->Get(treePath2);
  if(tauTree3 == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree3<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
         }
  //make the histograms
  //TH1F* tauPt,tauPt_t3,tauPt_t4,tauPt_t5,tauPt_t6;
  TH1F* PU0_beforefix = new TH1F("PU0_beforefix","PU0_beforefix",bins,low,high);  
  TH1F* PU0_afterfix = new TH1F("PU0_afterfix","PU0_afterfix",bins,low,high); 
  TH1F* PU200_beforefix = new TH1F("PU200_beforefix","PU200_beforefix",bins,low,high); 
  TH1F* PU200_afterfix = new TH1F("PU200_afterfix","PU200_afterfix",bins,low,high); 
  tauTree->Draw(variable+">>PU0_beforefix",basic);
  tauTree1->Draw(variable+">>PU0_afterfix",basic);
  tauTree2->Draw(variable+">>PU200_beforefix",basic);
  tauTree3->Draw(variable+">>PU200_afterfix",basic);
  
  PU0_beforefix->SetLineColor(color);
  PU0_afterfix->SetLineColor(color1);
  PU200_beforefix->SetLineColor(color2); 
  PU200_afterfix->SetLineColor(color3);

  //Normalize the histograms
  PU0_beforefix->Scale(1/(PU0_beforefix->Integral()));
  PU0_afterfix->Scale(1/(PU0_afterfix->Integral()));
  PU200_beforefix->Scale(1/(PU200_beforefix->Integral()));
  PU200_afterfix->Scale(1/(PU200_afterfix->Integral()));
  
  PU0_beforefix->GetXaxis()->SetTitle("PuCorrPtSum [GeV]");
  PU0_beforefix->GetXaxis()->SetTitleSize(0.03);

  PU0_beforefix->GetXaxis()->SetTitleOffset(1.4);
  //  DecayMode->GetYaxis()->SetTitle("Events");
//  DecayMode->GetYaxis()->SetTitleSize(0.03);
//  DecayMode->GetYaxis()->SetTitleOffset(1.9);

  //TVirtualPad *pad = Tcan->cd();
  //pad->SetLogy();
  PU0_beforefix->Draw("HIST");
  PU0_afterfix->Draw("HIST SAME");
  PU200_beforefix->Draw("HIST SAME");
  PU200_afterfix->Draw("HIST SAME");
  
  PU0_beforefix->SetLineWidth(2);  
  PU0_afterfix->SetLineWidth(2);  
  PU200_beforefix->SetLineWidth(2); 
  PU0_beforefix->SetMinimum(0.001);
  PU0_beforefix->SetMaximum(1.0);

  TString legLabel = "jet #rightarrow #tau_{h} fakes (QCDFlat)";
  TLegend *leg = new TLegend(.30, .70, .55, .90,legLabel);  
  leg->SetBorderSize(0);
  leg->SetShadowColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(0.03);
  leg->AddEntry(PU0_beforefix,"QCD 900_pre4(PU0)-BeforeFix");
  leg->AddEntry(PU0_afterfix,"QCD 900_pre4 (PU0)-AfterFix"); 
  leg->AddEntry(PU200_beforefix,"QCD 900_pre4 (PU140)-BeforeFix"); 
  leg->AddEntry(PU200_afterfix,"QCD 900_pre4 (PU140)-AfterFix");
  leg->Draw();
  
//  TLatex *texS = new TLatex(0.49,0.807173,"");
//  texS->SetNDC();
//  texS->SetTextFont(42);
//  texS->SetTextSize(0.040);
//  texS->Draw()TLatex *texS1 = new TLatex(0.16092,0.877173,"#bf{CMS} #it{Phase II Internal}");
  TLatex *texS1 = new TLatex(0.16092,0.877173,"#bf{CMS} #it{Internal}");
  texS1->SetTextFont(42);
  texS1->SetTextSize(0.040);
  texS1->Draw();

  Tcan->cd();
  Tcan->SaveAs("FakePlots/"+outFileName+".pdf");
}
