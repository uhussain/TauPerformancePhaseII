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
    //pad1->SetLogy();
    pad1->SetGrid(10,10);
      }
void plotTwoTauDistributions(){ 
  TString fileName = "/data/uhussain/TwoTausEff_July19_hadd/RelValTTbar_miniAOD_300.root";  
  TString fileName1 = "/data/uhussain/TwoTausEff_July19_hadd/RelValTTbar_miniAOD_3000.root"; 
  TString fileName2 = "/data/uhussain/TwoTausEff_July18_hadd/GluHtoTT_911_miniAOD_TwoTaus_pu0.root"; 
  TString treePath = "cutBased/OrgPFTaus"; 
  TString treePath2 = "cutBased/ModFixedStripTaus";
  int bins = 50;
  double low  = -1.0;
  double high = 1.0;

  //Plotting variables
  //String variable = "(genTauPt-genTauPt)/(genTauPt)";
  //TString variable = "dm";
   TString variable = "(genTauPt-tauPt)/(genTauPt)";
  //TString variable2 = "(genTauPt-pfChHadrPt)/(genTauPt)";

  TString basic = "abs(tauEta)>0 && abs(tauEta)<1.4";
  //TString t3 = "vtxIndex!=-1 && isPV_zt3==1 && dmf==10";
  TString t3 = "dm==0&&"+basic;
  //TString t4 = "dmf==1&&" + basic;
  //TString t5 = "dmf==10&&" + basic;
  //TString t6 = "tauIsPVt6==1";

  //Style
  TString xaxis = variable;

  Color_t color = TColor::GetColor("#283593");//dark blue color1
  Color_t color2 = TColor::GetColor("#F44336");//red color4 (Signal After fix)
  Color_t color1 = TColor::GetColor("#0288D1"); //green blue color2
  //Color_t colort3 = TColor::GetColor("#00FFFF");//Cyan
  //Color_t colort4 = TColor::GetColor("#F44336");//red color4
  //Color_t colort5 = TColor::GetColor("#0288D1"); //green blue color2
  Color_t color3 = TColor::GetColor("#FF00FF"); //magenta (Signal before fix)
  
  TString outFileName = "Res-Barrel_dm0-pu0_SignalVsRelVal";
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
  TFile *tauFile2    = new TFile(fileName2);

  if(!tauFile2->IsOpen()||tauFile2==0){
    std::cout<<"ERROR FILE "<< fileName2<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
           }
  TCanvas *Tcan= new TCanvas("Tcan","",100,20,600,600); Tcan->cd();  Tcan->SetFillColor(0);
  TPad* pad1 = new TPad("pad1","The pad",0,0,0.98,0.98);
  
  
  applyPadStyle(pad1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  TTree* tauTree = (TTree*)tauFile->Get(treePath2);
  if(tauTree == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
         }

  TTree* tauTree1 = (TTree*)tauFile1->Get(treePath2);
  if(tauTree == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree1<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
         }
  TTree* tauTree2 = (TTree*)tauFile2->Get(treePath2);
  if(tauTree == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree2<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
         }
  //Signal before fix
  TTree* tauTree3 = (TTree*)tauFile2->Get(treePath);
  if(tauTree == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree3<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
         }
  //make the histograms
  //TH1F* tauPt,tauPt_t3,tauPt_t4,tauPt_t5,tauPt_t6;
  TH1F* RelVal_300 = new TH1F("RelVal_300","RelVal_300",bins,low,high);  
  TH1F* RelVal_3000 = new TH1F("RelVal_3000","RelVal_3000",bins,low,high); 
  TH1F* Signal_AfterFix = new TH1F("Signal_AfterFix","Signal_AfterFix",bins,low,high); 
  TH1F* Signal_beforeFix = new TH1F("Signal_beforeFix","Signal_beforeFix",bins,low,high); 
  tauTree->Draw(variable+">>RelVal_300",t3);
  tauTree1->Draw(variable+">>RelVal_3000",t3);
  tauTree2->Draw(variable+">>Signal_AfterFix",t3);
  tauTree3->Draw(variable+">>Signal_beforeFix",t3);
  
  RelVal_300->SetLineColor(color);
  RelVal_3000->SetLineColor(color1);
  Signal_AfterFix->SetLineColor(color2); 
  Signal_beforeFix->SetLineColor(color3);

  //Normalize the histograms
  RelVal_300->Scale(1/(RelVal_300->Integral()));
  RelVal_3000->Scale(1/(RelVal_3000->Integral()));
  Signal_AfterFix->Scale(1/(Signal_AfterFix->Integral()));
  Signal_beforeFix->Scale(1/(Signal_beforeFix->Integral()));
  
  RelVal_300->GetXaxis()->SetTitle("1 - #tau P_{T}^{reco}/Gen #tau P_{T}^{vis}");
  RelVal_300->GetXaxis()->SetTitleSize(0.03);

  RelVal_300->GetXaxis()->SetTitleOffset(1.4);
  //  DecayMode->GetYaxis()->SetTitle("Events");
//  DecayMode->GetYaxis()->SetTitleSize(0.03);
//  DecayMode->GetYaxis()->SetTitleOffset(1.9);

  //TVirtualPad *pad = Tcan->cd();
  //pad->SetLogy();
  RelVal_300->Draw("HIST");
  RelVal_3000->Draw("HIST SAME");
  Signal_AfterFix->Draw("HIST SAME");
  Signal_beforeFix->Draw("HIST SAME");
  
  RelVal_300->SetLineWidth(2);  
  RelVal_3000->SetLineWidth(2);  
  Signal_AfterFix->SetLineWidth(2); 
  RelVal_300->SetMinimum(0.0);
  RelVal_300->SetMaximum(0.65);

  TString legLabel = "1prong0pi0 0 < |eta|<1.4";
  TLegend *leg = new TLegend(.30, .70, .55, .90,legLabel);  
  leg->SetBorderSize(0);
  leg->SetShadowColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(0.03);
  leg->AddEntry(RelVal_300,"ttbar relval (aging 300)-AfterFix");
  leg->AddEntry(RelVal_3000,"ttbar relval (aging 3000)-AfterFix"); 
  leg->AddEntry(Signal_AfterFix,"Signal GluHTT (PU0)-AfterFix"); 
  leg->AddEntry(Signal_beforeFix,"Signal GluHTT (PU0)-BeforeFix");
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
  Tcan->SaveAs("TwoTausPlots/"+outFileName+".pdf");
}
