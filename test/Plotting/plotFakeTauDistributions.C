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
#include "PlotStyles.C"
#include "CMS_lumi.C"
void applyPadStyle(TPad* pad1){
    pad1->SetFillColor(0);
    pad1->Draw();  pad1->cd();  pad1->SetLeftMargin(0.15);  pad1->SetBottomMargin(0.13); pad1->SetRightMargin(0.05);
    pad1->SetLogy();
    pad1->SetGrid(10,10);
      }
void plotFakeTauDistributions(){ 
  
  
  TString fileName = "/data/uhussain/TwoTausEff_Aug8_hadd/DYJetsToLL_PU0.root"; 
  TString fileName1 = "/data/uhussain/TwoTausEff_Aug8_hadd/DYJetsToLL_PU200.root";
  TString fileName2 = "/data/uhussain/TwoTausEff_Aug9_hadd/RelValTTbar_miniAOD_300.root";  
  TString fileName3 = "/data/uhussain/TwoTausEff_Aug9_hadd/RelValTTbar_miniAOD_3000.root"; 
  //TString fileName = "/data/uhussain/TwoTausEff_July26_hadd/QCD_pre4_miniAOD_TwoTaus_pu0.root";  
  //TString fileName1 = "/data/uhussain/TwoTausEff_July26_hadd/QCD_pre4_miniAOD_TwoTaus_pu140.root";  
  //TString fileName2 = "/data/uhussain/TwoTausEff_July18_hadd/GluHtoTT_911_miniAOD_TwoTaus_pu0.root"; 
  //TString treePath = "cutBased/jetOrgPFTaus"; 
  TString treePath2 = "cutBased/jetModFixedStripTaus";
  int bins = 50;
  double low  = 0;
  double high = 2.5;

  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "       Phase-2 Simulation";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_sqrtS = "14 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos = 0; //0 is out of the frame

  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  //Plotting variables
  //String variable = "(genTauPt-genTauPt)/(genTauPt)";
  //TString variable = "tauChargedIsoPtSum";
   //TString variable = "tauPuCorrPtSum";
  TString variable = "tauMass";
  //TString variable2 = "(genTauPt-pfChHadrPt)/(genTauPt)";

  TString basic = "abs(tauEta)>0 && abs(tauEta)<1.5 && taupfTausDiscriminationByDecayModeFinding==1";
  //TString t3 = "vtxIndex!=-1 && isPV_zt3==1 && dmf==10";
  //TString t3 = "genTauPt>30&&"+basic;
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
  
  TString outFileName = "tauMass-eta1.5-puall_Fakes";
  
  TString legLabel = "jet #rightarrow #tau_{h} (Fakes) from DYJets"; 
  TString legLabel1 = "jet #rightarrow #tau_{h} (Fakes) from TTbar production (PU0)";
  
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

  TFile *tauFile2    = new TFile(fileName2);

  if(!tauFile2->IsOpen()||tauFile2==0){
    std::cout<<"ERROR FILE "<< fileName2<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }
  TFile *tauFile3    = new TFile(fileName3);

  if(!tauFile3->IsOpen()||tauFile3==0){
    std::cout<<"ERROR FILE "<< fileName3<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }


  TCanvas *c1 = new TCanvas("c","c",800,800);
  c1->SetGrid();
  //c1->SetLogy();
  setCanvasStyle(c1);
  c1->cd();
  
  //PU0 Signal after fix
  TTree* tauTree = (TTree*)tauFile->Get(treePath2);
  if(tauTree == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
         }
  //PU200 Signal after fix
  TTree* tauTree1 = (TTree*)tauFile1->Get(treePath2);
  if(tauTree1 == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree1<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
         }

 //ttbar 300 
  TTree* tauTree2 = (TTree*)tauFile2->Get(treePath2);
  if(tauTree2 == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree2<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }
  //ttbar 3000
  TTree* tauTree3 = (TTree*)tauFile3->Get(treePath2);
  if(tauTree3 == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree3<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }
  //make the histograms
  //TH1F* tauPt,tauPt_t3,tauPt_t4,tauPt_t5,tauPt_t6;
  TH1F* PU0_afterfix = new TH1F("PU0_afterfix","PU0_afterfix",bins,low,high); 
  TH1F* PU200_afterfix = new TH1F("PU200_afterfix","PU200_afterfix",bins,low,high);  
  TH1F* ttbar_300 = new TH1F("ttbar_300","ttbar_300",bins,low,high); 
  TH1F* ttbar_3000 = new TH1F("ttbar_3000","ttbar_3000",bins,low,high); 

  tauTree->Draw(variable+">>PU0_afterfix",basic);
  tauTree1->Draw(variable+">>PU200_afterfix",basic);
  
  tauTree2->Draw(variable+">>ttbar_300",basic);
  tauTree3->Draw(variable+">>ttbar_3000",basic);
  
  PU0_afterfix->SetLineColor(color);
  PU200_afterfix->SetLineColor(color1);
  ttbar_300->SetLineColor(color2);
  ttbar_3000->SetLineColor(color3);

  //Normalize the histograms
  PU0_afterfix->Scale(1/(PU0_afterfix->Integral()));
  PU200_afterfix->Scale(1/(PU200_afterfix->Integral()));
  
  ttbar_300->Scale(1/(ttbar_300->Integral()));
  ttbar_3000->Scale(1/(ttbar_3000->Integral()));

  PU0_afterfix->GetXaxis()->SetTitle("Tau Mass");  
  PU0_afterfix->GetYaxis()->SetTitle("a.u.");
  PU0_afterfix->GetXaxis()->SetTitleSize(0.03);

  PU0_afterfix->GetXaxis()->SetTitleOffset(1.4);
  PU0_afterfix->Draw("HIST");
  PU200_afterfix->Draw("HIST SAME");
  
  ttbar_300->Draw("HIST SAME");
  ttbar_3000->Draw("HIST SAME");
  
  PU0_afterfix->SetLineWidth(2);  
  PU200_afterfix->SetLineWidth(2);
  ttbar_300->SetLineWidth(2); 
  ttbar_3000->SetLineWidth(2);

  PU0_afterfix->SetMinimum(0);
  PU0_afterfix->SetMaximum(0.3);

  TLegend *leg = new TLegend(.30,.844,.67,.925,legLabel,"nbNDC");
  //setLegendStyles(leg,legLabel, 1);
  leg->SetBorderSize(0);
  leg->SetShadowColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(0.025);
  leg->AddEntry(PU0_afterfix,"PU0","PL");
  leg->AddEntry(PU200_afterfix,"PU200","PL");
  TLegend *leg1 = new TLegend(.30,.735,.67,.842,legLabel1,"nbNDC");
  //setLegendStyles(leg,legLabel, 1);
  leg1->SetBorderSize(0);
  leg1->SetShadowColor(kWhite);
  leg1->SetFillColor(kWhite);
  leg1->SetTextSize(0.025); 
  leg1->AddEntry(ttbar_300,"300 fb^{-1} lumi","PL");
  leg1->AddEntry(ttbar_3000,"3000 fb^{-1} lumi","PL");
  leg->Draw();
  leg1->Draw();

  c1->cd();

  //Writing the lumi info
  CMS_lumi(c1, iPeriod, iPos );

  c1->Update();
  c1->RedrawAxis();
  c1->GetFrame()->Draw();
  c1->SaveAs("TDRPlots/"+outFileName+".pdf");
}
