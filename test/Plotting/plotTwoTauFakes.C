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
#include "tdrstyle.C"
#include "PlotStyles.C"
#include "CMS_lumi.C"
void applyPadStyle(TPad* pad1){
  pad1->SetFillColor(0);
  pad1->Draw();  pad1->cd();  pad1->SetLeftMargin(0.15);  pad1->SetBottomMargin(0.13); pad1->SetRightMargin(0.05);
  pad1->SetLogy();
  pad1->SetGrid(10,10);
}

void plotTwoTauFakes(){
   
  //TString fileName = "/data/uhussain/TwoTausEff_Aug8_hadd/DYJetsToLL_PU0.root"; 
  //TString fileName2 = "/data/uhussain/TwoTausEff_Aug8_hadd/DYJetsToLL_PU200.root";
  TString fileName = "/data/uhussain/Phase2Taus_Oct25_hadd/QCD_Flat_MiniAOD_PU0.root";  
  //TString fileName1 = "/data/uhussain/Phase2Taus_Oct25_hadd/QCD_Flat_MiniAOD_PU140.root"; 
  TString fileName2 = "/data/uhussain/Phase2Taus_Oct25_hadd/QCD_Flat_MiniAOD_PU200.root"; 
  //TString fileName2 = "/data/uhussain/TwoTausEff_Aug9_hadd/RelValTTbar_miniAOD_300.root";  
  //TString fileName3 = "/data/uhussain/TwoTausEff_Aug9_hadd/RelValTTbar_miniAOD_3000.root"; 
  //TString treePath = "cutBased/jetOrgPFTaus"; 
  TString treePath2 = "cutBased/jetModFixedStripTaus";
  //int bins = 30;
  //double low  = -3.0;
  //double high = 3.0;

  float binarray[]={0,20,40,60,80,100,120,160,200};
  int bins = sizeof(binarray)/sizeof(float) -1 ;

  //TString z1("abs(vtxZ) < 0.5"), 
  //  z2("abs(vtxZ) < 4.0 && abs(vtxZ) > 1.5"), 
  //  z3("abs(vtxZ) < 6.0 && abs(vtxZ) > 4.0"), 
  //  z4("abs(vtxZ) < 9.0 && abs(vtxZ) > 7.75");

  TString isoCut = "2";
  //Plotting Variables
  TString variable = "jetPt";
  TString GenCut= "jetPt> 20 && jetPt < 400 && genJetMatch > 0  && (jetEta) > 1.4 && (jetEta) < 3.0 && (dm!=5&&dm!=6 && dm > -1)";
  
  //TString GenCut1= "genTauPt > 22 && abs(genTauEta)> 2.3 && abs(genTauEta) <4.0 && (dmf!=5&&dmf!=6 && dmf > -1) && (dmf == 10) &&"+z3;
 
  //TString RecoCut= "tauPt > 20 && abs(tauEta)<2.3 && jetTauMatch==1 && taupfTausDiscriminationByDecayModeFinding==1 && tauChargedIsoPtSum<1.5 &&" + GenCut;

  TString RecoCut= "tauPt > 20 && (tauEta) > 1.4 && (tauEta)< 3.0 && taupfTausDiscriminationByDecayModeFinding==1 &&" + GenCut;

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
  //Style
  TString xaxis = "Gen Jet P_{T} (GeV)";
  int markerstyle = 20;

  Color_t color = TColor::GetColor("#283593");//dark blue color1
  Color_t color1 = TColor::GetColor("#F44336");//red color4 (Signal After fix)
  Color_t color3 = TColor::GetColor("#0288D1"); //green blue color2
  //Color_t colort3 = TColor::GetColor("#00FFFF");//Cyan
  //Color_t colort4 = TColor::GetColor("#F44336");//red color4
  //Color_t colort5 = TColor::GetColor("#0288D1"); //green blue color2
  Color_t color2 = TColor::GetColor("#FF00FF"); //magenta (Signal before fix)
  
  TString outFileName = "plot-Fakes-PT-HGCAL-puall_olddmf";

  //TString legLabel = "jet #rightarrow #tau_{h} (QCDFlat) oldDMF && chargedIso < 1.5 GeV";
  TString legLabel = "jet #rightarrow #tau_{h} (Fakes) from QCD_Pt-15To7000"; 
  
  //TString legLabel = "jet #rightarrow #tau_{h} (Fakes) from DYJets"; 
  //TString legLabel1 = "jet #rightarrow #tau_{h} (Fakes) from TTbar production (PU0)";
  
  TFile *tauFile    = new TFile(fileName);
  //PU0
  if(!tauFile->IsOpen()||tauFile==0){
    std::cout<<"ERROR FILE "<< fileName<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }
  //PU140
  //TFile *tauFile1    = new TFile(fileName1);

  //if(!tauFile1->IsOpen()||tauFile1==0){
  //  std::cout<<"ERROR FILE "<< fileName1<<" NOT FOUND; EXITING"<<std::endl;
  //  exit(0);
  //}

  TFile *tauFile2    = new TFile(fileName2);

  if(!tauFile2->IsOpen()||tauFile2==0){
    std::cout<<"ERROR FILE "<< fileName2<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }
  //TFile *tauFile3    = new TFile(fileName3);

  //if(!tauFile3->IsOpen()||tauFile3==0){
  //  std::cout<<"ERROR FILE "<< fileName3<<" NOT FOUND; EXITING"<<std::endl;
  //  exit(0);
  //}
  
  TCanvas *c1 = new TCanvas("c","c",800,800);
  //c1->SetGrid();
  //c1->SetLogy();
  setCanvasStyle(c1);
  c1->cd();
  
  //PU0 DYJets
  TTree* tauTree = (TTree*)tauFile->Get(treePath2);
  if(tauTree == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
         }
  //PU140 DYJets
  //TTree* tauTree1 = (TTree*)tauFile1->Get(treePath2);
  //if(tauTree1 == 0){
  //  std::cout<<"ERROR Tau Tree is "<< tauTree1<<" NOT FOUND; EXITING"<<std::endl;
  //  exit(0);
  //       } 
 //PU200 
  TTree* tauTree2 = (TTree*)tauFile2->Get(treePath2);
  if(tauTree2 == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree2<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }
  ////ttbar 3000
  //TTree* tauTree3 = (TTree*)tauFile3->Get(treePath2);
  //if(tauTree3 == 0){
  //  std::cout<<"ERROR Tau Tree is "<< tauTree3<<" NOT FOUND; EXITING"<<std::endl;
  //  exit(0);
  //}
  /// first
  TH1F* Denom;
  Denom = new TH1F("Denom","Denom",bins,binarray);
  Denom->Sumw2();
  tauTree->Draw(variable+">>+Denom",GenCut);

  TH1F* Num;
  Num = new TH1F("Num","Num",bins,binarray);
  tauTree->Draw(variable+">>+Num",RecoCut);
  std::cout<<"PU0 denominator: "<<Denom->Integral()<<std::endl;
  std::cout<<"PU0 Numerator: "<<Num->Integral()<<std::endl; 
  TGraphAsymmErrors *eff = new TGraphAsymmErrors();
  eff->Divide(Num,Denom,"cp");

  eff->SetMarkerStyle(markerstyle);
  eff->SetMarkerColor(color);
  eff->SetFillStyle(1001);
  eff->SetLineWidth((short)1.5);
  //Num->Divide(Denom);

  /// second
 // TH1F* Denom1;
 // Denom1 = new TH1F("Denom1","Denom1",bins,binarray);
 // Denom1->Sumw2();
 // tauTree1->Draw(variable+">>+Denom1",GenCut);

 // TH1F* Num1;
 // Num1 = new TH1F("Num1","Num1",bins,binarray);
 // tauTree1->Draw(variable+">>+Num1",RecoCut);

 // TGraphAsymmErrors *eff_1 = new TGraphAsymmErrors();
 //
 // //Num1->Divide(Denom1);
 // //for( unsigned ibin=1; ibin < Num1->GetNbinsX(); ++ibin )
 // // std::cout << ibin << " pass/total: " <<Num1->GetBinContent(ibin) << std::endl;
 // 
 // eff_1->Divide(Num1,Denom1,"cp");
 // 
 // //Num1->Divide(Denom1);
 // eff_1->SetMarkerStyle(markerstyle);
 // eff_1->SetMarkerColor(color1);

 // eff_1->SetFillStyle(1001);
 // eff_1->SetLineWidth((short)1.5);
  
  //Num1->Divide(Denom1);
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

  std::cout<<"PU200 denominator: "<<Denom2->Integral()<<std::endl;
  std::cout<<"PU200 Numerator: "<<Num2->Integral()<<std::endl;

  eff_2->Divide(Num2,Denom2,"cp");
  //Num2->Divide(Denom2);

  eff_2->SetMarkerStyle(markerstyle);
  eff_2->SetMarkerColor(color1);

  eff_2->SetFillStyle(1001);
  eff_2->SetLineWidth((short)1.5);
  //
  ////Num2->Divide(Denom2);

  //////
 /////fourth
  //TH1F* Denom3;
  //Denom3 = new TH1F("Denom3","Denom3",bins,binarray);
  //Denom3->Sumw2();
  //tauTree3->Draw(variable+">>+Denom3",GenCut);

  //TH1F* Num3;
  //Num3 = new TH1F("Num3","Num3",bins,binarray);
  //tauTree3->Draw(variable+">>+Num3",RecoCut);

  //TGraphAsymmErrors *eff_3 = new TGraphAsymmErrors();

  //eff_3->Divide(Num3,Denom3,"cp");
  ////Num3->Divide(Denom3);

  //eff_3->SetMarkerStyle(markerstyle);
  //eff_3->SetMarkerColor(color3);

  //eff_3->SetFillStyle(1001);
  //eff_3->SetLineWidth((short)1.5);
  //Num3->Divide(Denom3);

  gStyle->SetErrorX(0.5);

  TMultiGraph *mg = new TMultiGraph();
  //mg->SetTitle("Phase 2 Preliminary Studies");

  mg->Add(eff);
  //mg->Add(eff_1); 
  mg->Add(eff_2);
  //mg->Add(eff_3);

  mg->Draw("AP");
  
  mg->GetXaxis()->SetTitle(xaxis);
  mg->GetYaxis()->SetTitle("Tau Misidentification Probability");
  mg->GetYaxis()->SetTitleOffset(1.1);
  mg->SetMaximum(1.5);
  mg->SetMinimum(0);
  


  TLegend *leg = new TLegend(.30,.744,.67,.905,legLabel,"nbNDC");
  //setLegendStyles(leg,legLabel, 1);
  leg->SetBorderSize(0);
  leg->SetShadowColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(0.025);
  leg->AddEntry(eff,"PU0","PL"); 
  //leg->AddEntry(eff_1,"PU140","PL");
  leg->AddEntry(eff_2,"PU200","PL");
  //TLegend *leg1 = new TLegend(.30,.735,.67,.842,legLabel1,"nbNDC");
  ////setLegendStyles(leg,legLabel, 1);
  //leg1->SetBorderSize(0);
  //leg1->SetShadowColor(kWhite);
  //leg1->SetFillColor(kWhite);
  //leg1->SetTextSize(0.025); 
  //leg1->AddEntry(eff_2,"300 fb^{-1} lumi","PL");
  //leg1->AddEntry(eff_3,"3000 fb^{-1} lumi","PL");
  leg->Draw();
  //leg1->Draw();

  c1->cd();

  //Writing the lumi info
  CMS_lumi(c1, iPeriod, iPos );

  c1->Update();
  c1->RedrawAxis();
  c1->GetFrame()->Draw();
  c1->SaveAs("Oct25Plots/"+outFileName+".pdf");

}
