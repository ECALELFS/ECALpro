#include <TROOT.h>   
#include <TAttFill.h>
#include <TAxis.h>  
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>  
#include <TFile.h>
#include <TFitResult.h>  
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TVirtualFitter.h>

#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h                                
#include <cstdio>
#include <cmath>
#include <sstream> //to use ostringstream to convert numbers to string in c++          

#include "./ICfunctions.h"

using namespace std;

// macro to plot of ratio between two calibration maps in EB 
// pass:
// 1) the name of the folder where to store the output plots
// 2) the input file name where the maps are stored (the file is probably on EOS, use root://eoscms//eos/cms/...)
 
void realDrawMapRatio(const string& outDir = "",
		      const string& inputFile = "",
		      const string& inputFile2 = "",
		      const string& canvasSuffix = "",
		      const string& mapName1 = "",
		      const string& mapName2 = "",
		      const Double_t mapMin = 0.95,
		      const Double_t mapMax = 1.05
		      )   
{

  TH1::SetDefaultSumw2();

  gStyle->SetPalette(55, 0);  // 55:raibow palette ; 57: kBird (blue to yellow) ; 107 kVisibleSpectrum ; 77 kDarkRainBow          
  gStyle->SetNumberContours(100); // default is 20 

  TFile* f = TFile::Open(inputFile.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error opening file \"" << inputFile << "\".\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }

  TH2F *mapEB = NULL;
  mapEB = (TH2F*) f->Get(mapName1.c_str());
  if (!mapEB) {
    cout << "Error: could not get EB histogram. End of programme" << endl;
    exit(EXIT_FAILURE);
  }
  mapEB->SetDirectory(0);
  f->Close();


  TFile* f2 = TFile::Open(inputFile2.c_str(),"READ");
  if (!f2 || !f2->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error opening file \"" << inputFile2 << "\".\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }

  TH2F *mapEB2 = NULL;
  mapEB2 = (TH2F*) f2->Get(mapName2.c_str());
  if (!mapEB2) {
    cout << "Error: could not get EB histogram 2. End of programme" << endl;
    exit(EXIT_FAILURE);
  }
  mapEB2->SetDirectory(0);
  f2->Close();

  TH2F* mapEB_new = new TH2F("mapEB_new","", 360, 0.5, 360.5, 171, -85.5, 85.5);
  TH2F* mapEB2_new = new TH2F("mapEB2_new","", 360, 0.5, 360.5, 171, -85.5, 85.5);

  Double_t mapEB_binContent = 0.0;
  Int_t bin = 0;
  
  for (Int_t iphi = 1; iphi <= 360; iphi++) {

    for (Int_t ieta = -85; ieta <= 85; ieta++) {

      bin = mapEB->FindFixBin(iphi,ieta);
      mapEB_binContent = mapEB->GetBinContent( bin );
      mapEB_new->SetBinContent(bin, mapEB_binContent);
      mapEB_binContent = mapEB2->GetBinContent( bin );
      mapEB2_new->SetBinContent(bin, mapEB_binContent);

    }

  }

  TH1F* hratioDistr = new TH1F("hratioDistr","",50, 0.975,1.025);

  TH2F* hRatio = (TH2F*) mapEB_new->Clone("ratio");
  divideEBmap(hRatio,mapEB_new,mapEB2_new,true,0);
  //hRatio->Divide(mapEB2_new);
  for (Int_t i = 1; i <= hRatio->GetNbinsX(); i++) {
    for (Int_t j = 1; j <= hRatio->GetNbinsY(); j++) {
      if (j != 86) hratioDistr->Fill(hRatio->GetBinContent(i,j));
    }
  }

  TProfile* profileEB_phi_final = new TProfile("profileEB_phi_final","i#phi profile of IC ratio in EB",360, 0.5, 360.5);
  makeICprofileIphiFromMap(profileEB_phi_final, hRatio, true, true, false);
  drawDistribution(profileEB_phi_final, "i#phi", "mean IC", Form("calibMap_EB_ratio_%s_iphiProfile",canvasSuffix.c_str()), outDir, 0.5, 360.5, 700, 500, true);


  //EB
  Int_t xsizeCanvas = 1200;
  Int_t ysizeCanvas = 1.0 * xsizeCanvas * mapEB_new->GetNbinsY() / mapEB_new->GetNbinsX() + 0.1 *xsizeCanvas;

  TCanvas *cEB = new TCanvas("cEB","",xsizeCanvas,ysizeCanvas);
  // cEB->SetLeftMargin(0.16);
  // cEB->SetRightMargin(0.20);
  cEB->SetRightMargin(0.14);
  cEB->cd();
  hRatio->Draw("COLZ");
  hRatio->GetXaxis()->SetTitle("i #phi");
  hRatio->GetXaxis()->SetTitleSize(0.06);
  hRatio->GetXaxis()->SetTitleOffset(0.7);
  hRatio->GetYaxis()->SetTitle("i #eta");
  hRatio->GetYaxis()->SetTitleSize(0.06);
  hRatio->GetYaxis()->SetTitleOffset(0.8);
  hRatio->GetZaxis()->SetRangeUser(mapMin, mapMax);
  hRatio->SetStats(0);
  gPad->Update();
  cEB->SaveAs(Form("%s/calibMap_EB_ratio_%s.pdf",outDir.c_str(),canvasSuffix.c_str()));
  cEB->SaveAs(Form("%s/calibMap_EB_ratio_%s.png",outDir.c_str(),canvasSuffix.c_str()));
  delete cEB;

  //EB
  TCanvas *cRatio1D = new TCanvas("cRatio1D","");
  // cRatio1D->SetLeftMargin(0.16);
  // cRatio1D->SetRightMargin(0.20);
  cRatio1D->cd();
  hratioDistr->Draw("HIST");
  hratioDistr->GetXaxis()->SetTitle("ratio");
  hratioDistr->GetXaxis()->SetTitleSize(0.06);
  hratioDistr->GetXaxis()->SetTitleOffset(0.7);
  hratioDistr->GetYaxis()->SetTitle("Events");
  hratioDistr->GetYaxis()->SetTitleSize(0.06);
  hratioDistr->GetYaxis()->SetTitleOffset(0.8);
  //hratioDistr->SetStats(0);
  gPad->Update();
  cRatio1D->SaveAs(Form("%s/calibMap_EB_ratio_%s_1D.pdf",outDir.c_str(),canvasSuffix.c_str()));
  cRatio1D->SaveAs(Form("%s/calibMap_EB_ratio_%s_1D.png",outDir.c_str(),canvasSuffix.c_str()));
  delete cRatio1D;

  delete hratioDistr;

}


void makeICratio(//const string& outDir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/AlCaP0_Run2017_F_CCiter0/iter_6/2DMaps/ratio/",
		 const string& outDir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/ratioIC/AlCaP0_Run2018B_ext1_fromIter3_iter4__Over__AlCaP0_Run2018A_iter6/",
		 const string& canvasSuffix = "ratioIC",
		 const string& inputFile1 = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/AlCaP0_Run2018B_ext1_fromIter3/iter_4/2DMaps/ICmaps/IC_work/calibMap_EB_divided_foldSMafterNorm1eachModulePlusMinusSeparate_norm1etaRing.root",
		 const string& inputFile2 = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/AlCaP0_Run2018A/iter_6/2DMaps/ICmaps/IC_work/calibMap_EB_divided_foldSMafterNorm1eachModulePlusMinusSeparate_norm1etaRing.root",
		 const string& mapName1 = "mapEB_original_Over_norm1eachModuleFoldSMallEB_plusMinusSeparate_norm1etaRing",
		 const string& mapName2 = "mapEB_original_Over_norm1eachModuleFoldSMallEB_plusMinusSeparate_norm1etaRing",
		 const Double_t mapMin = 0.98,
		 const Double_t mapMax = 1.02
		 ) 
{

  system(Form("mkdir -p %s",outDir.c_str()));
  system(Form("cp /afs/cern.ch/user/m/mciprian/public/index.php %s",outDir.c_str()));

  realDrawMapRatio(outDir, inputFile1, inputFile2, canvasSuffix, mapName1, mapName2,mapMin,mapMax);
  

}
