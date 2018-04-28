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

using namespace std;

// macro to plot maps of containment correction from E/Etrue in MC
// pass:
// 1) the name of the folder where to store the output plots
// 2) the input file name where the maps are stored (the file is probably on EOS, use root://eoscms//eos/cms/...)
 
void realDrawEoverEtrueMaps(const string& outDir = "",
			    const string& inputFile = "",
			    const string& inputFile2 = "",
			    const Int_t nPhoton = 1, // 1 or 2
			    const string& canvasSuffix = ""
			    )   
{

  TH1::SetDefaultSumw2();

  gStyle->SetPalette(55, 0);  // 55:raibow palette ; 57: kBird (blue to yellow) ; 107 kVisibleSpectrum ; 77 kDarkRainBow                                               
  gStyle->SetNumberContours(50); // default is 20 

  TFile* f = TFile::Open(inputFile.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error opening file \"" << inputFile << "\".\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }

  TH2F *mapEB = NULL;
  mapEB = (TH2F*) f->Get((nPhoton == 1) ? "calibMap_EB" : "calibMap_EB_g2");
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
  mapEB2 = (TH2F*) f2->Get((nPhoton == 1) ? "calibMap_EB" : "calibMap_EB_g2");
  if (!mapEB2) {
    cout << "Error: could not get EB histogram 2. End of programme" << endl;
    exit(EXIT_FAILURE);
  }
  mapEB2->SetDirectory(0);
  f2->Close();

  TH1F* hratioDistr = new TH1F("hratioDistr","",50, 0.975,1.025);

  TH2F *mapEB_SM = new TH2F("mapEB_SM",Form("containment correction in SM - #gamma%d",nPhoton), 85, 0.5, 85.5 , 20, 0.5, 20.5);
  TH2F *mapEB2_SM = new TH2F("mapEB2_SM",Form("containment correction in SM - #gamma%d",nPhoton), 85, 0.5, 85.5 , 20, 0.5, 20.5);
  // the map in the calibMap.root file has ieta on x axis
  // Int_t nbinsX = mapEB->GetNbinsX(); // ieta
  // Int_t nbinsY = mapEB->GetNbinsY(); // iphi
    
  // now we copy 1 SM
  for (Int_t i = 1; i <= 85; i++) {
      
    for (Int_t j = 1; j <= 20; j++) {
	
      mapEB_SM->Fill(i,j,mapEB->GetBinContent(i+86,j));
      mapEB2_SM->Fill(i,j,mapEB2->GetBinContent(i+86,j));
	
    }
      
  }

  TH2F* hRatio = (TH2F*) mapEB_SM->Clone("ratio");
  hRatio->Divide(mapEB2_SM);
  for (Int_t i = 0; i < hRatio->GetNbinsX(); i++) {
    for (Int_t j = 0; j < hRatio->GetNbinsY(); j++) {
      hratioDistr->Fill(hRatio->GetBinContent(i,j));
    }
  }

  //EB
  TCanvas *cEB = new TCanvas("cEB","");
  // cEB->SetLeftMargin(0.16);
  // cEB->SetRightMargin(0.20);
  cEB->SetRightMargin(0.15);
  cEB->cd();
  hRatio->Draw("COLZ");
  hRatio->GetXaxis()->SetTitle("i #eta");
  hRatio->GetXaxis()->SetTitleSize(0.06);
  hRatio->GetXaxis()->SetTitleOffset(0.7);
  hRatio->GetYaxis()->SetTitle("i #phi");
  hRatio->GetYaxis()->SetTitleSize(0.06);
  hRatio->GetYaxis()->SetTitleOffset(0.8);
  hRatio->GetZaxis()->SetRangeUser(0.975,(nPhoton == 1) ? 1.025 : 1.025);
  hRatio->SetStats(0);
  gPad->Update();
  cEB->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_%s.pdf",outDir.c_str(),nPhoton,canvasSuffix.c_str()));
  cEB->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_%s.png",outDir.c_str(),nPhoton,canvasSuffix.c_str()));
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
  cRatio1D->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_%s_1D.pdf",outDir.c_str(),nPhoton,canvasSuffix.c_str()));
  cRatio1D->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_%s_1D.png",outDir.c_str(),nPhoton,canvasSuffix.c_str()));
  delete cRatio1D;

  delete hratioDistr;
  delete mapEB_SM;
  delete mapEB2_SM;

}


void makeMapRatio(const string& outDir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/CC_EoverEtrue/ratio_CC/oldMCv4_newMC/",
		  const string& canvasSuffix = "ratioMC_oldOverNew",
		  const string& inputFile1 = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/pi0Gun_MC_EoverEtrue_foldSM_v4/iter_0/pi0Gun_MC_EoverEtrue_foldSM_v4_calibMap.root",
		  const string& inputFile2 = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/pi0Gun_MCV2_EoverEtrue_foldSM/iter_0/pi0Gun_MCV2_EoverEtrue_foldSM_calibMap.root") 
{

  system(Form("mkdir -p %s",outDir.c_str()));
  system(Form("cp /afs/cern.ch/user/m/mciprian/public/index.php %s",outDir.c_str()));

  realDrawEoverEtrueMaps(outDir, inputFile1, inputFile2, 1, canvasSuffix);
  realDrawEoverEtrueMaps(outDir, inputFile1, inputFile2, 2, canvasSuffix);
  

}
