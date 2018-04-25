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
			    const string& inputFile = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/pi0Gun_MC_EoverEtrue_foldSM/iter_0/pi0Gun_MC_EoverEtrue_foldSM_calibMap.root",
			    const Int_t nPhoton = 1, // 1 or 2
			    const Double_t mapMin = 1.0, 
			    const Double_t mapMax = 1.12			    
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
  TH2F *h = NULL;

  mapEB = (TH2F*) f->Get((nPhoton == 1) ? "calibMap_EB" : "calibMap_EB_g2");
  if (!mapEB) {
    cout << "Error: could not get EB histogram. End of programme" << endl;
    exit(EXIT_FAILURE);
  }
  mapEB->SetDirectory(0);
  f->Close();

  TH2F *mapEB_SM = new TH2F("mapEB_SM",Form("containment correction in SM - #gamma%d",nPhoton), 85, 0.5, 85.5 , 20, 0.5, 20.5);
  TProfile * EB_ieta_profile = new TProfile("EB_ieta_profile",Form("containment correction in SM - i#eta profile - #gamma%d",nPhoton),85,0.5,85.5);
  TProfile * EB_iphi_profile = new TProfile("EB_iphi_profile",Form("containment correction in SM - i#phi profile - #gamma%d",nPhoton),20,0.5,20.5);

  // the map in the calibMap.root file has ieta on x axis
  // Int_t nbinsX = mapEB->GetNbinsX(); // ieta
  // Int_t nbinsY = mapEB->GetNbinsY(); // iphi
    
  // now we copy 1 SM

  for (Int_t i = 1; i <= 85; i++) {
      
    for (Int_t j = 1; j <= 20; j++) {
	
      mapEB_SM->Fill(i,j,mapEB->GetBinContent(i+86,j));
      EB_ieta_profile->Fill(i,mapEB->GetBinContent(i+86,j));
      EB_iphi_profile->Fill(j,mapEB->GetBinContent(i+86,j));
	
    }
      
  }

  //EB
  TCanvas *cEB = new TCanvas("cEB","");
  // cEB->SetLeftMargin(0.16);
  // cEB->SetRightMargin(0.20);
  cEB->cd();
  mapEB->Draw("COLZ");
  mapEB->GetXaxis()->SetTitle("i #eta");
  mapEB->GetXaxis()->SetTitleSize(0.06);
  mapEB->GetXaxis()->SetTitleOffset(0.7);
  mapEB->GetYaxis()->SetTitle("i #phi");
  mapEB->GetYaxis()->SetTitleSize(0.06);
  mapEB->GetYaxis()->SetTitleOffset(0.8);
  mapEB->GetZaxis()->SetRangeUser(mapMin,mapMax);
  mapEB->SetStats(0);
  gPad->Update();
  cEB->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_allSM.pdf",outDir.c_str(),nPhoton));
  cEB->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_allSM.png",outDir.c_str(),nPhoton));
  delete cEB;


  cEB = new TCanvas("cEB","");
  // cEB->SetLeftMargin(0.16);
  // cEB->SetRightMargin(0.20);
  cEB->cd();
  mapEB_SM->Draw("COLZ");
  mapEB_SM->GetXaxis()->SetTitle("i #eta");
  mapEB_SM->GetXaxis()->SetTitleSize(0.06);
  mapEB_SM->GetXaxis()->SetTitleOffset(0.7);
  mapEB_SM->GetYaxis()->SetTitle("i #phi");
  mapEB_SM->GetYaxis()->SetTitleSize(0.06);
  mapEB_SM->GetYaxis()->SetTitleOffset(0.8);
  mapEB_SM->GetZaxis()->SetRangeUser(mapMin,mapMax);
  mapEB_SM->SetStats(0);
  gPad->Update();
  cEB->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_SingleSM.pdf",outDir.c_str(),nPhoton));
  cEB->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_SingleSM.png",outDir.c_str(),nPhoton));
  delete cEB;

  
  TCanvas *cEB_ietaProfile = new TCanvas("cEB_ietaProfile","");
  EB_ieta_profile->Draw("HE");
  EB_ieta_profile->GetXaxis()->SetTitle("i #eta");
  EB_ieta_profile->GetXaxis()->SetTitleSize(0.06);
  EB_ieta_profile->GetXaxis()->SetTitleOffset(0.7);
  EB_ieta_profile->GetYaxis()->SetTitle("Containment correction");
  EB_ieta_profile->GetYaxis()->SetTitleSize(0.06);
  EB_ieta_profile->GetYaxis()->SetTitleOffset(0.8);
  // Double_t maxY = EB_ieta_profile->GetBinContent(EB_ieta_profile->GetMaximumBin());
  // Double_t scale_factor = 1.1;
  // Double_t minY = 999.9; // minimum would be 0, corresponding to ieta = 0; look for minimum excluding ieta = 0
  // for (Int_t ieta = 1; ieta<= 85; ieta++) {
  //   minY = (EB_ieta_profile->GetBinContent(ieta) < minY) ? EB_ieta_profile->GetBinContent(ieta) : minY;
  // }
  // Double_t offset = scale_factor * (maxY -minY); 
  //EB_ieta_profile->GetYaxis()->SetRangeUser(0.98, 1.1 * EB_ieta_profile->GetBinContent(EB_ieta_profile->GetMaximumBin()));
  //EB_ieta_profile->GetYaxis()->SetRangeUser(0.89,0.99);
  EB_ieta_profile->SetStats(0);
  gPad->Update();
  cEB_ietaProfile->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_SingleSM_ietaProfile.pdf",outDir.c_str(),nPhoton));
  cEB_ietaProfile->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_SingleSM_ietaProfile.png",outDir.c_str(),nPhoton));
  delete cEB_ietaProfile;

  TCanvas *cEB_iphiProfile = new TCanvas("cEB_iphiProfile","");
  EB_iphi_profile->Draw("HE");
  EB_iphi_profile->GetXaxis()->SetTitle("i #phi");
  EB_iphi_profile->GetXaxis()->SetTitleSize(0.06);
  EB_iphi_profile->GetXaxis()->SetTitleOffset(0.7);
  EB_iphi_profile->GetYaxis()->SetTitle("Containment correction");
  EB_iphi_profile->GetYaxis()->SetTitleSize(0.06);
  EB_iphi_profile->GetYaxis()->SetTitleOffset(0.8);
  // maxY = EB_iphi_profile->GetBinContent(EB_iphi_profile->GetMaximumBin());
  // minY = EB_iphi_profile->GetBinContent(EB_iphi_profile->GetMinimumBin()); 
  // offset = scale_factor * (maxY -minY); 
  //EB_iphi_profile->GetYaxis()->SetRangeUser(0.98, 1.1 * EB_iphi_profile->GetBinContent(EB_iphi_profile->GetMaximumBin()));
  //EB_iphi_profile->GetYaxis()->SetRangeUser(0.91,0.97);
  EB_iphi_profile->SetStats(0);
  gPad->Update();
  cEB_iphiProfile->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_SingleSM_iphiProfile.pdf",outDir.c_str(),nPhoton));
  cEB_iphiProfile->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_SingleSM_iphiProfile.png",outDir.c_str(),nPhoton));
  delete cEB_iphiProfile;


  delete mapEB_SM;
  delete EB_ieta_profile;
  delete EB_iphi_profile;

}


void drawEoverEtrueMaps(const string& outDir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/CC_EoverEtrue/pi0Gun_MCV2_EoverEtrue_foldSM_EoverEtrueCC_iter1/",
			const string& inputFile = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/pi0Gun_MCV2_EoverEtrue_foldSM_EoverEtrueCC_iter1/iter_0/pi0Gun_MCV2_EoverEtrue_foldSM_EoverEtrueCC_iter1_calibMap.root",
			const Double_t mapMin = 0.98,
			const Double_t mapMax = 1.02) 
{

  system(Form("mkdir -p %s",outDir.c_str()));
  system(Form("cp /afs/cern.ch/user/m/mciprian/public/index.php %s",outDir.c_str()));

  realDrawEoverEtrueMaps(outDir, inputFile, 1, mapMin, mapMax);
  realDrawEoverEtrueMaps(outDir, inputFile, 2, mapMin, mapMax);
  

}
