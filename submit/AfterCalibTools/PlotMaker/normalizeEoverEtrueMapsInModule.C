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

// macro to normalize EB maps to 1 for crystals inside modules without using gaps
// pass:
// 1) the name of the folder where to store the output plots
// 2) the input file name where the maps are stored (the file is probably on EOS, use root://eoscms//eos/cms/...)
 
void realNormalizeEoverEtrueMapsInModul(const string& outDir = "./",
					const string& inputFile = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/pi0Gun_MC_EoverEtrue_foldSM/iter_0/pi0Gun_MC_EoverEtrue_foldSM_calibMap.root",
					const Int_t nPhoton = 1, // 1 or 2
					const Double_t mapMin = 1.0, 
					const Double_t mapMax = -1			    
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

  TH2F *mapEB_SM = new TH2F("mapEB_SM",Form("cont. corr in SM (mean = 1 in each module) - #gamma%d",nPhoton), 85, 0.5, 85.5 , 20, 0.5, 20.5);
  vector<Float_t> meanInModuleNoGaps(4); // mean in a module without using crystals close to gaps

  // the map in the calibMap.root file has ieta on x axis
  // Int_t nbinsX = mapEB->GetNbinsX(); // ieta
  // Int_t nbinsY = mapEB->GetNbinsY(); // iphi
    
  // now we copy 1 SM


  for (Int_t j = 1; j <= 20; j++) {

    for (Int_t i = 1; i <= 85; i++) {
	
      mapEB_SM->SetBinContent(i,j,mapEB->GetBinContent(i+86,j));
      if (j != 1 && j != 20) {
	if      (i > 1  && i < 25) meanInModuleNoGaps[0] += mapEB->GetBinContent(i+86,j);
	else if (i > 26 && i < 45) meanInModuleNoGaps[1] += mapEB->GetBinContent(i+86,j);
	else if (i > 46 && i < 65) meanInModuleNoGaps[2] += mapEB->GetBinContent(i+86,j);
	else if (i > 66 && i < 85) meanInModuleNoGaps[3] += mapEB->GetBinContent(i+86,j);
      }

    }
      
  }

  // divide by number of crystals in module excluding borders
  meanInModuleNoGaps[0] /= (18*23);
  meanInModuleNoGaps[1] /= (18*18);
  meanInModuleNoGaps[2] /= (18*18);
  meanInModuleNoGaps[3] /= (18*18);

  for (Int_t j = 1; j <= 20; j++) {

    for (Int_t i = 1; i <= 85; i++) {
	
      if (i <= 25)      mapEB_SM->SetBinContent(i,j,mapEB_SM->GetBinContent(i,j)/meanInModuleNoGaps[0]);
      else if (i <= 45) mapEB_SM->SetBinContent(i,j,mapEB_SM->GetBinContent(i,j)/meanInModuleNoGaps[1]);
      else if (i <= 65) mapEB_SM->SetBinContent(i,j,mapEB_SM->GetBinContent(i,j)/meanInModuleNoGaps[2]);
      else              mapEB_SM->SetBinContent(i,j,mapEB_SM->GetBinContent(i,j)/meanInModuleNoGaps[3]);

    }

  }

  cout << "==================" << endl;
  cout << "Photon " << nPhoton << endl;
  for (UInt_t j = 0; j < meanInModuleNoGaps.size(); j++) {
    cout << "module " << j+1 << ": mean = " << meanInModuleNoGaps[j] << " (no gaps)" << endl; 
  }
  cout << "==================" << endl;

  string outHistName = (nPhoton == 1) ? "calibMap_EB" : "calibMap_EB_g2";
  TH2F *map_norm_allEB = new TH2F(outHistName.c_str(),Form("cont. corr in EB (mean = 1 in each module, no gaps) - #gamma%d",nPhoton), 171, -85-5, 85.5 , 360, 0.5, 360.5);
  // save in output file the new map on all EB
  string outputFile = outDir + "contCorrEoverEtrueNormTo1inEachModule.root";
  string outFileMode = (nPhoton == 1) ? "RECREATE" : "UPDATE";
  TFile* fout = TFile::Open(outputFile.c_str(),outFileMode.c_str());
  if (!fout || !fout->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error opening file \"" << outputFile << "\".\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }
  
  for (Int_t j = 1; j <= 360; j++) {

    for (Int_t i = 1; i <= 171; i++) {

      Int_t i_SM = 0;
      if (i > 86)      i_SM = i - 86; 
      else if (i < 86) i_SM = 86 - i;
      Int_t j_SM = (j - 1) % 20 + 1;
      map_norm_allEB->SetBinContent(i, j, mapEB_SM->GetBinContent(i_SM, j_SM));

    }

  }
  map_norm_allEB->Write();
  fout->Close();


  TCanvas* cEB = new TCanvas("cEB","");
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
  if (mapMin < mapMax) mapEB_SM->GetZaxis()->SetRangeUser(mapMin,mapMax);
  mapEB_SM->SetStats(0);
  gPad->Update();
  cEB->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_SingleSM_normalizedTo1inEachModule.pdf",outDir.c_str(),nPhoton));
  cEB->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_SingleSM_normalizedTo1inEachModule.png",outDir.c_str(),nPhoton));
  delete cEB;

  // all EB
  cEB = new TCanvas("cEB","");
  // cEB->SetLeftMargin(0.16);
  // cEB->SetRightMargin(0.20);
  cEB->cd();
  map_norm_allEB->Draw("COLZ");
  map_norm_allEB->GetXaxis()->SetTitle("i #eta");
  map_norm_allEB->GetXaxis()->SetTitleSize(0.06);
  map_norm_allEB->GetXaxis()->SetTitleOffset(0.7);
  map_norm_allEB->GetYaxis()->SetTitle("i #phi");
  map_norm_allEB->GetYaxis()->SetTitleSize(0.06);
  map_norm_allEB->GetYaxis()->SetTitleOffset(0.8);
  if (mapMin < mapMax) map_norm_allEB->GetZaxis()->SetRangeUser(mapMin,mapMax);
  map_norm_allEB->SetStats(0);
  gPad->Update();
  cEB->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_allSM_normalizedTo1inEachModule.pdf",outDir.c_str(),nPhoton));
  cEB->SaveAs(Form("%s/calibMap_EB_g%d_EoverEtrue_allSM_normalizedTo1inEachModule.png",outDir.c_str(),nPhoton));
  delete cEB;


  delete mapEB_SM;
  delete map_norm_allEB;

}


void normalizeEoverEtrueMapsInModule(const string& outDir = "./",
				     const string& inputFile = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/pi0Gun_MC_EoverEtrue_foldSM/iter_0/pi0Gun_MC_EoverEtrue_foldSM_calibMap.root",
				     const Double_t mapMin = 0,
				     const Double_t mapMax = -1) 
{

  system(Form("mkdir -p %s",outDir.c_str()));
  system(Form("cp /afs/cern.ch/user/m/mciprian/public/index.php %s",outDir.c_str()));

  realNormalizeEoverEtrueMapsInModul(outDir, inputFile, 1, mapMin, mapMax);
  realNormalizeEoverEtrueMapsInModul(outDir, inputFile, 2, mapMin, mapMax);
  

}
