#define calibAnaEcal_cxx

#include "calibAnaEcal.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TPaletteAxis.h>

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

#ifdef calibAnaEcal_cxx

calibAnaEcal::calibAnaEcal(TTree *tree) : calibAnaEcal_base(tree) {

  ////////////////////////////
  //initializing data members
  ///////////////////////////

  //////////////////////////
  // protected data members

  EBorEE = "";
  Pi0orEta = "";
  dirName = "";
  iterNumber = "";
  wwwPath = "";

  ////////////////////////////////
  // public data members

  normalizedS = 0.0;
  normalizedB = 0.0;

  Init(tree);  // calling Init() of this class
  

}

//===============================================

void calibAnaEcal::setHistograms() {

  hSignal = new TH2D("hSignal",Form("Signal in %s",EBorEE.c_str()),NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);
  hBackground = new TH2D("hBackground",Form("Background in %s",EBorEE.c_str()),NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);
  SoverB = new TH2D("SoverB",Form("S/B in %s",EBorEE.c_str()),NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);     
  SoverSqrtSplusB = new TH2D("SoverSqrtSplusB",Form("S/sqrt(S+B) in %s",EBorEE.c_str()),NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);
  SigmaMeanOverMean = new TH2D("SigmaMeanOverMean",Form("sigma(fit_mean)/fit_mean in %s",EBorEE.c_str()),NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap); //mean is the fit mean (should be the pi0 peak mass)  
  mean = new TH2D("mean",Form("fit_mean in %s",EBorEE.c_str()),NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);
  sigma = new TH2D("sigma",Form("fit_sigma in %s",EBorEE.c_str()),NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);
  chisquare = new TH2D("chisquare",Form("#Chi^{2} in %s",EBorEE.c_str()),NbinsX_2Dmap,lowerX_2Dmap,upperX_2Dmap,NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap);

  if (EBorEE == "EB")
    chisquare_vs_etaring = new TH2D("chisquare_vs_etaring",Form("#Chi^{2} in %s vs #eta-ring",EBorEE.c_str()),NbinsY_2Dmap,lowerY_2Dmap,upperY_2Dmap,100,0,100);
  else {
    chisquare_vs_etaring = new TH2D("chisquare_vs_etaring",Form("#Chi^{2} in %s vs #eta-ring",EBorEE.c_str()),NbinsX_etaProfile,lowerX_etaProfile,upperX_etaProfile,100,0,100);
  }

  th2dVector.push_back(hSignal);
  th2dVector.push_back(hBackground);
  th2dVector.push_back(SoverB);
  th2dVector.push_back(SoverSqrtSplusB);
  th2dVector.push_back(SigmaMeanOverMean);
  th2dVector.push_back(mean);
  th2dVector.push_back(sigma);
  th2dVector.push_back(chisquare);

  hSignal_etaProfile = new TProfile("hSignal_etaProfile",Form("Signal profile in %s",EBorEE.c_str()),NbinsX_etaProfile,lowerX_etaProfile,upperX_etaProfile);
  hBackground_etaProfile = new TProfile("hBackground_etaProfile",Form("Background profile in %s",EBorEE.c_str()),NbinsX_etaProfile,lowerX_etaProfile,upperX_etaProfile);
  SoverB_etaProfile = new TProfile("SoverB_etaProfile",Form("S/B profile in %s",EBorEE.c_str()),NbinsX_etaProfile,lowerX_etaProfile,upperX_etaProfile);
  SoverSqrtSplusB_etaProfile = new TProfile("SoverSqrtSplusB_etaProfile",Form("S/sqrt(S+B) profile in %s",EBorEE.c_str()),NbinsX_etaProfile,lowerX_etaProfile,upperX_etaProfile);
  SigmaMeanOverMean_etaProfile = new TProfile("SigmaMeanOverMean_etaProfile",Form("sigma(fit_mean)/fit_mean profile in %s",EBorEE.c_str()),NbinsX_etaProfile,lowerX_etaProfile,upperX_etaProfile);
  mean_etaProfile = new TProfile("mean_etaProfile",Form("fit_mean profile in %s",EBorEE.c_str()),NbinsX_etaProfile,lowerX_etaProfile,upperX_etaProfile);
  sigma_etaProfile = new TProfile("sigma_etaProfile",Form("fit_sigma profile in %s",EBorEE.c_str()),NbinsX_etaProfile,lowerX_etaProfile,upperX_etaProfile);
  chisquare_etaProfile = new TProfile("chisquare_etaProfile",Form("#Chi^{2} profile in %s",EBorEE.c_str()),NbinsX_etaProfile,lowerX_etaProfile,upperX_etaProfile);

  profileEtaVector.push_back(hSignal_etaProfile);
  profileEtaVector.push_back(hBackground_etaProfile);
  profileEtaVector.push_back(SoverB_etaProfile);
  profileEtaVector.push_back(SoverSqrtSplusB_etaProfile);
  profileEtaVector.push_back(SigmaMeanOverMean_etaProfile);
  profileEtaVector.push_back(mean_etaProfile);
  profileEtaVector.push_back(sigma_etaProfile);
  profileEtaVector.push_back(chisquare_etaProfile);

  profileYaxisTitle.push_back("events");
  profileYaxisTitle.push_back("events");
  profileYaxisTitle.push_back("S/B");
  profileYaxisTitle.push_back("S/#sqrt{S+B}");
  profileYaxisTitle.push_back("#sigma(mean)/mean");
  profileYaxisTitle.push_back("mean [GeV]");
  profileYaxisTitle.push_back("#sigma [GeV]");
  profileYaxisTitle.push_back("#Chi^{2}");

}


//===============================================

void calibAnaEcal::draw2Dmap(TH2D* hist2d) {

  gStyle->SetPalette(55, 0);  // 55:raibow palette  ; 57: kBird (from blue to yellow) ; 107: kVisibleSpectrum
  gStyle->SetNumberContours(50); // default is 20

  string canvasName(hist2d->GetName());
  canvasName = "c_" + canvasName;
  TCanvas *c = new TCanvas(canvasName.c_str(),canvasName.c_str());
  string name = wwwPath + hist2d->GetName() + "_" + EBorEE;  // name  (with path) of file to save canvas: EBorEE can be "EB" or "EEp" or "EEm" 
  c->SetRightMargin(0.18);

  hist2d->Draw("COLZ");
  if (EBorEE == "EB") {
    hist2d->GetXaxis()->SetTitle("i #phi");
    hist2d->GetYaxis()->SetTitle("i #eta");
  } else {
    hist2d->GetXaxis()->SetTitle("iX");
    hist2d->GetYaxis()->SetTitle("iY");
  } 
  hist2d->GetXaxis()->SetTitleSize(0.06);
  hist2d->GetXaxis()->SetTitleOffset(0.7);
  hist2d->GetYaxis()->SetTitleSize(0.06);
  hist2d->GetYaxis()->SetTitleOffset(0.7);
  hist2d->SetStats(0);
  hist2d->Draw("COLZ");
  // after drawing, fix the palette                                                                           

  gPad->Update();
  // TPaletteAxis *palette = (TPaletteAxis*)hist2d->GetListOfFunctions()->FindObject("palette");
  //// //cout << "hist2d = " << hist2d << "    palette = " << palette << endl;
  // if (!palette || palette == NULL) {
  //   cout << "Error in function draw2Dmap(): palette not found. ABORT" << endl;
  //   exit(EXIT_FAILURE);
  // }
  //// the following lines move the palette. Choose the values you need for the position.                                         
  // palette->SetX1NDC(0.91);
  // palette->SetX2NDC(0.94);
  // //palette->SetY1NDC(0.2);                                                                                                                                          
  // //palette->SetY2NDC(0.8);                                                                                                                                          
  gPad->Modified();
  gPad->Update();
  // end of palette fixes                                                                                                                                             
  c->SaveAs((name + ".pdf").c_str());
  c->SaveAs((name + ".png").c_str());

}

//===============================================

void calibAnaEcal::drawProfile(TProfile *profile, const string& yAxisName) {

  string canvasName(profile->GetName());
  canvasName = "c_" + canvasName;
  TCanvas *c = new TCanvas(canvasName.c_str(),canvasName.c_str());
  c->SetTickx(1);
  c->SetTicky(1);
  c->SetGrid();
  string name = wwwPath + profile->GetName() + "_" + EBorEE;  // name  (with path) of file to save canvas: EBorEE can be "EB" or "EEp" or "EEm" 

  profile->Draw("HE");
  if (EBorEE == "EB") {
    profile->GetXaxis()->SetTitle("#eta");
  } else {
    profile->GetXaxis()->SetTitle("#eta Ring");
  }
  profile->GetXaxis()->SetTitleSize(0.06);
  profile->GetXaxis()->SetTitleOffset(0.7);
  profile->GetYaxis()->SetTitle( yAxisName.c_str() );
  profile->GetYaxis()->SetTitleSize(0.055);
  profile->GetYaxis()->SetTitleOffset(0.8);
  profile->SetStats(0);
  profile->Draw("HE");
  c->SaveAs((name + ".pdf").c_str());
  c->SaveAs((name + ".png").c_str());


}

//===============================================

void calibAnaEcal::drawChisquare(TH2D* hist2d = NULL, const Bool_t drawProfileX = true) {
  
  gStyle->SetPalette(57, 0);  // 55:raibow palette  ; 57: kBird (blue to yellow, default one) ; 107: kVisibleSpectrum
  gStyle->SetNumberContours(50); // default is 20
  gStyle->SetOptStat(1000111110);


  string canvasName(hist2d->GetName());
  canvasName = "c_" + canvasName;
  TCanvas *c = new TCanvas(canvasName.c_str(),canvasName.c_str());
  string name = wwwPath + hist2d->GetName() + "_" + EBorEE;  // name  (with path) of file to save canvas: EBorEE can be "EB" or "EEp" or "EEm" 
  c->SetRightMargin(0.18);

  hist2d->Draw("COLZ");
  // hist2d->Draw("SCAT");
  // hist2d->SetMarkerStyle(3);
  // hist2d->SetMarkerColor(kRed+1);
  if (EBorEE == "EB") {
    hist2d->GetXaxis()->SetTitle("i #eta");
    hist2d->GetYaxis()->SetTitle("#Chi^{2}");
  } else {
    hist2d->GetXaxis()->SetTitle("#eta-ring");
    hist2d->GetYaxis()->SetTitle("#Chi^{2}");
  } 
  hist2d->GetXaxis()->SetTitleSize(0.06);
  hist2d->GetXaxis()->SetTitleOffset(0.7);
  hist2d->GetYaxis()->SetTitleSize(0.06);
  hist2d->GetYaxis()->SetTitleOffset(0.7);
  hist2d->Draw("COLZ");
  //hist2d->Draw("SCAT");
  // after drawing, fix the palette                                                                           

  gPad->Update();
  // TPaletteAxis *palette = (TPaletteAxis*)hist2d->GetListOfFunctions()->FindObject("palette");
  // //cout << "hist2d = " << hist2d << "    palette = " << palette << endl;
  // if (!palette || palette == NULL) {
  //   cout << "Error in function draw2Dmap(): palette not found. ABORT" << endl;
  //   exit(EXIT_FAILURE);
  // }
  // the following lines move the palette. Choose the values you need for the position.                                         
  // palette->SetX1NDC(0.91);
  // palette->SetX2NDC(0.94);
  // //palette->SetY1NDC(0.2);                                                                                                                                          
  // //palette->SetY2NDC(0.8);                                                                                                                                          
  gPad->Modified();
  gPad->Update();

  TProfile* h2DProfile = NULL;
  if (drawProfileX) {
    h2DProfile = hist2d->ProfileX(Form("%s_pfx",hist2d->GetName()));
    h2DProfile->SetMarkerColor(kBlack);
    h2DProfile->SetMarkerStyle(20);
    h2DProfile->SetMarkerSize(1);
    h2DProfile->Draw("EPsame");
  }
  TLegend leg(0.35,0.12,0.80,0.18);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  //leg.SetTextFont(62);
  if (drawProfileX) leg.AddEntry(h2DProfile,"mean #Chi^{2}","PLE");
  leg.Draw("same");


  // end of palette fixes                                                                                                                                             
  c->SaveAs((name + ".pdf").c_str());
  c->SaveAs((name + ".png").c_str());


}

//===============================================

void calibAnaEcal::setVerticalRangeInHisto() {

  for (UInt_t i = 0; i < th2dVector.size(); i++) {

    if (th2dMinZaxisVector[i] < 0.0) th2dMinZaxisVector[i] = 0.0;
    th2dVector[i]->SetMinimum(th2dMinZaxisVector[i]);
    //profileIetaVector[i]->SetMinimum(th2dMinZaxisVector[i]);                                                                                                         
    profileEtaVector[i]->SetMinimum(th2dMinZaxisVector[i]);

    // if the maximum choosen by user is bigger than default, don't do anything, otherwise set the user value as the maximum                                           
    if (th2dMaxZaxisVector[i] < th2dVector[i]->GetBinContent(th2dVector[i]->GetMaximumBin())) {

      th2dVector[i]->SetMaximum(th2dMaxZaxisVector[i]);
      profileEtaVector[i]->SetMaximum(th2dMaxZaxisVector[i]);

    }

  }


}

//===============================================

void calibAnaEcal::Init(TTree *tree) {

  calibAnaEcal_base::Init(tree);

}

void calibAnaEcal::Loop()
{
//   In a ROOT session, you can do:
//      root> .L calibAnaEcal.C
//      root> calibAnaEcal t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;                                                                                                                                   
    // if (jentry % 100000 == 0) cout << jentry << endl;

  }


}

#endif
