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

void drawICmap(const string& wwwPath = "",
	       const string& eosPath = "", 
	       const string& dirName = "", 
	       const string& iterNumber = "", 
	       const string& tagName = "",
	       const string& ECALdetToSkip = "") 
{

  gStyle->SetPalette(107, 0);  // 1:raibow palette ; 107 kVisibleSpectrum                                               
  gStyle->SetNumberContours(50); // default is 20 

  string filename = "root://eoscms//eos/cms" + eosPath + dirName + "/" + iterNumber + "/" + tagName + "calibMap.root";

  TFile* f = TFile::Open(filename.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<filename<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  TH2F *mapEB = NULL;
  TH2F *h = NULL;
  TH2F *mapEEp = NULL;
  TH2F *mapEEm = NULL;

  if (ECALdetToSkip != "EB") {
    mapEB = (TH2F*) f->Get("calibMap_EB");
    if (!mapEB) {
      cout << "Error: could not get EB histogram. End of programme" << endl;
      exit(EXIT_FAILURE);
    }
  }
  if (ECALdetToSkip != "EE") {
    h = (TH2F*) f->Get("calibMap_EEp");
    mapEEp = (TH2F*) h->Clone();
    if (!mapEEp) {
      cout << "Error: could not get EE+ histogram. End of programme" << endl;
      exit(EXIT_FAILURE);
    }
    h = (TH2F*) f->Get("calibMap_EEm");
    mapEEm = (TH2F*) h->Clone();
    if (!mapEEm) {
      cout << "Error: could not get EE- histogram. End of programme" << endl;
      exit(EXIT_FAILURE);
    }
  }    


  TH2F *mapEB_new = new TH2F("mapEB_new","EB calib coefficients", 360, 0.5, 360.5, 171,-85.5,85.5 );

  // profile along ieta. ieta goea from -85 to 85, exluding 0, for a total of 170 non empty bins (they are 171 including ieta = 0 which is actually empty)
  // in the profile, ieta = 30 is the bin with x from 29.5 to 30.5
  // simila logic for profile along iphi
  TProfile * EB_ieta_profile = new TProfile("EB_ieta_profile","EB calib coefficients - i#eta profile",171,-85.5,85.5);
  TProfile * EB_iphi_profile = new TProfile("EB_iphi_profile","EB calib coefficients - i#phi profile",360,0.5,360.5);


  if (ECALdetToSkip != "EB") {

    Int_t nbinsX = mapEB->GetNbinsX(); // ieta
    Int_t nbinsY = mapEB->GetNbinsY(); // iphi
    
    for (Int_t i = 1; i <= nbinsX; i++) {
      
      for (Int_t j = 1; j <= nbinsY; j++) {
	
	mapEB_new->Fill(j,(i-86.0),mapEB->GetBinContent(i,j));
	EB_ieta_profile->Fill((i-86.0),mapEB->GetBinContent(i,j));
	EB_iphi_profile->Fill(j,mapEB->GetBinContent(i,j));
	
      }
      
    }
    
  }

  string wwwAllPath = wwwPath + dirName + "/" + iterNumber + "/2DMaps/";
  string name = "";
  TPaletteAxis *palette = NULL;

  if (ECALdetToSkip != "EB") {

    //EB
    TCanvas *cEB = new TCanvas("cEB","IC map EB");
    mapEB_new->Draw("COLZ");
    mapEB_new->GetXaxis()->SetTitle("i #phi");
    mapEB_new->GetXaxis()->SetTitleSize(0.06);
    mapEB_new->GetXaxis()->SetTitleOffset(0.7);
    mapEB_new->GetYaxis()->SetTitle("i #eta");
    mapEB_new->GetYaxis()->SetTitleSize(0.06);
    mapEB_new->GetYaxis()->SetTitleOffset(0.8);
    mapEB_new->GetZaxis()->SetRangeUser(0.9,1.1);
    mapEB_new->SetStats(0);
    gPad->Update();
    palette = (TPaletteAxis*)mapEB_new->GetListOfFunctions()->FindObject("palette");
    // the following lines move the palette. Choose the values you need for the position.                                                                              
    palette->SetX1NDC(0.91);
    palette->SetX2NDC(0.94);
    gPad->Modified();
    gPad->Update();
    // end of palette fixes                                                                                                                                             
    name = wwwAllPath + "Barrel/IC_calibMapEB";
    cEB->SaveAs((name + ".pdf").c_str());
    cEB->SaveAs((name + ".png").c_str());

    TCanvas *cEB_ietaProfile = new TCanvas("cEB_ietaProfile","IC map EB - i#eta profile");
    EB_ieta_profile->Draw("HIST");
    EB_ieta_profile->GetXaxis()->SetTitle("i #eta");
    EB_ieta_profile->GetXaxis()->SetTitleSize(0.06);
    EB_ieta_profile->GetXaxis()->SetTitleOffset(0.7);
    EB_ieta_profile->GetYaxis()->SetTitle("IC");
    EB_ieta_profile->GetYaxis()->SetTitleSize(0.06);
    EB_ieta_profile->GetYaxis()->SetTitleOffset(0.8);
    Double_t maxY = EB_ieta_profile->GetBinContent(EB_ieta_profile->GetMaximumBin());
    Double_t scale_factor = 1.1;
    Double_t minY = 999.9; // minimum would be 0, corresponding to ieta = 0; look for minimum excluding ieta = 0
    for (Int_t ieta = -85; ieta<= 85; ieta++) {
      if (ieta == 0) continue;
      minY = (EB_ieta_profile->GetBinContent(ieta+86) < minY) ? EB_ieta_profile->GetBinContent(ieta+86) : minY;
    }
    Double_t offset = scale_factor * (maxY -minY); 
    EB_ieta_profile->GetYaxis()->SetRangeUser(minY - offset, maxY + offset);
    //EB_ieta_profile->GetYaxis()->SetRangeUser(0.89,0.99);
    EB_ieta_profile->SetStats(0);
    gPad->Update();
    name = wwwAllPath + "Barrel/IC_calibMapEB_ietaProfile";
    cEB_ietaProfile->SaveAs((name + ".pdf").c_str());
    cEB_ietaProfile->SaveAs((name + ".png").c_str());


    TCanvas *cEB_iphiProfile = new TCanvas("cEB_iphiProfile","IC map EB - i#phi profile");
    EB_iphi_profile->Draw("HIST");
    EB_iphi_profile->GetXaxis()->SetTitle("i #phi");
    EB_iphi_profile->GetXaxis()->SetTitleSize(0.06);
    EB_iphi_profile->GetXaxis()->SetTitleOffset(0.7);
    EB_iphi_profile->GetYaxis()->SetTitle("IC");
    EB_iphi_profile->GetYaxis()->SetTitleSize(0.06);
    EB_iphi_profile->GetYaxis()->SetTitleOffset(0.8);
    maxY = EB_iphi_profile->GetBinContent(EB_iphi_profile->GetMaximumBin());
    minY = EB_iphi_profile->GetBinContent(EB_iphi_profile->GetMinimumBin()); 
    offset = scale_factor * (maxY -minY); 
    EB_iphi_profile->GetYaxis()->SetRangeUser(minY - offset, maxY + offset);
    //EB_iphi_profile->GetYaxis()->SetRangeUser(0.91,0.97);
    EB_iphi_profile->SetStats(0);
    gPad->Update();
    name = wwwAllPath + "Barrel/IC_calibMapEB_iphiProfile";
    cEB_iphiProfile->SaveAs((name + ".pdf").c_str());
    cEB_iphiProfile->SaveAs((name + ".png").c_str());

  }

  if (ECALdetToSkip != "EE") {

    //EE+
    TCanvas *cEEp = new TCanvas("cEEp","IC map EE+");
    mapEEp->Draw("COLZ");
    mapEEp->GetXaxis()->SetTitle("iX");
    mapEEp->GetXaxis()->SetTitleSize(0.06);
    mapEEp->GetXaxis()->SetTitleOffset(0.7);
    mapEEp->GetYaxis()->SetTitle("iY");
    mapEEp->GetYaxis()->SetTitleSize(0.06);
    mapEEp->GetYaxis()->SetTitleOffset(0.8);
    mapEEp->GetZaxis()->SetRangeUser(0.75,1.25);
    mapEEp->SetStats(0);
    gPad->Update();
    palette = (TPaletteAxis*)mapEEp->GetListOfFunctions()->FindObject("palette");
    // the following lines move the palette. Choose the values you need for the position.                    
    palette->SetX1NDC(0.91);
    palette->SetX2NDC(0.94);
    gPad->Modified();
    gPad->Update();
    // end of palette fixes                                    
    name = wwwAllPath + "Endcap/EEp/IC_calibMapEEp";
    cEEp->SaveAs((name + ".pdf").c_str());
    cEEp->SaveAs((name + ".png").c_str());

    //EE-
    TCanvas *cEEm = new TCanvas("cEEm","IC map EE-");
    mapEEm->Draw("COLZ");
    mapEEm->GetXaxis()->SetTitle("iX");
    mapEEm->GetXaxis()->SetTitleSize(0.06);
    mapEEm->GetXaxis()->SetTitleOffset(0.7);
    mapEEm->GetYaxis()->SetTitle("iY");
    mapEEm->GetYaxis()->SetTitleSize(0.06);
    mapEEm->GetYaxis()->SetTitleOffset(0.8);
    mapEEm->GetZaxis()->SetRangeUser(0.75,1.25);
    mapEEm->SetStats(0);
    gPad->Update();
    palette = (TPaletteAxis*)mapEEm->GetListOfFunctions()->FindObject("palette");
    // the following lines move the palette. Choose the values you need for the position.                    
    palette->SetX1NDC(0.91);
    palette->SetX2NDC(0.94);
    gPad->Modified();
    gPad->Update();
    // end of palette fixes                                    
    name = wwwAllPath + "Endcap/EEm/IC_calibMapEEm";
    cEEm->SaveAs((name + ".pdf").c_str());
    cEEm->SaveAs((name + ".png").c_str());

  }  


}
