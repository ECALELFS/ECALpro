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

void drawICmap(string eosPath, string dirName, string iterNumber, string tagName) {

  string filename = "root://eoscms//eos/cms" + eosPath + dirName + "/" + iterNumber + "/" + tagName + "calibMap.root";

  TFile* f = TFile::Open(filename.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<filename<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  TH2F *mapEB = (TH2F*) f->Get("calibMap_EB");
  TH2F *h = NULL;
  h = (TH2F*) f->Get("calibMap_EEp");
  TH2F *mapEEp = (TH2F*) h->Clone();
  h = (TH2F*) f->Get("calibMap_EEm");
  TH2F *mapEEm = (TH2F*) h->Clone();

  if (!mapEB || !mapEEp || !mapEEm) {
    cout << "Error: could not get one or more histograms. End of programme" << endl;
    exit(EXIT_FAILURE);
  }

  TH2F *mapEB_new = new TH2F("mapEB_new","EB calib coefficients", 360, 0.5, 360.5, 171,-85.5,85.5 );

  Int_t nbinsX = mapEB->GetNbinsX(); // ieta
  Int_t nbinsY = mapEB->GetNbinsY(); // iphi

  for (Int_t i = 1; i <= nbinsX; i++) {

      for (Int_t j = 1; j <= nbinsY; j++) {
	
	mapEB_new->Fill(j,(i-86.0),mapEB->GetBinContent(i,j));

      }

  }

  string wwwPath = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/" + dirName + "/" + iterNumber + "/2DMaps/";
  string name = "";
  TPaletteAxis *palette = NULL;

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
  name = wwwPath + "Barrel/IC_calibMapEB";
  cEB->SaveAs((name + ".pdf").c_str());
  cEB->SaveAs((name + ".png").c_str());

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
  name = wwwPath + "Endcap/EEp/IC_calibMapEEp";
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
  name = wwwPath + "Endcap/EEm/IC_calibMapEEm";
  cEEm->SaveAs((name + ".pdf").c_str());
  cEEm->SaveAs((name + ".png").c_str());

  


}
