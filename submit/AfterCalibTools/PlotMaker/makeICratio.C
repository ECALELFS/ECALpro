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
#include "./CMS_lumi.h"

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
		      const Double_t mapMax = 1.05,
		      const string& canvasLabel1 = "IC set 1",
		      const string& canvasLabel2 = "IC set 2"
		      )   
{

  TH1::SetDefaultSumw2();
  TH1::StatOverflows(kTRUE);

  gStyle->SetPalette(55, 0);  // 55:raibow palette ; 57: kBird (blue to yellow) ; 107 kVisibleSpectrum ; 77 kDarkRainBow          
  gStyle->SetNumberContours(101); // default is 20 

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
  //Int_t bin = 0;
  Int_t binx = 0;
  Int_t biny = 0;  

  for (Int_t iphi = 1; iphi <= 360; iphi++) {

    for (Int_t ieta = -85; ieta <= 85; ieta++) {

      //bin = mapEB->FindFixBin(iphi,ieta);
      binx = mapEB->GetXaxis()->FindFixBin(iphi);
      biny = mapEB->GetYaxis()->FindFixBin(ieta);
      //mapEB_binContent = mapEB->GetBinContent( bin );
      mapEB_binContent = mapEB->GetBinContent( binx, biny );
      mapEB_new->SetBinContent(binx, biny, mapEB_binContent);
      //mapEB_binContent = mapEB2->GetBinContent( bin );
      mapEB_binContent = mapEB2->GetBinContent( binx, biny );
      mapEB2_new->SetBinContent(binx, biny, mapEB_binContent);

    }

  }

  TH1F* hratioDistr = new TH1F("hratioDistr","",51, 0.975,1.025);

  TH2F* hRatio = (TH2F*) mapEB_new->Clone("ratio");
  Double_t ratioVal = 0.0;
  divideEBmap(hRatio,mapEB_new,mapEB2_new,true,0);
  //hRatio->Divide(mapEB2_new);
  for (Int_t i = 1; i <= hRatio->GetNbinsX(); i++) {
    for (Int_t j = 1; j <= hRatio->GetNbinsY(); j++) {
      if (j != 86) {
	ratioVal = hRatio->GetBinContent(i,j);
	if (ratioVal > 0.00001) hratioDistr->Fill(ratioVal);
      }
    }
  }

  TProfile* profileEB_phi_final = new TProfile("profileEB_phi_final","i#phi profile of IC ratio in EB",360, 0.5, 360.5);
  makeICprofileIphiFromMap(profileEB_phi_final, hRatio, true, true, false);
  drawDistribution(profileEB_phi_final, "i#phi", "mean IC", Form("calibMap_EB_ratio_%s_iphiProfile",canvasSuffix.c_str()), outDir, 0.5, 360.5, 700, 500, true);


  //EB
  Int_t xsizeCanvas = 1200;
  Int_t ysizeCanvas = 1.0 * xsizeCanvas * mapEB_new->GetNbinsY() / mapEB_new->GetNbinsX() + 0.1 *xsizeCanvas;

  TCanvas *cEB = new TCanvas("cEB","",xsizeCanvas,ysizeCanvas);
  cEB->SetLeftMargin(0.08);
  // cEB->SetRightMargin(0.20);
  cEB->SetRightMargin(0.14);
  cEB->cd();
  hRatio->Draw("COLZ");
  hRatio->GetXaxis()->SetTitle("i #phi");
  hRatio->GetXaxis()->SetTitleSize(0.06);
  hRatio->GetXaxis()->SetTitleOffset(0.7);
  hRatio->GetYaxis()->SetTitle("i #eta");
  hRatio->GetYaxis()->SetTitleSize(0.06);
  hRatio->GetYaxis()->SetTitleOffset(0.65);
  hRatio->GetZaxis()->SetRangeUser(mapMin, mapMax);
  hRatio->GetZaxis()->SetTitle("IC ratio");
  hRatio->GetZaxis()->SetTitleSize(0.05);
  hRatio->GetZaxis()->SetTitleOffset(0.9);
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
  cRatio1D->SetTickx(1);
  cRatio1D->SetTicky(1);
  cRatio1D->cd();
  cRatio1D->SetGridx(1);
  cRatio1D->SetGridy(1);
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

  vector<TH1*> dispersionIC_EBmod;
  vector<TH1*> dispersionIC2_EBmod;
  for (Int_t i = 1; i <=4; i++) {
    dispersionIC_EBmod.push_back(  new TH1D(Form("dispersionIC_EBmod%d",i), Form("IC dispersion in EB Module %d",i),51, 0.95,1.05) );
    dispersionIC2_EBmod.push_back( new TH1D(Form("dispersionIC2_EBmod%d",i),Form("IC dispersion in EB Module %d",i),51, 0.95,1.05) );
    // do not use Under/Overflow bins for these plots
    dispersionIC_EBmod.back()->StatOverflows(0);
    dispersionIC2_EBmod.back()->StatOverflows(0);
  }

  Double_t imodule = 0;
  Int_t fabsieta = -1;
  for (Int_t i = 1; i <= mapEB->GetNbinsX(); i++) {
    for (Int_t j = 1; j <= mapEB->GetNbinsY(); j++) {
      if (j != 86) {	
	fabsieta = fabs(j - 86); 
	if      (fabsieta <= 25) imodule = 1;
	else if (fabsieta <= 45) imodule = 2;
	else if (fabsieta <= 65) imodule = 3;
	else                     imodule = 4;
	mapEB_binContent = mapEB->GetBinContent(i,j);
	if (isGoodIC(mapEB_binContent)) {
	  dispersionIC_EBmod[imodule-1]->Fill(mapEB_binContent);
	}
	mapEB_binContent = mapEB2->GetBinContent(i,j);
	if (isGoodIC(mapEB_binContent)) {
	  dispersionIC2_EBmod[imodule-1]->Fill(mapEB_binContent);
	}
      }
    }
  }

  TCanvas* canvas = new TCanvas("canvas","",1800,600);
  TLatex* tex = new TLatex();
  setTDRStyle();
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  Double_t clm = 0.05;
  Double_t crm = 0.02;
  canvas->SetBottomMargin(0.1);
  canvas->SetRightMargin(crm);
  canvas->SetLeftMargin(clm);
  //canvas->Divide(4,1,0.0,0.0);
  Double_t maxy = 0.0;
  for (Int_t i = 0; i < 4; i++) {
    Double_t tmpmax = std::max(dispersionIC_EBmod[i]->GetBinContent(dispersionIC_EBmod[i]->GetMaximumBin()), 
			       dispersionIC2_EBmod[i]->GetBinContent(dispersionIC2_EBmod[i]->GetMaximumBin()));
    maxy = std::max(maxy, tmpmax);
  }

  Double_t pw = (1. - clm - crm) / 4.;

  TLatex* hinfo = new TLatex();
  hinfo->SetTextSize(0.035);
  hinfo->SetTextFont(42);

  vector<TPad*> pads;
  for (Int_t i = 0; i < 4; i++) {
    // Double_t padlowx = (i == 0) ? 0 : i*pw + clm;
    // pads.push_back( new TPad(Form("pad%d",i),Form("pad%d",i), padlowx, 0., pw*(i+1) + clm, 0.97) );
    // pads[i]->SetLeftMargin((i == 0) ? 0.18 : 0.0);
    // pads[i]->SetRightMargin(0.0);

    pads.push_back( new TPad(Form("pad%d",i),Form("pad%d",i), 0, 0., 1, 0.97) );
    pads[i]->SetLeftMargin(clm + pw*i);
    pads[i]->SetRightMargin(crm + pw*(3-i));
    pads[i]->Draw();
    pads[i]->cd();
    //pads[i]->SetTopMargin(0.7);
    pads[i]->SetFillColor(0);
    pads[i]->SetGridy(1);
    pads[i]->SetGridx(1);
    pads[i]->SetFillStyle(0);
    pads[i]->SetTickx(1);
    pads[i]->SetTicky(1);

    dispersionIC_EBmod[i]->SetStats(0);
    dispersionIC2_EBmod[i]->SetStats(0);    
    dispersionIC_EBmod[i]->SetLineColor(kBlack);
    dispersionIC_EBmod[i]->GetYaxis()->SetRangeUser(0, 1.2*maxy);
    dispersionIC_EBmod[i]->Draw("HIST");
    dispersionIC_EBmod[i]->GetXaxis()->SetTitle("IC value");
    dispersionIC_EBmod[i]->GetXaxis()->SetTitleSize(0.06);
    dispersionIC_EBmod[i]->GetXaxis()->SetLabelSize(0.035);
    dispersionIC_EBmod[i]->GetXaxis()->SetTitleOffset(0.7);
    if (i == 0) {
      dispersionIC_EBmod[i]->GetYaxis()->SetTitle("Events");
      dispersionIC_EBmod[i]->GetYaxis()->SetTitleSize(0.055);
      dispersionIC_EBmod[i]->GetYaxis()->SetLabelSize(0.04);
      dispersionIC_EBmod[i]->GetYaxis()->SetTitleOffset(0.48);
    } else {
      dispersionIC_EBmod[i]->GetYaxis()->SetTitleSize(0.0);
      dispersionIC_EBmod[i]->GetYaxis()->SetLabelSize(0.0);
    }
    dispersionIC_EBmod[i]->GetYaxis()->SetTickSize(0.01);
    dispersionIC_EBmod[i]->GetXaxis()->SetTickSize(0.05);
    dispersionIC2_EBmod[i]->SetLineColor(kRed+2);
    dispersionIC2_EBmod[i]->Draw("HIST SAME");
    pads[i]->RedrawAxis("sameaxis");
    pads[i]->Update();
    canvas->cd();
    tex->SetTextSize(0.04);
    tex->SetTextFont(62);
    tex->SetTextColor(kGreen+2);
    tex->DrawLatex(i*pw + clm + 0.02 ,0.83, Form("Module %d",i+1));
    Double_t pass = 0.04;
    Double_t yh = 0.75;
    hinfo->SetTextColor(kBlack);
    hinfo->DrawLatex(i*pw + clm + 0.02, yh,         Form("Mean %.3f",dispersionIC_EBmod[i]->GetMean()));
    hinfo->DrawLatex(i*pw + clm + 0.02, yh-pass,    Form("RMS  %.3f",dispersionIC_EBmod[i]->GetStdDev()));
    hinfo->DrawLatex(i*pw + clm + 0.02, yh-2.*pass, Form("Integral %.0f",dispersionIC_EBmod[i]->Integral()));    
    hinfo->DrawLatex(i*pw + clm + 0.02, yh-3.*pass, Form("Overflow %.0f",dispersionIC_EBmod[i]->GetBinContent(dispersionIC_EBmod[i]->GetNbinsX()+1)));
    hinfo->DrawLatex(i*pw + clm + 0.02, yh-4.*pass, Form("Underflow  %.0f",dispersionIC_EBmod[i]->GetBinContent(0)));
    hinfo->SetTextColor(kRed+2);
    hinfo->DrawLatex(i*pw + clm + 0.03 + pw/2., yh,         Form("Mean %.3f",dispersionIC2_EBmod[i]->GetMean()));
    hinfo->DrawLatex(i*pw + clm + 0.03 + pw/2., yh-pass,    Form("RMS  %.3f",dispersionIC2_EBmod[i]->GetStdDev()));
    hinfo->DrawLatex(i*pw + clm + 0.03 + pw/2., yh-2.*pass, Form("Integral %.0f",dispersionIC2_EBmod[i]->Integral()));
    hinfo->DrawLatex(i*pw + clm + 0.03 + pw/2., yh-3.*pass, Form("Overflow %.0f",dispersionIC2_EBmod[i]->GetBinContent(dispersionIC_EBmod[i]->GetNbinsX()+1)));
    hinfo->DrawLatex(i*pw + clm + 0.03 + pw/2., yh-4.*pass, Form("Underflow  %.0f",dispersionIC2_EBmod[i]->GetBinContent(0)));
    dispersionIC_EBmod[i]->GetXaxis()->SetRangeUser(0.002+mapMin, mapMax-0.002); // to avoid overlapping labels
  }
  TLatex* tex2 = new TLatex();
  tex2->SetTextSize(0.05);
  tex2->SetTextFont(42);  
  tex2->DrawLatex(0.1,0.95, canvasLabel1.c_str());
  tex2->SetTextColor(kRed+2);
  tex2->DrawLatex(0.4,0.95, canvasLabel2.c_str());
  canvas->SaveAs(Form("%s/calibMap_EB_%s_dispersionIC_module.pdf",outDir.c_str(),canvasSuffix.c_str()));
  canvas->SaveAs(Form("%s/calibMap_EB_%s_dispersionIC_module.png",outDir.c_str(),canvasSuffix.c_str()));   

  // vector<TPad*> pads;
  // for (Int_t i = 0; i < 4; i++) {
  //   pads.push_back(new TPad());
  //   pads[i] = (TPad*) canvas->cd(i+1);    
  //   // pad->SetLeftMargin(0.01);
  //   // pad->SetRightMargin(0.02);
  //   pads[i]->SetTickx(1);
  //   pads[i]->SetTicky(1);
  //   dispersionIC_EBmod[i]->SetStats(0);
  //   dispersionIC2_EBmod[i]->SetStats(0);
  //   dispersionIC_EBmod[i]->SetLineColor(kBlack);
  //   dispersionIC_EBmod[i]->GetYaxis()->SetRangeUser(0, 1.2*maxy);
  //   dispersionIC_EBmod[i]->Draw("HIST");
  //   dispersionIC_EBmod[i]->GetXaxis()->SetTitle("IC value");
  //   dispersionIC_EBmod[i]->GetXaxis()->SetTitleSize(0.06);
  //   dispersionIC_EBmod[i]->GetXaxis()->SetLabelSize(0.035);
  //   dispersionIC_EBmod[i]->GetXaxis()->SetTitleOffset(0.7);
  //   if (i == 0) {
  //     dispersionIC_EBmod[i]->GetYaxis()->SetTitle("Events");
  //     dispersionIC_EBmod[i]->GetYaxis()->SetTitleSize(0.06);
  //     dispersionIC_EBmod[i]->GetYaxis()->SetLabelSize(0.04);
  //     dispersionIC_EBmod[i]->GetYaxis()->SetTitleOffset(1.2);
  //   }
  //   dispersionIC2_EBmod[i]->SetLineColor(kRed+2);
  //   dispersionIC2_EBmod[i]->Draw("HIST SAME");
  //   tex->SetTextSize(0.025);
  //   tex->SetTextFont(42);
  //   tex->DrawLatex(i*0.25 + 0.1,0.85, Form("Module %d",i+1));
  //   pads[i]->RedrawAxis("sameaxis");
  //   pads[i]->Update();
  // }
  // canvas->SaveAs(Form("%s/calibMap_EB_%s_dispersionIC_module.pdf",outDir.c_str(),canvasSuffix.c_str()));
  // canvas->SaveAs(Form("%s/calibMap_EB_%s_dispersionIC_module.png",outDir.c_str(),canvasSuffix.c_str()));   


  // // prepare spread of IC for each module
  // vector<TProfile*> profileEB_phi_mod;
  // vector<TProfile*> profileEB2_phi_mod;
  // for (Int_t i = 1; i <=4; i++) {
  //   profileEB_phi_mod.push_back( new TProfile(Form("profileEB_phi_mod%d",i),"i#phi profile of IC in EB",360, 0.5, 360.5) );
  //   makeICprofileIphiFromMap(profileEB_phi_mod.back(), mapEB, true, true, false, i);
  //   profileEB2_phi_mod.push_back( new TProfile(Form("profileEB2_phi_mod%d",i),"i#phi profile of IC in EB",360, 0.5, 360.5) );
  //   makeICprofileIphiFromMap(profileEB2_phi_mod.back(), mapEB2, true, true, false, i);
  // }

  // for (Int_t i = 0; i < 4; i++) {
    
  //   TLatex* tex = new TLatex();
  //   tex->SetTextSize(0.025);
  //   tex->SetTextFont(42);
  //   tex->DrawLatex(0.15,0.8, Form("Module %d",i+1));
  //   cRatio1D->SaveAs(Form("%s/calibMap_EB_ratio_%s_1D_module%d.pdf",outDir.c_str(),canvasSuffix.c_str(),i+1));
  //   cRatio1D->SaveAs(Form("%s/calibMap_EB_ratio_%s_1D_module%d.png",outDir.c_str(),canvasSuffix.c_str(),i+1));   
  //   delete tex;
  // }


  delete hratioDistr;
  delete tex;

}

//================================================

void realDrawMapRatioEE(const string& outDir = "",
			const string& inputFile = "",
			const string& inputFile2 = "",
			const string& canvasSuffix = "",
			const string& mapName1 = "",
			const string& mapName2 = "",
			const Double_t mapMin = 0.95,
			const Double_t mapMax = 1.05,
			const Int_t is_EEp1_EEm2 = 1
			)   
{

  TH1::SetDefaultSumw2();
  TH1::StatOverflows(kTRUE);

  gStyle->SetPalette(55, 0);  // 55:raibow palette ; 57: kBird (blue to yellow) ; 107 kVisibleSpectrum ; 77 kDarkRainBow          
  gStyle->SetNumberContours(101); // default is 20 

  TFile* f = TFile::Open(inputFile.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error opening file \"" << inputFile << "\".\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }

  TH2F *mapEE = NULL;
  mapEE = (TH2F*) f->Get(mapName1.c_str());
  if (!mapEE) {
    cout << "Error: could not get EE histogram. End of programme" << endl;
    exit(EXIT_FAILURE);
  }
  mapEE->SetDirectory(0);
  f->Close();


  TFile* f2 = TFile::Open(inputFile2.c_str(),"READ");
  if (!f2 || !f2->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error opening file \"" << inputFile2 << "\".\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }

  TH2F *mapEE2 = NULL;
  mapEE2 = (TH2F*) f2->Get(mapName2.c_str());
  if (!mapEE2) {
    cout << "Error: could not get EE histogram 2. End of programme" << endl;
    exit(EXIT_FAILURE);
  }
  mapEE2->SetDirectory(0);
  f2->Close();

  TH2F* mapEE_new = new TH2F("mapEE_new","", 100, 0.5, 100.5, 100, 0.5, 100.5);
  TH2F* mapEE2_new = new TH2F("mapEE2_new","", 100, 0.5, 100.5, 100, 0.5, 100.5);

  Double_t mapEE_binContent = 0.0;
  Int_t bin = 0;
  
  for (Int_t iphi = 1; iphi <= 100; iphi++) {

    for (Int_t ieta = 1; ieta <= 100; ieta++) {

      bin = mapEE->FindFixBin(iphi,ieta);
      mapEE_binContent = mapEE->GetBinContent( bin );
      mapEE_new->SetBinContent(bin, mapEE_binContent);
      mapEE_binContent = mapEE2->GetBinContent( bin );
      mapEE2_new->SetBinContent(bin, mapEE_binContent);

    }

  }

  TH1F* hratioDistr = new TH1F("hratioDistr","",51, 0.9,1.1);

  TH2F* hRatio = (TH2F*) mapEE_new->Clone("ratio");
  Double_t ratioVal = 0.0;
  divideEEmap(hRatio,mapEE_new,mapEE2_new,true,0);
  //hRatio->Divide(mapEE2_new);
  for (Int_t i = 1; i <= hRatio->GetNbinsX(); i++) {
    for (Int_t j = 1; j <= hRatio->GetNbinsY(); j++) {
      ratioVal = hRatio->GetBinContent(i,j);
      if (ratioVal > 0.00001) hratioDistr->Fill(ratioVal);
    }
  }

  //TProfile* profileEB_phi_final = new TProfile("profileEB_phi_final","i#phi profile of IC ratio in EB",360, 0.5, 360.5);
  //makeICprofileIphiFromMap(profileEB_phi_final, hRatio, true, true, false);
  //drawDistribution(profileEB_phi_final, "i#phi", "mean IC", Form("calibMap_EB_ratio_%s_iphiProfile",canvasSuffix.c_str()), outDir, 0.5, 360.5, 700, 500, true);


  //EB
  Int_t xsizeCanvas = 800;
  Int_t ysizeCanvas = 800;

  string detIdName = (is_EEp1_EEm2 == 1) ? "EEp" :"EEm";

  TCanvas *cEB = new TCanvas("cEB","",xsizeCanvas,ysizeCanvas);
  // cEB->SetLeftMargin(0.16);
  // cEB->SetRightMargin(0.20);
  cEB->SetRightMargin(0.14);
  cEB->cd();
  hRatio->Draw("COLZ");
  hRatio->GetXaxis()->SetTitle("iX");
  hRatio->GetXaxis()->SetTitleSize(0.06);
  hRatio->GetXaxis()->SetTitleOffset(0.7);
  hRatio->GetYaxis()->SetTitle("iY");
  hRatio->GetYaxis()->SetTitleSize(0.06);
  hRatio->GetYaxis()->SetTitleOffset(0.8);
  hRatio->GetZaxis()->SetRangeUser(mapMin, mapMax);
  hRatio->SetStats(0);
  gPad->Update();
  cEB->SaveAs(Form("%s/calibMap_%s_ratio_%s.pdf",outDir.c_str(),detIdName.c_str(),canvasSuffix.c_str()));
  cEB->SaveAs(Form("%s/calibMap_%s_ratio_%s.png",outDir.c_str(),detIdName.c_str(),canvasSuffix.c_str()));
  delete cEB;

  //EB
  TCanvas *cRatio1D = new TCanvas("cRatio1D","");
  // cRatio1D->SetLeftMargin(0.16);
  // cRatio1D->SetRightMargin(0.20);
  cRatio1D->cd();
  cRatio1D->SetTickx(1);
  cRatio1D->SetTicky(1);
  cRatio1D->cd();
  cRatio1D->SetGridx(1);
  cRatio1D->SetGridy(1);
  hratioDistr->Draw("HIST");
  hratioDistr->GetXaxis()->SetTitle("ratio");
  hratioDistr->GetXaxis()->SetTitleSize(0.06);
  hratioDistr->GetXaxis()->SetTitleOffset(0.7);
  hratioDistr->GetYaxis()->SetTitle("Events");
  hratioDistr->GetYaxis()->SetTitleSize(0.06);
  hratioDistr->GetYaxis()->SetTitleOffset(0.8);
  //hratioDistr->SetStats(0);
  gPad->Update();
  cRatio1D->SaveAs(Form("%s/calibMap_%s_ratio_%s_1D.pdf",outDir.c_str(),detIdName.c_str(),canvasSuffix.c_str()));
  cRatio1D->SaveAs(Form("%s/calibMap_%s_ratio_%s_1D.png",outDir.c_str(),detIdName.c_str(),canvasSuffix.c_str()));
  delete cRatio1D;

  delete hratioDistr;

}

//===============================================

void makeICratio(//const string& outDir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/AlCaP0_Run2017_F_CCiter0/iter_6/2DMaps/ratio/",
		 const string& outDir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot_Legacy/ratioIC/AlCaP0_AllRun2017_condor_iter1__Over__AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_ext1_fromIter6_iter6/",
		 const string& canvasSuffix = "ratioIC",
		 const string& inputFile1 = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot_Legacy/AlCaP0_AllRun2017_condor/iter_1/2DMaps/ICmaps/IC_work_test/calibMap_EEm_norm1etaRing.root",
		 const string& inputFile2 = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_ext1_fromIter6/iter_6/2DMaps/ICmaps/IC_work_test/calibMap_EEm_norm1etaRing.root",
		 const string& mapName1 = "mapEEm_IC_norm1etaRing",  // mapEEm_IC_norm1etaRing, mapEB_original_Over_norm1eachModuleFoldSMallEB_plusMinusSeparate_norm1etaRing
		 const string& mapName2 = "mapEEm_IC_norm1etaRing",  // mapEEm_IC_norm1etaRing
		 const Double_t mapMin = 0.9,
		 const Double_t mapMax = 1.1,
		 const Int_t is_EB0_EEp1_EEm2 = 2, 
		 const string& canvasLabel1 = "IC set 1",
		 const string& canvasLabel2 = "IC set 2"
		 ) 
{

  system(Form("mkdir -p %s",outDir.c_str()));
  system(Form("cp /afs/cern.ch/user/m/mciprian/public/index.php %s",outDir.c_str()));

  if (is_EB0_EEp1_EEm2 > 0) {
    realDrawMapRatioEE(outDir, inputFile1, inputFile2, canvasSuffix, mapName1, mapName2,mapMin,mapMax, is_EB0_EEp1_EEm2);
  } else 
    realDrawMapRatio(outDir, inputFile1, inputFile2, canvasSuffix, mapName1, mapName2,mapMin,mapMax, canvasLabel1, canvasLabel2);

}
