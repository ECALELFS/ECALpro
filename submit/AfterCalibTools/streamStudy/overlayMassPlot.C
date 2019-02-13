#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
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

#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "./CMS_lumi.h"

using namespace std;
using namespace RooFit;

static const string PhpToCopy = "/afs/cern.ch/user/m/mciprian/www/index.php";

//==========================================================

void createPlotDirAndCopyPhp(const string& outputDIR) {

  if (outputDIR != "./") {
    system(("mkdir -p " + outputDIR).c_str());
    system(("cp "+ PhpToCopy + " " + outputDIR).c_str());
  }

}


//====================================================

void overlayMassPlot(const string& outdir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/massComparison/"
		     ) {


  createPlotDirAndCopyPhp(outdir);  
  /////////////////////////////////////
  // it seems that the first time CMS_lumi is used the settings are screwed up
  // produce a dummy plot (either do not save it or remove it)
  //double lumi = 0.18; //in fb-1
  TCanvas*ctmp = new TCanvas("ctmp","");
  ctmp->cd();
  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  htmp1->Fill(0.5);
  htmp1->Draw("H");
  CMS_lumi(ctmp,"35.9",false,false);
  setTDRStyle();
  delete htmp1;
  delete ctmp;

  ///////////////////////////
  // some configs

  Double_t lumi = 12.9;

  Bool_t isPi0 = true;
  Bool_t isEB = false;
  
  string subdet = isEB ? "Barrel" : "Endcap";
  string detector = Form("ECAL %s", subdet.c_str());
  Double_t eta = isEB ? -0.82 : 1.63;  
  Double_t sigmaOverM_multifit = isEB ? 9.57 : 11.9;
  Double_t sigmaOverM_weight   = isEB ? 10.6 : 14.8;

  string id1 = isEB ? "14007" : "14018";
  string id2 = isEB ? "14007" : "14018";

  string input1 = Form("/afs/cern.ch/user/m/mciprian/www/pi0calib/testForEmanuele/AlCaP0_Run2018D_goldenJson_13_09_2018/multifit/iter_0/fitResPlots/%s/pi0Mass_singleXtal_index_%s_Rooplot.root",subdet.c_str(), id1.c_str());
  string input2 = Form("/afs/cern.ch/user/m/mciprian/www/pi0calib/testForEmanuele/AlCaP0_Run2018D_goldenJson_13_09_2018_weight/weight/iter_0/fitResPlots/%s/pi0Mass_singleXtal_index_%s_Rooplot.root",subdet.c_str(), id2.c_str());

  string name1 = "Fit_n_" + id1 + "_attempt0_rp";
  string name2 = "Fit_n_" + id2 + "_attempt0_rp";

  string canvasname = Form("massComparison__%s_%s_%s", isEB ? "EB" : "EE", id1.c_str(), id2.c_str()); 

  string legentry1 = "multifit";
  string legentry2 = "weight";

  //////////////////////////////

  // file 1

  TFile* f1 = TFile::Open(input1.c_str(),"READ");
  if (!f1 || !f1->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<input1<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
  RooPlot * xframe1 = (RooPlot*) f1->Get(name1.c_str());
  if (!xframe1) {
    cout << "Warning: RooPlot object with name " << name1 << " not found in file " << input1 << ". Exit" <<endl;
    exit(EXIT_FAILURE);
  }
  //f1->Close();

  // file 2 
  TFile* f2 = TFile::Open(input2.c_str(),"READ");
  if (!f2 || !f2->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<input2<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
  RooPlot * xframe2 = (RooPlot*) f2->Get(name2.c_str());
  if (!xframe2) {
    cout << "Warning: RooPlot object with name " << name2 << " not found in file " << input2 << ". Exit" <<endl;
    exit(EXIT_FAILURE);
  }
  //f2->Close();

  TCanvas *canvas = new TCanvas("canvas","",700,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  // canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);
  canvas->SetLeftMargin(0.18);

  RooHist* data1 = xframe1->getHist("data"); // can take fit model with (RooCurve*) xframe->getCurve("model")
  RooCurve* model1 = xframe1->getCurve("model");
  RooHist* data2 = xframe2->getHist("data"); // can take fit model with (RooCurve*) xframe->getCurve("model")
  RooCurve* model2 = xframe2->getCurve("model");

  Double_t maxY = 1.3 * std::max(data1->getYAxisMax(),data2->getYAxisMax());
  xframe1->GetYaxis()->SetRangeUser(0,maxY);
  xframe1->GetYaxis()->SetTitle("Number of #gamma#gamma pairs");
  xframe1->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV)");

  // to add a dummy legend
  TH1D* h1 = new TH1D("h1","",1,0,1);
  TH1D* h2 = new TH1D("h2","",1,0,1);
  TH1D* h3 = new TH1D("h3","",1,0,1);
  TH1D* h4 = new TH1D("h4","",1,0,1);

  // data
  h1->SetStats(0);
  h1->SetLineColor(kBlack);
  h1->SetMarkerColor(kBlack);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  // S+B
  h2->SetStats(0);
  h2->SetLineColor(kBlue);
  h2->SetLineWidth(2);
  // B only
  h3->SetStats(0);
  h3->SetLineColor(kRed);
  h3->SetLineWidth(2);
  h3->SetLineStyle(2);
  // S only
  h4->SetStats(0);
  h4->SetLineColor(kGreen+2);
  h4->SetLineWidth(2);
  h4->SetLineStyle(2);

  xframe1->GetXaxis()->SetLabelSize(0.04);
  xframe1->GetXaxis()->SetTitleSize(0.05);
  xframe1->GetYaxis()->SetTitleOffset(1.45);
  xframe1->GetYaxis()->SetTitleSize(0.055);
  xframe1->Draw();
  data2->SetMarkerStyle(22);
  data2->SetMarkerColor(kGray+2);
  data2->SetMarkerSize(1);
  data2->Draw("p SAME");
  model2->SetLineColor(kOrange+1);
  model2->Draw("L SAME");

  // Double_t x = 0;
  // Double_t y = 0; 
  // Double_t yerr = 0;
  // data2->GetPoint(0,x,y);
  // yerr = data2->GetErrorY(0);
  // cout << x << ": y = " << y << " +/- " << yerr << endl;

  TLegend *leg = NULL;
  if (isPi0) {
    if (isEB) leg = new TLegend(0.6,0.66,0.95,0.9);
    else leg = new TLegend(0.60,0.21,0.95,0.45);
    //else leg = new TLegend(0.60,0.68,0.95,0.9);
  } else {
    leg = new TLegend(0.50,0.25,0.95,0.5);
  }
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h1,"data (multifit)","PLE");
  //leg->AddEntry(h2,"S + B","LF");
  leg->AddEntry(h2,"Fit model","LF");
  leg->AddEntry(h4,"Signal","LF");
  leg->AddEntry(h3,"Background","LF");
  leg->AddEntry(data2,"data (weight)","PLE");
  leg->AddEntry(model2,"Model (weight)","L");
  leg->Draw("same");


  // or TLatex
  TLatex lat;
  std::string line = "";
  lat.SetNDC();
  lat.SetTextSize(0.045);
  lat.SetTextFont(42);
  lat.SetTextColor(1);
  float xmin(0.22), yhi(0.85), ypass(0.05);
  line = Form("%s",detector.c_str());
  lat.DrawLatex(xmin,yhi, line.c_str());
  //line = Form("#eta = %.2g",eta);
  //lat.DrawLatex(xmin,yhi-ypass, line.c_str());

  TLatex lat2;
  lat2.SetNDC();
  lat2.SetTextSize(0.038);
  lat2.SetTextFont(42);
  lat2.SetTextColor(1);
  float xmin2 = isEB ? 0.22 : 0.62;
  float yhi2  = isEB ? 0.8 : 0.85;
  //line = "method: #sigma / #mu";
  //  line = Form("multifit: #sigma/M = %.3g%%",sigmaOverM_multifit);
  // lat2.DrawLatex(xmin2,yhi2, line.c_str());
  line = Form("#sigma/M_{multifit} = %.1f%%",sigmaOverM_multifit);
  lat2.DrawLatex(xmin2,yhi2, line.c_str());
  line = Form("#sigma/M_{weight} = %.1f%%",sigmaOverM_weight);
  lat2.DrawLatex(xmin2,yhi2-ypass, line.c_str());


  canvas->RedrawAxis("sameaxis");

  if (lumi < 1.0) CMS_lumi(canvas,Form("%.2f",lumi),true,true);
  else CMS_lumi(canvas,Form("%.1f",lumi),true,true);
  setTDRStyle();
  
  canvas->SaveAs((outdir + canvasname + ".pdf").c_str());
  canvas->SaveAs((outdir + canvasname + ".png").c_str());
  canvas->SaveAs((outdir + canvasname + ".C").c_str());
  canvas->SaveAs((outdir + canvasname + ".root").c_str());

  delete canvas;
  delete h1; 
  delete h2; 
  delete h3; 
  delete leg;
  

}
