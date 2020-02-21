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

void realComparePeak(const string& outdir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/massComparison_Run2/",
		     const Bool_t isEB = true,		    
		     const Bool_t isPi0 = true
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

  Double_t lumi = -1.0; // if < 0, nothing is printed on top right
  
  string subdet = isEB ? "Barrel" : "Endcap";
  string detector = Form("ECAL %s", subdet.c_str());
  Double_t eta = isEB ? 0.10 : 1.83;  
  Double_t sigmaOverM_2016 = 0.0;
  Double_t sigmaOverM_2017 = 0.0;
  Double_t sigmaOverM_2018 = 0.0;

  // EB: 30003, eta = -0.03
  // EB: 32429, eta = 0.10
  string id1 = isEB ? "32429" : "8155";
  string id2 = isEB ? "32429" : "8155";
  string id3 = isEB ? "32429" : "8155";

  string input1 = "";
  string input2 = "";
  string input3 = "";

  if (isPi0) {

    input1 = Form("/afs/cern.ch/user/m/mciprian/www/pi0calib/plot_approve_UL2016data/AlCaP0_2016_ULrereco_from0/iter_0/fitResPlots/%s/pi0Mass_singleXtal_index_%s_Rooplot.root",subdet.c_str(), id1.c_str());
    input2 = Form("/afs/cern.ch/user/m/mciprian/www/pi0calib/plot_approve_full2017data_Pi0_legacyReRecoCalib/AlCaP0_AllRun2017_condor_fixEBm16/iter_0/fitResPlots/%s/pi0Mass_singleXtal_index_%s_Rooplot.root",subdet.c_str(), id2.c_str());
    input3 = Form("/afs/cern.ch/user/m/mciprian/www/pi0calib/plot_approve_UL2018data/AlCaP0_2018_ULrereco_1every2/iter_0/fitResPlots/%s/pi0Mass_singleXtal_index_%s_Rooplot.root",subdet.c_str(), id2.c_str());

    // 30003
    sigmaOverM_2016 = isEB ? 8.74 : 8.64;
    sigmaOverM_2017 = isEB ? 10.5 : 9.25;
    sigmaOverM_2018 = isEB ? 10.97 : 8.91;
    // 32429
    sigmaOverM_2016 = isEB ? 7.72 : 8.64;
    sigmaOverM_2017 = isEB ? 9.20 : 9.25;
    sigmaOverM_2018 = isEB ? 9.85 : 8.91;

  } else {

    input1 = Form("/afs/cern.ch/user/m/mciprian/www/pi0calib/plot_approve_UL2016data_Eta/AlCaEta_2016_ULrereco/iter_0/fitResPlots/%s/etaMass_singleXtal_index_%s_Rooplot.root",subdet.c_str(), id1.c_str());
    input2 = Form("/afs/cern.ch/user/m/mciprian/www/pi0calib/plot_approve_UL2017data_Eta/AlCaEta_2017_ULrereco_all2017data/iter_0/fitResPlots/%s/etaMass_singleXtal_index_%s_Rooplot.root",subdet.c_str(), id2.c_str());
    input3 = Form("/afs/cern.ch/user/m/mciprian/www/pi0calib/plot_approve_UL2018data_Eta/AlCaEta_2018_ULrereco_all2018data/iter_0/fitResPlots/%s/etaMass_singleXtal_index_%s_Rooplot.root",subdet.c_str(), id2.c_str());

    // 30003
    sigmaOverM_2016 = isEB ? 4.15 : 6.34;
    sigmaOverM_2017 = isEB ? 4.44 : 6.16;
    sigmaOverM_2018 = isEB ? 4.70 : 5.79;
    //32429
    sigmaOverM_2016 = isEB ? 3.78 : 6.34;
    sigmaOverM_2017 = isEB ? 4.07 : 6.16;
    sigmaOverM_2018 = isEB ? 4.55 : 5.79;

  }


  string name1 = "Fit_n_" + id1 + "_attempt0_rp";
  string name2 = "Fit_n_" + id2 + "_attempt0_rp";
  string name3 = "Fit_n_" + id3 + "_attempt0_rp";

  string canvasname = Form("massComparison__%s_%s_%s_%s", isEB ? "EB" : "EE", id1.c_str(), id2.c_str(), id3.c_str()); 

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

  // file 3 
  TFile* f3 = TFile::Open(input3.c_str(),"READ");
  if (!f3 || !f3->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<input3<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }
  RooPlot * xframe3 = (RooPlot*) f3->Get(name3.c_str());
  if (!xframe3) {
    cout << "Warning: RooPlot object with name " << name3 << " not found in file " << input3 << ". Exit" <<endl;
    exit(EXIT_FAILURE);
  }
  //f3->Close();


  TCanvas *canvas = new TCanvas("canvas","",700,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  // canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);
  canvas->SetLeftMargin(0.18);

  RooHist* rhdata1 = xframe1->getHist("data"); // can take fit model with (RooCurve*) xframe->getCurve("model")
  RooCurve* model1 = xframe1->getCurve("model");
  RooHist* rhdata2 = xframe2->getHist("data"); // can take fit model with (RooCurve*) xframe->getCurve("model")
  RooCurve* model2 = xframe2->getCurve("model");
  RooHist* rhdata3 = xframe3->getHist("data"); // can take fit model with (RooCurve*) xframe->getCurve("model")
  RooCurve* model3 = xframe3->getCurve("model");

  Int_t npoints = rhdata1->GetN();
  Double_t xxmin = 0.0;
  Double_t xxmax = 0.0;
  Double_t yatxxmin = 0.0;
  Double_t yatxxmax = 0.0;
  rhdata2->GetPoint(0,xxmin,yatxxmin);
  rhdata2->GetPoint(npoints-1,xxmax,yatxxmax);
  xxmin -= rhdata2->GetErrorX(0);
  xxmax += rhdata2->GetErrorX(npoints-1);
  cout << "xlow " << xxmin << "   xhigh " << xxmax << endl;
  cout << "y " << yatxxmin << "   yerr " << rhdata2->GetErrorY(0) << endl;
  TH1D* data1 = new TH1D("data1","",npoints,xxmin,xxmax);
  TH1D* data2 = new TH1D("data2","",npoints,xxmin,xxmax);
  TH1D* data3 = new TH1D("data3","",npoints,xxmin,xxmax);
  for (Int_t i = 0; i < npoints; ++i) {
    Double_t xx = 0.0;
    Double_t yy = 0.0;
    rhdata1->GetPoint(i,xx,yy);
    data1->SetBinContent(i+1, yy);
    data1->SetBinError(  i+1, rhdata1->GetErrorY(i));
    rhdata2->GetPoint(i,xx,yy);
    data2->SetBinContent(i+1, yy);
    data2->SetBinError(  i+1, rhdata2->GetErrorY(i));
    rhdata3->GetPoint(i,xx,yy);
    data3->SetBinContent(i+1, yy);
    data3->SetBinError(  i+1, rhdata3->GetErrorY(i));
  }

  data1->Scale(1./data1->Integral());
  data2->Scale(1./data2->Integral());
  data3->Scale(1./data3->Integral());

  // for RooHist, which is a graph
  //Double_t maxY = std::max(data1->getYAxisMax(),data2->getYAxisMax());
  //maxY = 1.3 * std::max(data3->getYAxisMax(),maxY);
  // for TH1
  Double_t maxY = std::max(data1->GetMaximum(),data2->GetMaximum());
  maxY = 1.3 * std::max(data3->GetMaximum(),maxY);

  data1->GetXaxis()->SetRangeUser(xxmin,xxmax);
  data1->GetYaxis()->SetRangeUser(0,maxY);
  //data1->GetYaxis()->SetTitle("Number of #gamma#gamma pairs");
  data1->GetYaxis()->SetTitle("Arbitrary units");
  data1->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV)");

  data1->GetXaxis()->SetLabelSize(0.04);
  data1->GetXaxis()->SetTitleSize(0.05);
  data1->GetYaxis()->SetTitleOffset(1.45);
  data1->GetYaxis()->SetTitleSize(0.055);
  data1->SetMarkerStyle(20);
  data1->SetMarkerColor(kBlack);
  data1->SetLineColor(kBlack);
  data1->SetMarkerSize(1);
  data1->Draw("PE0");
  
  data2->SetMarkerStyle(22);
  data2->SetMarkerColor(kRed+2);
  data2->SetLineColor(kRed+2);
  data2->SetMarkerSize(1);
  data2->Draw("PE0 SAME");

  data3->SetMarkerStyle(23);
  data3->SetMarkerColor(kGreen+2);
  data3->SetLineColor(kGreen+2);
  data3->SetMarkerSize(1);
  data3->Draw("PE0 SAME");
  

  // Double_t x = 0;
  // Double_t y = 0; 
  // Double_t yerr = 0;
  // data2->GetPoint(0,x,y);
  // yerr = data2->GetErrorY(0);
  // cout << x << ": y = " << y << " +/- " << yerr << endl;

  TLegend *leg = NULL;
  if (isPi0) {
    if (isEB) leg = new TLegend(0.68,0.41,0.93,0.65);
    else leg = new TLegend(0.60,0.21,0.95,0.45);
    //else leg = new TLegend(0.60,0.68,0.95,0.9);
  } else {
    //if (isEB) leg = new TLegend(0.50,0.16,0.95,0.41);
    //else leg = new TLegend(0.50,0.25,0.95,0.5);
    if (isEB) leg = new TLegend(0.20,0.16,0.95,0.36);
    else leg = new TLegend(0.20,0.3,0.95,0.5);
  }
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  if (isPi0) {
    leg->AddEntry(data1,"2016","PLE");
    leg->AddEntry(data2,"2017","PLE");
    leg->AddEntry(data3,"2018","PLE");
  } else {
    leg->AddEntry(data1,Form("2016:  #sigma/M = %.1f%%",sigmaOverM_2016),"PLE");
    leg->AddEntry(data2,Form("2017:  #sigma/M = %.1f%%",sigmaOverM_2017),"PLE");
    leg->AddEntry(data3,Form("2018:  #sigma/M = %.1f%%",sigmaOverM_2018),"PLE");    
  }
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
  line = Form("#eta = %.3g",eta);
  lat.DrawLatex(xmin,yhi-ypass, line.c_str());

  TLatex lat2;
  lat2.SetNDC();
  lat2.SetTextSize(0.038);
  lat2.SetTextFont(42);
  lat2.SetTextColor(1);
  //line = "method: #sigma / #mu";
  //  line = Form("multifit: #sigma/M = %.3g%%",sigmaOverM_multifit);
  // lat2.DrawLatex(xmin2,yhi2, line.c_str());
  if (isPi0) {
    float xmin2 = isEB ? 0.6 : 0.62;
    float yhi2  = isEB ? 0.85 : 0.85;
    line = Form("#sigma/M_{2016} = %.1f%%",sigmaOverM_2016);
    lat2.DrawLatex(xmin2,yhi2, line.c_str());
    line = Form("#sigma/M_{2017} = %.1f%%",sigmaOverM_2017);
    lat2.DrawLatex(xmin2,yhi2-ypass, line.c_str());
    line = Form("#sigma/M_{2018} = %.1f%%",sigmaOverM_2018);
    lat2.DrawLatex(xmin2,yhi2-2.*ypass, line.c_str());
  }

  canvas->RedrawAxis("sameaxis");

  if (lumi < 0.0) CMS_lumi(canvas,"",true,true);
  else if (lumi < 1.0) CMS_lumi(canvas,Form("%.2f",lumi),true,true);
  else CMS_lumi(canvas,Form("%.1f",lumi),true,true);
  setTDRStyle();
  
  canvas->SaveAs((outdir + canvasname + ".pdf").c_str());
  canvas->SaveAs((outdir + canvasname + ".png").c_str());
  canvas->SaveAs((outdir + canvasname + ".C").c_str());
  canvas->SaveAs((outdir + canvasname + ".root").c_str());

  delete canvas;
  delete leg;
  

}


void comparePeak(const string& outdir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/massComparison_Run2/",
		 const int do_all0_pi0Only1_etaOnly2 = 1) 
{

  // pi0
  if (do_all0_pi0Only1_etaOnly2 != 2) {
    string outdirPi0 = outdir + "pi0/";
    realComparePeak(outdirPi0, true);
    realComparePeak(outdirPi0, false);
  }
  // eta
  if (do_all0_pi0Only1_etaOnly2 != 1) {
    string outdirEta = outdir + "eta/";
    realComparePeak(outdirEta, true,  false);
    realComparePeak(outdirEta, false, false);
  }

}
