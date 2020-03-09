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

//==========================================================


void getRooplotWithIndex(const string& fitResFileOnEos = "", 
			 const bool isEB = true, 
			 const string& outputDIR = "./", 
			 const Int_t rooplotIndex = 1, 
			 const string& outputFileName = ""
			 ) {

  TFile* f = TFile::Open(fitResFileOnEos.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<fitResFileOnEos<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  // stringstream ss_index;
  // ss_index << rooplotIndex;
  // string str_index = ss_index.str();

  RooPlot * xframe = (RooPlot*) f->Get(Form("Fit_n_%d_attempt0_rp",rooplotIndex));
  if (!xframe) {
    cout << "Warning: RooPlot object with name \"Fit_n_" << rooplotIndex << "_attempt0_rp\" not found in file " << fitResFileOnEos<< ". Exit" <<endl;
    exit(EXIT_FAILURE);
  } else {
    // check if there are other attempts
    for (Int_t irp = 3; irp >= 1; irp--) {
      string plotInFile = Form("Fit_n_%d_attempt%d_rp",rooplotIndex,irp);
      RooPlot * xframe_irp = (RooPlot*) f->Get(plotInFile.c_str());
      if (xframe_irp) {
	cout << "Found object with name " << plotInFile << ". Will plot this one" << endl;
	xframe = (RooPlot*) f->Get(plotInFile.c_str());
      }
    }
  }

  if (xframe) {
    xframe->GetYaxis()->SetTitle("#gamma#gamma pairs / 0.004 GeV/c^{2}");
    xframe->GetXaxis()->SetTitle("#gamma#gamma invariant mass [GeV/c^{2}]");
  }

  TFile* outputFile = new TFile((outputDIR + outputFileName).c_str(),"UPDATE");
  if (!outputFile || outputFile->IsZombie()) {
    cout << "Error: file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }
  outputFile->cd();
  
  xframe->Write();
  outputFile->Close();
  delete outputFile;

  f->Close();
  delete f;

}

//===============================================

void drawRooPlotFromFile(const string& inputDir = "", 
			 const bool isEB = true, 
			 const string &inputFileName = "", 
			 const Int_t rooplotIndex = 1, 
			 const bool isPi0 = true,
			 const double lumi= 0.18,
			 const double eta = 1.6, 
			 const Int_t year = 2017
			 ) 
{

  //double lumi = 0.18; //in fb-1

  TFile* f = TFile::Open((inputDir+inputFileName).c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<< inputDir+inputFileName <<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  string plotInFile(Form("Fit_n_%d_attempt0_rp",rooplotIndex));

  RooPlot * xframe = (RooPlot*) f->Get(plotInFile.c_str());
  if (!xframe) {
    cout << "Warning: RooPlot object with name \"Fit_n_" << rooplotIndex << "_attempt0_rp\" not found in file";
    cout << inputFileName;
    cout << "Abort" << endl;
    exit(EXIT_FAILURE);
  } else {
    // check if there are other attempts
    for (Int_t irp = 3; irp >= 1; irp--) {
      plotInFile = Form("Fit_n_%d_attempt%d_rp",rooplotIndex,irp);
      RooPlot * xframe_irp = (RooPlot*) f->Get(plotInFile.c_str());
      if (xframe_irp) {
	cout << "Found object with name " << plotInFile << ". Will plot this one" << endl;
	xframe = (RooPlot*) f->Get(plotInFile.c_str());
      }
    }      
  }

  string canvasname = isEB ? "pi0MassEBxtal" : "pi0MassEExtal";
  if (not isPi0) canvasname = isEB ? "etaMassEBxtal" : "etaMassEExtal"; 
  canvasname += Form("_index_%d",rooplotIndex);
  TCanvas *canvas = new TCanvas("canvas","",700,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  // canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);
  canvas->SetLeftMargin(0.18);

  RooHist* data = xframe->getHist("data");
  Double_t maxY = 1.2 * data->getYAxisMax();
  xframe->GetYaxis()->SetRangeUser(0,maxY);
  xframe->GetYaxis()->SetTitle("Number of #gamma#gamma pairs");
  xframe->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV)");

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

  xframe->GetXaxis()->SetLabelSize(0.04);
  xframe->GetXaxis()->SetTitleSize(0.05);
  xframe->GetYaxis()->SetTitleOffset(1.45);
  xframe->GetYaxis()->SetTitleSize(0.055);
  xframe->Draw();

  TLegend *leg = NULL;
  if (isPi0) {
    if (isEB) leg = new TLegend(0.6,0.7,0.95,0.9);
    else leg = new TLegend(0.60,0.3,0.95,0.5);
    //else leg = new TLegend(0.60,0.35,0.95,0.55);
    //else leg = new TLegend(0.60,0.7,0.95,0.9);
  } else {
    if (isEB) leg = new TLegend(0.63,0.66,0.93,0.91); // new TLegend(0.2,0.25,0.5,0.5); 
    else      leg = new TLegend(0.50,0.25,0.95,0.5);
  }
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h1,"data","PLE");
  //leg->AddEntry(h2,"S + B","LF");
  leg->AddEntry(h2,"Fit model","LF");
  leg->AddEntry(h4,"Signal","LF");
  leg->AddEntry(h3,"Background","LF");
  leg->Draw("same");

  // // for "ECAL Barrel"
  // // try TLegend
  TLegend *leg2 = NULL;
  // if (isEB) leg2 = new TLegend(0.12,0.75,0.48,0.90);
  // else leg2 = new TLegend(0.12,0.75,0.5,0.90);
  // // if (isEB) leg2 = new TLegend(0.50,0.1,0.99,0.3);
  // // //else  leg2 = new TLegend(0.31,0.11,0.8,0.31);
  // // else  leg2 = new TLegend(0.4,0.15,0.9,0.35);
  // leg2->SetFillColor(0);
  // leg2->SetFillStyle(0);
  // leg2->SetBorderSize(0);
  // if (isEB) {
  //   leg2->AddEntry((TObject*) 0, Form("ECAL Barrel"),""); // remove Crystal
  //   leg2->AddEntry((TObject*) 0, Form("#eta = %.1g",eta),""); // remove Crystal
  // } else {
  //   leg2->AddEntry((TObject*) 0, Form("ECAL Endcap"),""); // remove Crystal
  //   leg2->AddEntry((TObject*) 0, Form("#eta = %.2g ",eta),""); // remove Crystal
  // }
  // leg2->Draw("same");

  // or TLatex
  TLatex lat;
  std::string line = "";
  lat.SetNDC();
  lat.SetTextSize(0.045);
  lat.SetTextFont(42);
  lat.SetTextColor(1);
  float xmin(0.22), yhi(0.85), ypass(0.05);
  string detector = isEB ? "Barrel" : "Endcap";
  line = Form("ECAL %s",detector.c_str());
  lat.DrawLatex(xmin,yhi, line.c_str());
  line = Form("#eta = %.2f",eta);
  lat.DrawLatex(xmin,yhi-ypass, line.c_str());
  lat.DrawLatex(xmin,yhi-2.0*ypass, Form("%s#rightarrow#gamma#gamma",isPi0 ? "#pi^{0}" : "#eta^{0}"));

  canvas->RedrawAxis("sameaxis");

  if (lumi < 0.0) {
    CMS_lumi(canvas,"",true,false,0,0,0,year);
  } else {
    if (lumi < 1.0) CMS_lumi(canvas,Form("%.2f",lumi),true,false,0,0,0,year);
    else CMS_lumi(canvas,Form("%.1f",lumi),true,false,0,0,0,year);
  }
  setTDRStyle();
  
  canvas->SaveAs((inputDir + canvasname + ".pdf").c_str());
  canvas->SaveAs((inputDir + canvasname + ".png").c_str());
  canvas->SaveAs((inputDir + canvasname + ".root").c_str());
  canvas->SaveAs((inputDir + canvasname + ".C").c_str());

  TFile* outputFile = new TFile((inputDir + canvasname + ".root").c_str(),"RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    cout << "Error: file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }
  outputFile->cd();  
  canvas->Write();
  outputFile->Close();
  delete outputFile;


  delete canvas;
  delete h1; 
  delete h2; 
  delete h3; 
  delete leg;
  if (leg2 != 0) delete leg2;

  f->Close();
  delete f;

}

//===============================================

void printSignificanceInFile(const string& calibMapFile = "", 
			     const bool isEB = true, 
			     const int xtalIndex = 0, 
			     const string& outputDirEB = "", 
			     const string& significanceFileName = ""
			     ) {

  string treeName = isEB ? "calibEB" : "calibEE";

  // TFile* inputFile = new TFile(calibMapFile.c_str(),"READ");
  // if (!inputFile || inputFile->IsZombie()) {
  //   cout << "Error: file not opened. Exit" << endl;
  //   exit(EXIT_FAILURE);
  // }

  //TTreeReader reader(treeName.c_str(), inputFile);

  TChain* chain = new TChain(treeName.c_str());
  chain->Add(calibMapFile.c_str());

  TTreeReader reader(chain);

  string ieta_or_ix_str = isEB ? "ieta_" : "ix_";
  string iphi_or_iy_str = isEB ? "iphi_" : "iy_";

  TTreeReaderValue<Int_t> hashedIndex (reader,"hashedIndex_");
  TTreeReaderValue<Int_t> ieta_or_ix (reader,ieta_or_ix_str.c_str());
  TTreeReaderValue<Int_t> iphi_or_iy (reader,iphi_or_iy_str.c_str());
  TTreeReaderValue<Int_t>* iz = NULL; if (not isEB) iz = new TTreeReaderValue<Int_t>(reader,"zside_");
  TTreeReaderValue<Float_t> Signal (reader,"Signal_");
  TTreeReaderValue<Float_t> Backgr (reader,"Backgr_");
  TTreeReaderValue<Float_t> fit_Snorm (reader,"fit_Snorm_");
  TTreeReaderValue<Float_t> fit_Bnorm (reader,"fit_Bnorm_");
  TTreeReaderValue<Float_t> fit_mean (reader,"fit_mean_");
  TTreeReaderValue<Float_t> fit_mean_err (reader,"fit_mean_err_");
  TTreeReaderValue<Float_t> fit_sigma (reader,"fit_sigma_");

  bool entryFound = false;

  while(reader.Next() && not entryFound) {

    // look for the specific crystal
    if (*hashedIndex != xtalIndex) continue;
    else {

      /////////////////////////////////
      // print in file some information
      ofstream significanceFile(significanceFileName.c_str(),ios::out);

      if ( !significanceFile.is_open() ) {

	cout<<"Error: unable to open file " << significanceFileName <<" !"<<endl;
	exit(EXIT_FAILURE);
     
      } else {      

	double sig = *Signal * *fit_Snorm;
	double bkg = *Backgr * *fit_Bnorm;
	double SoverB = sig/bkg;
	double SoverSqrtSplusB = sig / sqrt(sig+bkg);
	if (isEB) {
	  significanceFile << "crystal in EB" << endl;
	  significanceFile << "hashedIndex\t" << *hashedIndex << endl;
	  significanceFile << "ieta\t" << *ieta_or_ix << endl;
	  significanceFile << "iphi\t" << *iphi_or_iy << endl;
	} else {
	  significanceFile << "crystal in EE" << endl;
	  significanceFile << "hashedIndex\t" << *hashedIndex << endl;
	  significanceFile << "iX\t" << *ieta_or_ix << endl;
	  significanceFile << "iY\t" << *iphi_or_iy << endl;
	  significanceFile << "iZ\t" << **iz << endl;
	}
	significanceFile << setprecision(3) << "S\t" << sig << endl;
	significanceFile << setprecision(3) << "B\t" << bkg << endl;
	significanceFile << setprecision(3) << "S/B\t" << SoverB << endl;
	significanceFile << setprecision(3) << "S/sqrt(S+B)\t" << SoverSqrtSplusB << endl;
	significanceFile << setprecision(3) << "mean(fit)\t" << *fit_mean << endl;
	significanceFile << setprecision(3) << "mean_err(fit)\t" << *fit_mean_err << endl;
	significanceFile << setprecision(3) << "mean_err/mean(fit)\t" << 100. * (*fit_mean_err/(*fit_mean)) << "%" << endl;
	significanceFile << setprecision(3) << "sigma(fit)\t" << *fit_sigma << endl;
	
	significanceFile.close();

      }
      /////////////////////////////////

      entryFound = true;

    }

  }

  // inputFile->Close();
  // delete inputFile;

}


//===============================================


void manageRooPlotFromFile(const string& dirName = "AlCaEta_2016_ULrereco", 
			   const string& outDirName = "plot_approve_UL2016data_Eta", 
			   const bool usePi0 = false, 
			   //const string& dirName = "AlCaP0_AllRun2017_condor_fixEBm16", 
			   //const string& outDirName = "plot_approve_full2017data_Pi0_legacyReRecoCalib", 
			   //const bool usePi0 = true, 
			   // const string& dirName = "AlCaP0_Run2018D_goldenJson_13_09_2018", 
			   // const string& outDirName = "plot_approve_2018D_pi0", 
			   // const bool usePi0 = true, 
			   const Int_t skip_EB1_EE2 = 0, 
			   const double lumi = -10.0, 
			   const int whichIteration = 0, 
			   const string& subdirTag = "",
			   const Int_t year = 2016,
			   const string& eosPath = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero_Run2/mciprian/" 
			   ) {

  // intLumi is in /fb, use 2 digits after .

  // gROOT->ProcessLine(".L saveRooPlotFromFile.C++");
  // gROOT->ProcessLine(".L makeRooPlotFromFile.C++");

  bool isPi0 = usePi0;
  // safety check in case user makes mistakes with dirName and isPi0
  if (dirName.find("AlCaP0") != string::npos) isPi0 = true;
  if (dirName.find("AlCaEta") != string::npos) isPi0 = false;

  int EBxtalIndex = 30003; //32429;
  string EBfitFileIndex = "15";//"16"; // need to find a way to derive it from EBxtalIndex
  double etaEB = -0.03; // 1.0;// would be negative but ok
  //int EExtalIndex = 12001; //12001;
  //string EEfitFileIndex = "6"; //"6"; // need to find a way to derive it from EExtalIndex
  //double etaEE = 2.5;
  int EExtalIndex = 8155; //8000;     //14018; //8155; //12001;
  string EEfitFileIndex = "4";  //4"; // "7"; // 4//"6"; // need to find a way to derive it from EExtalIndex
  double etaEE = 1.83;// 1.63;// 1.83;

  // if (not isPi0) {
  //   EBxtalIndex = 30107;
  //   EExtalIndex = 6397;
  // }

  string iter = string(Form("%d",whichIteration));

  // string outputDirEB = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/" + dirName + "/iter_" + iter + "/fitResPlots/Barrel/";
  // string outputDirEE = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/" + dirName + "/iter_" + iter + "/fitResPlots/Endcap/";
  string outputDirEB = "/afs/cern.ch/user/m/mciprian/www/pi0calib/" + outDirName + "/" + dirName + "/" + subdirTag + "/iter_" + iter + "/fitResPlots/Barrel/";
  string outputDirEE = "/afs/cern.ch/user/m/mciprian/www/pi0calib/" + outDirName + "/" + dirName + "/" + subdirTag + "/iter_" + iter + "/fitResPlots/Endcap/";

  if (skip_EB1_EE2 != 1) createPlotDirAndCopyPhp(outputDirEB);
  if (skip_EB1_EE2 != 2) createPlotDirAndCopyPhp(outputDirEE);

  // the following file is actually the output file from getRooplotWithIndex()
  // it is then used as input by drawRooPlotFromFile()
  string fileTagName = "pi0Mass_singleXtal";
  if (not isPi0) fileTagName = "etaMass_singleXtal";
  string inputFileNameEB = fileTagName + Form("_index_%d_Rooplot.root",EBxtalIndex);
  string inputFileNameEE = fileTagName + Form("_index_%d_Rooplot.root",EExtalIndex);

  string inputFitEosDIR_EB = eosPath + dirName + "/iter_" + iter + "/" + dirName + "_Barrel_" + EBfitFileIndex + "_fitRes.root";
  string inputFitEosDIR_EE = eosPath + dirName + "/iter_" + iter + "/" + dirName + "_Endcap_" + EEfitFileIndex + "_fitRes.root";

  // read file from eos, take fit for a given crystal and save the RooPlot
  if (skip_EB1_EE2 != 1) getRooplotWithIndex(inputFitEosDIR_EB, true, outputDirEB, EBxtalIndex, inputFileNameEB);
  if (skip_EB1_EE2 != 2) getRooplotWithIndex(inputFitEosDIR_EE, false, outputDirEE, EExtalIndex, inputFileNameEE);

  // it seems that the first time CMS_lumi is used the settings are screwed up
  // produce a dummy plot (either do not save it or remove it)
  //double lumi = 0.18; //in fb-1
  TCanvas*ctmp = new TCanvas("ctmp","");
  ctmp->cd();
  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  htmp1->Fill(0.5);
  htmp1->Draw("H");
  CMS_lumi(ctmp,Form("%.2f",lumi),false,false);
  setTDRStyle();
  delete htmp1;
  delete ctmp;

  // here we go with the real part
  if (skip_EB1_EE2 != 1) drawRooPlotFromFile(outputDirEB, true, inputFileNameEB, EBxtalIndex, isPi0, lumi, etaEB, year);
  if (skip_EB1_EE2 != 2) drawRooPlotFromFile(outputDirEE, false, inputFileNameEE, EExtalIndex, isPi0, lumi, etaEE, year);

  string calibMapFile = eosPath + dirName + "/iter_" + iter + "/" + dirName + "_calibMap.root";
  string significanceFileNameEB = fileTagName + Form("_index_%d_significance.txt",EBxtalIndex);
  string significanceFileNameEE = fileTagName + Form("_index_%d_significance.txt",EExtalIndex);

  if (skip_EB1_EE2 != 1) printSignificanceInFile(calibMapFile, true, EBxtalIndex, outputDirEB, outputDirEB + significanceFileNameEB);
  if (skip_EB1_EE2 != 2) printSignificanceInFile(calibMapFile, false, EExtalIndex, outputDirEB, outputDirEE + significanceFileNameEE);

}
