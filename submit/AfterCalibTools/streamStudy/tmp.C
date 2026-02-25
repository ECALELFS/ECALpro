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
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"

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

  TCanvas * xframe = (TCanvas*) f->Get(Form("Fit_n_%d_attempt0_c",rooplotIndex));
  if (!xframe) {
    cout << "Warning: RooPlot object with name \"Fit_n_" << rooplotIndex << "_attempt0_c\" not found in file. Exit" <<endl;
    exit(EXIT_FAILURE);
  }

  // if (xframe) {
  //   xframe->GetYaxis()->SetTitle("#gamma#gamma pairs / 0.004 GeV/c^{2}");
  //   xframe->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");
  // }

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

void drawRooPlotFromFile(const string& inputDir = "", const bool isEB = true, const string &inputFileName = "", const Int_t rooplotIndex = 1, const bool isPi0 = true,
			 const double lumi= 0.18) {

  //double lumi = 0.18; //in fb-1

  TFile* f = TFile::Open((inputDir+inputFileName).c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<< inputDir+inputFileName <<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  string plotInFile(Form("Fit_n_%d_attempt0_c",rooplotIndex));

  TCanvas* canvas = (TCanvas*) f->Get(plotInFile.c_str());
  if (!canvas) {
    cout << "Warning: Canvas object with name \"" << plotInFile <<"\" not found in file";
    cout << inputFileName;
    cout << "Abort" << endl;
    exit(EXIT_FAILURE);
  }

  string canvasname = isEB ? "pi0MassEBxtal" : "pi0MassEExtal";
  if (not isPi0) canvasname = isEB ? "etaMassEBxtal" : "etaMassEExtal"; 
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetRightMargin(0.06);

  RooPlot * xframe = NULL;
  TKey *key = NULL;
  // TIter next(canvas->GetListOfPrimitives());
  // canvas->GetListOfPrimitives(->Print();
  // while ( (key = (TKey*)next()) ){
  //   TClass *cl = gROOT->GetClass(key->GetClassName());
  //   cout << "name " << cl->GetName() << endl;
  //   if (cl->InheritsFrom("RooPlot")) {
  //     xframe = (RooPlot*) key->ReadObj();
  //     if (!xframe) {
  // 	cout << "Warning: RooPlot object not found in file. Skipping and going on with next object" <<endl;
  // 	exit(EXIT_FAILURE);
  //     } else {
  // 	cout << "name " << xframe->GetName() << endl;
  //     }
  //   }
  // }

  cout << "Abort" << endl;
  exit(EXIT_FAILURE);


  //  RooPlot* xframe = (RooPlot*) canvas->GetPrimitive(Form("Fit_n_%d_attempt0_rp",rooplotIndex));
  xframe->GetYaxis()->SetTitle("#gamma#gamma pairs / 0.004 GeV/c^{2}");
  xframe->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");

  // to add a dummy legend
  TH1D* h1 = new TH1D("h1","",1,0,1);
  TH1D* h2 = new TH1D("h2","",1,0,1);
  TH1D* h3 = new TH1D("h3","",1,0,1);
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

  xframe->GetXaxis()->SetLabelSize(0.04);
  xframe->GetXaxis()->SetTitleSize(0.05);
  xframe->GetYaxis()->SetTitleOffset(1.1);
  xframe->GetYaxis()->SetTitleSize(0.05);
  xframe->Draw();

  TLegend *leg = NULL;
  if (isPi0) {
    if (isEB) leg = new TLegend(0.55,0.7,0.95,0.9);
    //    else leg = new TLegend(0.50,0.25,0.95,0.5);
    else leg = new TLegend(0.50,0.65,0.95,0.9);
  } else {
    leg = new TLegend(0.50,0.25,0.95,0.5);
  }
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h1,"data","PLE");
  leg->AddEntry(h2,"Signal + background","LF");
  leg->AddEntry(h3,"background only","LF");
  leg->Draw("same");
  TLegend *leg2 = NULL;
  if (isEB) leg2 = new TLegend(0.50,0.1,0.99,0.3);
  else  leg2 = new TLegend(0.50,0.1,0.99,0.3);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  if (isEB) leg2->SetHeader("ECAL Barrel Crystal");
  else leg2->SetHeader("ECAL Endcap Crystal");
  leg2->Draw("same");
  canvas->RedrawAxis("sameaxis");

  if (lumi < 1.0) CMS_lumi(canvas,Form("%.2f",lumi),true,false);
  else CMS_lumi(canvas,Form("%.1f",lumi),true,false);
  setTDRStyle();
  
  canvas->SaveAs((inputDir + canvasname + ".pdf").c_str());
  canvas->SaveAs((inputDir + canvasname + ".png").c_str());

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
  delete leg2;

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


void tmp(const string& dirName = "AlCaP0_Run2017A_runs296966to296980_v2", const bool usePi0 = true, const Int_t skip_EB1_EE2 = 0, const double lumi = 3.94, const int whichIteration = 0) {

  // intLumi is in /fb, use 2 digits after .

  // gROOT->ProcessLine(".L saveRooPlotFromFile.C++");
  // gROOT->ProcessLine(".L makeRooPlotFromFile.C++");

  string eosPath = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/";

  bool isPi0 = usePi0;
  // safety check in case user makes mistakes with dirName and isPi0
  if (dirName.find("AlCaP0") != string::npos) isPi0 = true;
  if (dirName.find("AlCaEta") != string::npos) isPi0 = false;

  int EBxtalIndex = 30003;
  string EBfitFileIndex = "15"; // need to find a way to derive it from EBxtalIndex
  int EExtalIndex = 6397;
  string EEfitFileIndex = "3"; // need to find a way to derive it from EExtalIndex

  if (not isPi0) {
    EBxtalIndex = 30107;
    EExtalIndex = 6397;
  }

  string iter = string(Form("%d",whichIteration));

  // string outputDirEB = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/" + dirName + "/iter_" + iter + "/fitResPlots/Barrel/";
  // string outputDirEE = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/" + dirName + "/iter_" + iter + "/fitResPlots/Endcap/";
  string outputDirEB = "/afs/cern.ch/user/m/mciprian/www/pi0calib/plot_EPS_2017/" + dirName + "/iter_" + iter + "/fitResPlots/Barrel/";
  string outputDirEE = "/afs/cern.ch/user/m/mciprian/www/pi0calib/plot_EPS_2017/" + dirName + "/iter_" + iter + "/fitResPlots/Endcap/";

  if (skip_EB1_EE2 != 1) createPlotDirAndCopyPhp(outputDirEB);
  if (skip_EB1_EE2 != 2) createPlotDirAndCopyPhp(outputDirEE);

  // the following file is actually the output file from getRooplotWithIndex()
  // it is then used as input by drawRooPlotFromFile()
  string fileTagName = "pi0Mass_singleXtal";
  if (not isPi0) fileTagName = "etaMass_singleXtal";
  string inputFileName = fileTagName + "_Rooplot.root";

  string inputFitEosDIR_EB = eosPath + dirName + "/iter_" + iter + "/" + dirName + "_Barrel_" + EBfitFileIndex + "_fitRes.root";
  string inputFitEosDIR_EE = eosPath + dirName + "/iter_" + iter + "/" + dirName + "_Endcap_" + EEfitFileIndex + "_fitRes.root";


  // read file from eos, take fit for a given crystal and save the RooPlot
  if (skip_EB1_EE2 != 1) getRooplotWithIndex(inputFitEosDIR_EB, true, outputDirEB, EBxtalIndex, inputFileName);
  if (skip_EB1_EE2 != 2) getRooplotWithIndex(inputFitEosDIR_EE, false, outputDirEE, EExtalIndex, inputFileName);

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
  if (skip_EB1_EE2 != 1) drawRooPlotFromFile(outputDirEB, true, inputFileName, EBxtalIndex, isPi0, lumi);
  if (skip_EB1_EE2 != 2) drawRooPlotFromFile(outputDirEE, false, inputFileName, EExtalIndex, isPi0, lumi);

  string calibMapFile = eosPath + dirName + "/iter_" + iter + "/" + dirName + "_calibMap.root";
  string significanceFileName = fileTagName + "_significance.txt";

  if (skip_EB1_EE2 != 1) printSignificanceInFile(calibMapFile, true, EBxtalIndex, outputDirEB, outputDirEB + significanceFileName);
  if (skip_EB1_EE2 != 2) printSignificanceInFile(calibMapFile, false, EExtalIndex, outputDirEB, outputDirEE + significanceFileName);

}
