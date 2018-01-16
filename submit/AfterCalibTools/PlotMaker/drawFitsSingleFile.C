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
s#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

using namespace std;
using namespace RooFit;

int Xtal_ID[170][360]={0};
int Xtal_Ieta[61200]={0};
int Xtal_Iphi[61200]={0};

int Xtal_Ix[14648]={0};
int Xtal_Iy[14648]={0};
int Xtal_Iz[14648]={0};

// WARNING: GETTING IX, IY AND ZSIDE FROM DETID HASHEDINDEX MUST BE TESTED
// IT LOOKS LIKE SOME INDICES ARE OVERWRITTEN
// Actually it looks like some objects are written twice in the file. In the following, a check is made to avoid drawing the same object twice

void drawFitsSingleFile(const string& fitResFileOnEos = "", const string& BarrelOrEndcap = "Barrel", const string& outputDIR = "./", const Int_t nFitsToPlot = 10, 
			const Int_t fitIndexToPlot = -1) {

  // if fitIndexToPlot >= 0 we just look for plot with that index and plot that one
  // otherwise just plot nFitsToPlot plots from fitResFileOnEos file

  // fitResFileOnEos is the file on EOS you want to run on. It must start with --> root://eoscms//eos/cms/store/....
  // outputDIR is the output directory were plotted fits are stored. It is created if not existing
  // nFitsToPlot is the number of fits to plot. A loop on objects in file is made and the first nFitsToPlot objects are drawn

  // EB
  if (BarrelOrEndcap == "Barrel") {

    for(int i = 0; i < 61200; i++)
      {

	EBDetId ebseed(EBDetId::detIdFromDenseIndex);
        int ieta = ebseed.ieta();
        int iphi = ebseed.iphi();		
	Xtal_Ieta[i] = ieta;
	Xtal_Iphi[i] = iphi;
	//	cout<<i<<"   "<<ieta<<"  "<<iphi<<endl;
      }

  } else {

    // EE
    for(int i = 0; i < 14648; i++)
      {

	// TO BE TESTED
	EEDetId eeseed(EEDetId::detIdFromDenseIndex);
        int ix = eeseed.ix();
        int iy = eeseed.iy();		       
	int iz = eeseed.zside();		
	Xtal_Ix[i] = ix;
	Xtal_Iy[i] = iy;
	Xtal_Iz[i] = iz;
	//	cout<<i<<"   "<<ix<<"  "<<iy<<endl;
      }

  }

  system(("mkdir -p " + outputDIR).c_str());
  

  TFile* f = TFile::Open(fitResFileOnEos.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<fitResFileOnEos<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }


  TCanvas *c = NULL;

  TIter next(f->GetListOfKeys());
  TKey *key = NULL;
  Int_t iloop = 0;
  string previousObjectName = "";

  while ((key = (TKey*)next()) && iloop < nFitsToPlot) {

    // it seems there are duplicated objects in these files, so check that next object name is different from previous one
    // could just use entries for which iloop is even, but then if we fix the bug we should remember to modify this patch

    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("RooPlot")) continue;

    RooPlot * xframe = (RooPlot*) key->ReadObj();
    if (!xframe) {
      cout << "Warning: RooPlot object not found in file. Skipping and going on with next object" <<endl;
      continue;
    }

    string rooplotname(xframe->GetName());
    if (fitIndexToPlot >= 0 && (rooplotname.find(Form("%d",fitIndexToPlot)) == string::npos)) continue;
    string rooplotTitle = "";
    string canvasname = "";

    if (iloop == 0) {
      
      previousObjectName = rooplotname;

    } else {

      if (rooplotname == previousObjectName) continue;
      else previousObjectName = rooplotname;

    }

    // get crystal index number from RooPlot name in file (name looks like Fit_n_<Number>_attempt<N>_rp, where N = 0,1,2... indicates the fit attempt 
    // (if first fails another one with modified parameter is made))
    string rooplotnameTag = "Fit_n_";
    string fitIndexStr = ""; 
    fitIndexStr.assign(rooplotname, rooplotnameTag.size(), string::npos);
    Int_t fitIndex = std::stoi(fitIndexStr);

    string fitAttemptNumber = ""; 
    fitAttemptNumber.assign(rooplotname, rooplotname.find("attempt")+7, rooplotname.find("attempt")+8);
    Int_t fitAttemptNumber_int = atoi(fitAttemptNumber.c_str());
    //cout << "fit attempt: " << fitAttemptNumber_int << endl;

    if (BarrelOrEndcap == "Barrel") {

      stringstream os_ieta;
      stringstream os_iphi;
      
      os_ieta << Xtal_Ieta[fitIndex]; 
      os_iphi << Xtal_Iphi[fitIndex]; 
      
      string ss_iR = fitIndexStr;
      string ss_ieta = os_ieta.str();
      string ss_iphi = os_iphi.str();
      
      //rooplotTitle = "iR = " + ss_iR + " (iEta = " + ss_ieta + "  iPhi = " + ss_iphi + ")";
      rooplotTitle = "i#eta = " + ss_ieta + "  i#phi = " + ss_iphi;
      c = new TCanvas("c",rooplotname.c_str());
      canvasname = rooplotname + "_ieta" + ss_ieta + "_iphi" + ss_iphi + ".png";
    
    } else {

      // if (fitIndex < 6390 || fitIndex > 6399) continue;
      
      stringstream os_ix;
      stringstream os_iy;
      stringstream os_iz;
      
      os_ix << Xtal_Ix[fitIndex]; 
      os_iy << Xtal_Iy[fitIndex]; 
      os_iz << Xtal_Iz[fitIndex]; 
      
      string ss_iR = fitIndexStr;
      string ss_ix = os_ix.str();
      string ss_iy = os_iy.str();
      string ss_iz = os_iz.str();
      
      //rooplotTitle = "iR = " + ss_iR + " (iX = " + ss_ix + "  iY = " + ss_iy + "  iZ = " + ss_iz + ")";
      if (Xtal_Iz[fitIndex] > 0) rooplotTitle = "iX = " + ss_ix + "  iY = " + ss_iy + "  EE+";
      else rooplotTitle = "iX = " + ss_ix + "  iY = " + ss_iy + "  EE-";
      c = new TCanvas("c",rooplotname.c_str());
      canvasname = rooplotname + "_ix" + ss_ix + "_iy" + ss_iy + "_iz" + ss_iz + ".png";

    }
      
    if (xframe) {
      c->SetTickx(1);
      c->SetTicky(1);
      xframe->SetTitle(rooplotTitle.c_str());
      xframe->GetYaxis()->SetTitle("#gamma#gamma pairs / 0.004 GeV/c^{2}");
      xframe->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");
      xframe->Draw();
      c->SaveAs((outputDIR + canvasname).c_str());
    }


    // TFile* outputFile = NULL;

    // if ( (BarrelOrEndcap == "Barrel" && fitIndex == 30003) || (BarrelOrEndcap == "Endcap" && fitIndex == 6397) ) {
      
    //   outputFile = new TFile((outputDIR + "pi0Mass_singleXtal_Rooplot.root").c_str(),"UPDATE");
    //   if (!outputFile || outputFile->IsZombie()) {
    // 	cout << "Error: file not opened. Exit" << endl;
    // 	exit(EXIT_FAILURE);
    //   }
    //   outputFile->cd();
      
    //   xframe->Write();
    //   outputFile->Close();
    //   delete outputFile;
    
    // }

    delete c;

    // if for a given crystal we have more than one fit, do not increase loop counter when evaluating the attempts
    // iloop should refer to the nuber of crystals to plot, not actual fits to plot in total
    if (fitAttemptNumber_int == 0) iloop++;

    if (fitIndexToPlot >= 0 && (rooplotname.find(Form("%d",fitIndexToPlot)) != string::npos)) break;

  }


}
