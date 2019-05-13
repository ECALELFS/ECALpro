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
#include <TH2D.h>
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
#include <utility>

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#define N_FIT_IN_FILE 2000    // says how many fits per job were made, the number of the fit file containing the fit index is fitFile = fit_index/N_FIT_IN_FILE

using namespace std;


void doPlotTCanvas(const string& filePath = "", 
		   const Int_t& iterNum = 0, 
		   const string& dirName = "", 
		   const int ieta = 1, 
		   const int iphi = 1,
		   const Bool_t isEB = true,
		   const Bool_t isMC_EoverEtrue = true,
		   const string& outDir = "",
		   const Int_t nFitPerFile = 50,
		   const bool foldSM = true
		   ) {


  // it is assumed the fit attempt is 0 (so it is the first fit for that crystal)

  // EBDetId ebseed(EBDetId::detIdFromDenseIndex(fitIndex));
  // int ieta = ebseed.ieta();
  // int iphi = ebseed.iphi();		

  Int_t fitIndex = -1;
  EBDetId ebseed(ieta,iphi,0);
  if (foldSM) {
    fitIndex = ebseed.ic()-1;
    //printf("EBDetId.ieta() = %d    EBDetId.iphi() = %d    EBDetId.iphiSM() = %d\n",ebseed.ieta(),ebseed.iphi(),ebseed.iphiSM());
  } else {
    fitIndex = ebseed.hashedIndex();
  }
  //Int_t fitFileIndex = (Int_t) fitIndex/N_FIT_IN_FILE;
  Int_t fitFileIndex = (Int_t) fitIndex/nFitPerFile;

  string filename= filePath + dirName + Form("/iter_%d/",iterNum) + dirName + Form("_%s_%d_fitRes.root",(isEB ? "Barrel" : "Endcap"),fitFileIndex);
  //string filename= filePath + Form("_%s_%d_fitRes.root",(isEB ? "Barrel" : "Endcap"),fitFileIndex);

  TFile* f = TFile::Open(filename.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<filename<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  string canvasNameInFile = "";
  string canvasNameInFile_g2 = "";
  if (isMC_EoverEtrue) {
    canvasNameInFile = string(Form("Fit_n_%d_attempt0_g1_c",fitIndex));
    canvasNameInFile_g2 = string(Form("Fit_n_%d_attempt0_g2_c",fitIndex));
  } else {
    canvasNameInFile = string(Form("Fit_n_%d_attempt0_c",fitIndex));
  }
  

  TCanvas* canvas = (TCanvas*) f->Get(canvasNameInFile.c_str());
  if (!canvas) {
    cout << "Could not get canvas " << canvasNameInFile << " (ieta,iphi = " << ieta << "," << iphi << ") " << endl;
    cout << "Maybe the crystal is a dead one, for which there is no fit. Skipping this crystal" << endl;
    return; //exit(EXIT_FAILURE);
  }
  canvas->SaveAs(Form("%s/%s_ieta%d_iphi%d.png",outDir.c_str(),canvas->GetName(),ieta,iphi));
  canvas->SaveAs(Form("%s/%s_ieta%d_iphi%d.pdf",outDir.c_str(),canvas->GetName(),ieta,iphi));

  TCanvas* canvas_g2 = nullptr;
  if (isMC_EoverEtrue) {
    canvas_g2 = (TCanvas*) f->Get(canvasNameInFile_g2.c_str());
    if (!canvas) {
      cout << "Could not get canvas " << canvasNameInFile_g2 << " (ieta,iphi = " << ieta << "," << iphi << ") " << endl;
      cout << "Maybe the crystal is a dead one, for which there is no fit. Skipping this crystal" << endl;
      return; //exit(EXIT_FAILURE);
    }
    canvas_g2->SaveAs(Form("%s/%s_ieta%d_iphi%d.png",outDir.c_str(),canvas_g2->GetName(),ieta,iphi));
    canvas_g2->SaveAs(Form("%s/%s_ieta%d_iphi%d.pdf",outDir.c_str(),canvas_g2->GetName(),ieta,iphi));
  }
  

  f->Close();
  delete f;

}


void plotTCanvas(const string& dirName = "pi0CC_2017_EoverEtrue_foldSM_nFit10_onlyEB_testNewFitsMay2019", 
		 const Int_t iterNum = 0, 
		 const Bool_t isMC_EoverEtrue = true,
		 const string& outDir_base = "/afs/cern.ch/user/m/mciprian/www/pi0calib/",
		 const string& filePath = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero_Run2/mciprian/",
		 const string& EoEtrueFolderName = "CC_EoverEtrue_2017",
		 const Int_t nFitPerFile = 10,
		 const bool foldSM = true
		 ) {

  // is foldSM is true, it means there where only 1700 fits (1 SM)

  //const string& filePath = "/afs/cern.ch/work/m/mciprian/myEcalElf/2017_ECALpro/calib2017/CMSSW_9_4_1/src/CalibCode/submit/tmp_rootFiles_EoverEtrue_foldSM/pi0Gun_MC_EoverEtrue_foldSM";

  //const string& filePath = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/";
  const Bool_t isEB = true;
  string outDir = outDir_base + EoEtrueFolderName + "/" + dirName + "/fits/";
  if (not isMC_EoverEtrue) outDir = outDir_base + "ICplot/" + dirName + Form("/iter_%d/",iterNum) + "fitResPlots/" + Form("%s/", isEB ? "Barrel" : "Endcap");

  system(Form("mkdir -p %s",outDir.c_str()));
  system(Form("cp /afs/cern.ch/user/m/mciprian/public/index.php %s",outDir.c_str()));

  vector< std::pair<Int_t, Int_t> > xtal_ieta_iphi;

  if (isMC_EoverEtrue) {

    xtal_ieta_iphi.push_back( std::make_pair( 1, 1) );
    xtal_ieta_iphi.push_back( std::make_pair(65, 1) );
    xtal_ieta_iphi.push_back( std::make_pair(85, 1) );
    xtal_ieta_iphi.push_back( std::make_pair(15, 1) );
    xtal_ieta_iphi.push_back( std::make_pair(35, 1) );
    xtal_ieta_iphi.push_back( std::make_pair(55, 1) );
    xtal_ieta_iphi.push_back( std::make_pair(75, 1) );
    xtal_ieta_iphi.push_back( std::make_pair(85,10) );
    xtal_ieta_iphi.push_back( std::make_pair(84, 4) );
    xtal_ieta_iphi.push_back( std::make_pair(15,10) );
    xtal_ieta_iphi.push_back( std::make_pair(35,10) );
    xtal_ieta_iphi.push_back( std::make_pair(55,10) );
    xtal_ieta_iphi.push_back( std::make_pair(75,10) );

    xtal_ieta_iphi.push_back( std::make_pair( 8, 7) );
    xtal_ieta_iphi.push_back( std::make_pair( 8, 8) );
    xtal_ieta_iphi.push_back( std::make_pair( 8, 9) );
    xtal_ieta_iphi.push_back( std::make_pair( 7, 7) );
    xtal_ieta_iphi.push_back( std::make_pair( 7, 8) );
    xtal_ieta_iphi.push_back( std::make_pair( 7, 9) );
    xtal_ieta_iphi.push_back( std::make_pair( 9, 7) );
    xtal_ieta_iphi.push_back( std::make_pair( 9, 8) );
    xtal_ieta_iphi.push_back( std::make_pair( 9, 9) );

    xtal_ieta_iphi.push_back( std::make_pair(83, 18) );
    xtal_ieta_iphi.push_back( std::make_pair(81, 3) );
    xtal_ieta_iphi.push_back( std::make_pair(46, 5) );

  } else {

    xtal_ieta_iphi.push_back( std::make_pair(55, 187) );
    xtal_ieta_iphi.push_back( std::make_pair(55, 188) );
    xtal_ieta_iphi.push_back( std::make_pair(55, 189) );
    xtal_ieta_iphi.push_back( std::make_pair(54, 187) );
    xtal_ieta_iphi.push_back( std::make_pair(54, 188) );
    xtal_ieta_iphi.push_back( std::make_pair(54, 189) );

    xtal_ieta_iphi.push_back( std::make_pair(55, 195) );
    xtal_ieta_iphi.push_back( std::make_pair(54, 195) );
    xtal_ieta_iphi.push_back( std::make_pair(53, 195) );
    xtal_ieta_iphi.push_back( std::make_pair(55, 196) );
    xtal_ieta_iphi.push_back( std::make_pair(54, 196) );
    xtal_ieta_iphi.push_back( std::make_pair(53, 196) );

  }

  for (UInt_t i = 0; i < xtal_ieta_iphi.size(); ++i) {
    doPlotTCanvas(filePath, iterNum, dirName, xtal_ieta_iphi[i].first, xtal_ieta_iphi[i].second, isEB, isMC_EoverEtrue, outDir, nFitPerFile, foldSM);
    //doPlotTCanvas(filePath, xtal_ieta_iphi[i].first, xtal_ieta_iphi[i].second, isEB, isMC_EoverEtrue, outDir);
  }


}
