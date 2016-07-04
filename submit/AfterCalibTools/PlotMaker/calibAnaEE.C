#define calibAnaEE_cxx
#include "calibAnaEE.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TPaletteAxis.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h                                                                                                                                       
#include <cstdio>
#include <cmath>
#include <sstream> //to use ostringstream to convert numbers to string in c++                                                                                          
#include <algorithm> // to use max

using namespace std;

Int_t getEtaRing(Int_t &ix, Int_t &iy, Int_t &zside) {

  Int_t etaRing = -1;

  string fileName = "eerings_modified.dat";
  ifstream inputFile(fileName.c_str());

  Int_t a,b,c,d;  // file format is --> a b c d                                                                                                 
  if (inputFile.is_open()) {

    while ((etaRing == -1) && (inputFile >> a >> b >> c >> d)) {

      if (c == zside && a == ix && b == iy) {  // first check zside, because it will remove half of the line to check
	etaRing = d;
      }

    }

  } else {
    std::cout << "Error: could not open file " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }


  inputFile.close();

  return etaRing;

}

void calibAnaEE::Loop(string whichEE, string dirName, string iterNumber)
{
//   In a ROOT session, you can do:
//      root> .L calibAnaEE.C
//      root> calibAnaEE t
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

  fChain->SetBranchStatus("*",0);  // disable all branches                                                                
  fChain->SetBranchStatus("*",1);

  TH2D *hSignal = new TH2D("hSignal","Signal in EE",102,0,102,102,0,102);
  TH2D *hBackground = new TH2D("hBackground","Background in EE",102,0,102,102,0,102);
  TH2D *SoverB = new TH2D("SoverB","S/B in EE",102,0,102,102,0,102);                     
  TH2D *SoverSqrtSplusB = new TH2D("SoverSqrtSplusB","S/sqrt(S+B) in EE",102,0,102,102,0,102);
  TH2D *SigmaMeanOverMean = new TH2D("SigmaMeanOverMean","sigma(fit_mean)/fit_mean in EE",102,0,102,102,0,102); //mean is the fit mean (should be the pi0 peak mass)  
  TH2D *mean = new TH2D("mean","fit_mean in EE",102,0,102,102,0,102);
  TH2D *sigma = new TH2D("sigma","fit_sigma in EE",102,0,102,102,0,102);

  // TProfile *hSignal_iRProfile = new TProfile("hSignal_iRProfile","Signal profile in EE",52,0,52);
  // TProfile *hBackground_iRProfile = new TProfile("hBackground_iRProfile","Background profile in EE",52,0,52);
  // TProfile *SoverB_iRProfile = new TProfile("SoverB_iRProfile","S/B profile in EE",52,0,52);
  // TProfile *SoverSqrtSplusB_iRProfile = new TProfile("SoverSqrtSplusB_iRProfile","S/sqrt(S+B) profile in EE",52,0,52);
  // TProfile *SigmaMeanOverMean_iRProfile = new TProfile("SigmaMeanOverMean_iRProfile","sigma(fit_mean)/fit_mean profile in EE",52,0,52);
  // TProfile *mean_iRProfile = new TProfile("mean_iRProfile","fit_mean profile in EE",52,0,52);
  // TProfile *sigma_iRProfile = new TProfile("sigma_iRProfile","fit_sigma profile in EE",52,0,52);

  TProfile *hSignal_etaRingProfile = new TProfile("hSignal_etaRingProfile","Signal profile in EE",40,0,40);
  TProfile *hBackground_etaRingProfile = new TProfile("hBackground_etaRingProfile","Background profile in EE",40,0,40);
  TProfile *SoverB_etaRingProfile = new TProfile("SoverB_etaRingProfile","S/B profile in EE",40,0,40);
  TProfile *SoverSqrtSplusB_etaRingProfile = new TProfile("SoverSqrtSplusB_etaRingProfile","S/sqrt(S+B) profile in EE",40,0,40);
  TProfile *SigmaMeanOverMean_etaRingProfile = new TProfile("SigmaMeanOverMean_etaRingProfile","sigma(fit_mean)/fit_mean profile in EE",40,0,40);
  TProfile *mean_etaRingProfile = new TProfile("mean_etaRingProfile","fit_mean profile in EE",40,0,40);
  TProfile *sigma_etaRingProfile = new TProfile("sigma_etaRingProfile","fit_sigma profile in EE",40,0,40);


  vector<TH2D*> th2dVector;
  th2dVector.push_back(hSignal);
  th2dVector.push_back(hBackground);
  th2dVector.push_back(SoverB);
  th2dVector.push_back(SoverSqrtSplusB);
  th2dVector.push_back(SigmaMeanOverMean);
  th2dVector.push_back(mean);
  th2dVector.push_back(sigma);

  vector<Double_t> th2dMinZaxisVector;
  vector<Double_t> th2dMaxZaxisVector;  // code will use default when this is lower than the maximum 

  // vector<TProfile*> profileIRVector;
  // profileIRVector.push_back(hSignal_iRProfile);
  // profileIRVector.push_back(hBackground_iRProfile);
  // profileIRVector.push_back(SoverB_iRProfile);
  // profileIRVector.push_back(SoverSqrtSplusB_iRProfile);
  // profileIRVector.push_back(SigmaMeanOverMean_iRProfile);
  // profileIRVector.push_back(mean_iRProfile);
  // profileIRVector.push_back(sigma_iRProfile);

  vector<TProfile*> profileEtaRingVector;
  profileEtaRingVector.push_back(hSignal_etaRingProfile);
  profileEtaRingVector.push_back(hBackground_etaRingProfile);
  profileEtaRingVector.push_back(SoverB_etaRingProfile);
  profileEtaRingVector.push_back(SoverSqrtSplusB_etaRingProfile);
  profileEtaRingVector.push_back(SigmaMeanOverMean_etaRingProfile);
  profileEtaRingVector.push_back(mean_etaRingProfile);
  profileEtaRingVector.push_back(sigma_etaRingProfile);

  vector<string> profileYaxisTitle;
  profileYaxisTitle.push_back("events");
  profileYaxisTitle.push_back("events");
  profileYaxisTitle.push_back("S/B");
  profileYaxisTitle.push_back("S/#sqrt{S+B}");
  profileYaxisTitle.push_back("#sigma(mean)/mean");
  profileYaxisTitle.push_back("mean [GeV]");
  profileYaxisTitle.push_back("#sigma [GeV]");

  th2dMinZaxisVector.push_back(0.0);
  th2dMinZaxisVector.push_back(0.0); 
  th2dMinZaxisVector.push_back(0.0);
  th2dMinZaxisVector.push_back(0.0);
  th2dMinZaxisVector.push_back(0.0);
  th2dMinZaxisVector.push_back(0.13);
  th2dMinZaxisVector.push_back(0.005);

  //  string whichEE = "EEp";  // can be EEp or EEm, change this to select which Endcap to use

  Double_t normalizedS = 0.0;
  Double_t normalizedB = 0.0;

  Int_t etaRing = -1;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
 
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;                                                                                                                                  
    if (jentry % 50000 == 0) cout << jentry << endl;

    if (whichEE == "EEp") {
      if (zside < 0) continue;
    } else {
      if (zside > 0) continue;
    }

    if ((abs(Backgr) > 0.00001) && (abs(Signal) > 0.00001)) { // avoid empty crystals due to masked towers or whatever               

      normalizedS = Signal * fit_Snorm;
      normalizedB = Backgr * fit_Bnorm;

      // Double_t ixshift = ix -51;
      // Double_t iyshift = iy -51;

      etaRing = getEtaRing(ix,iy,zside);  // will return -1 if eta ring is not found for any reason

      //      Double_t iR = sqrt(ixshift * ixshift + iyshift * iyshift);

      hSignal->Fill((Double_t)ix,(Double_t)iy,max(th2dMinZaxisVector[0],(Double_t)normalizedS));
      hBackground->Fill((Double_t)ix,(Double_t)iy,max(th2dMinZaxisVector[1],(Double_t)normalizedB));
      SoverB->Fill((Double_t)ix,(Double_t)iy,max(th2dMinZaxisVector[2],(Double_t)normalizedS/normalizedB));
      SoverSqrtSplusB->Fill((Double_t)ix,(Double_t)iy, max(th2dMinZaxisVector[3],(Double_t)normalizedS/sqrt(normalizedS + normalizedB)));
      SigmaMeanOverMean->Fill((Double_t)ix,(Double_t)iy, max(th2dMinZaxisVector[4],(Double_t)fit_mean_err/fit_mean));
      mean->Fill((Double_t)ix,(Double_t)iy, max(th2dMinZaxisVector[5],(Double_t)fit_mean));
      sigma->Fill((Double_t)ix,(Double_t)iy, max(th2dMinZaxisVector[6],(Double_t)fit_sigma));

      // hSignal_iRProfile->Fill((Double_t)iR,normalizedS);
      // hBackground_iRProfile->Fill((Double_t)iR,normalizedB);
      // SoverB_iRProfile->Fill(iR,normalizedS/normalizedB);
      // SoverSqrtSplusB_iRProfile->Fill(iR, normalizedS/sqrt(normalizedS + normalizedB));
      // SigmaMeanOverMean_iRProfile->Fill(iR, fit_mean_err/fit_mean);
      // mean_iRProfile->Fill(iR, fit_mean);
      // sigma_iRProfile->Fill(iR, fit_sigma);

      hSignal_etaRingProfile->Fill((Double_t)etaRing,normalizedS);
      hBackground_etaRingProfile->Fill((Double_t)etaRing,normalizedB);
      SoverB_etaRingProfile->Fill(etaRing,normalizedS/normalizedB);
      SoverSqrtSplusB_etaRingProfile->Fill(etaRing, normalizedS/sqrt(normalizedS + normalizedB));
      SigmaMeanOverMean_etaRingProfile->Fill(etaRing, fit_mean_err/fit_mean);
      mean_etaRingProfile->Fill(etaRing, fit_mean);
      sigma_etaRingProfile->Fill(etaRing, fit_sigma);
 
    }

  }

  th2dMaxZaxisVector.push_back(hSignal->GetBinContent(hSignal->GetMaximumBin()));
  th2dMaxZaxisVector.push_back(hBackground->GetBinContent(hBackground->GetMaximumBin()));
  th2dMaxZaxisVector.push_back(10e9);
  th2dMaxZaxisVector.push_back(10e9);
  th2dMaxZaxisVector.push_back(0.1);  // when this value is very large (bigger than the default) use the default to plot axis
  th2dMaxZaxisVector.push_back(0.16);
  th2dMaxZaxisVector.push_back(0.020); // 15 MeV

  for (UInt_t i = 0; i < th2dVector.size(); i++) {

    if (th2dMinZaxisVector[i] < 0.0) th2dMinZaxisVector[i] = 0.0;
    th2dVector[i]->SetMinimum(th2dMinZaxisVector[i]);
    //    profileIRVector[i]->SetMinimum(th2dMinZaxisVector[i]);
    profileEtaRingVector[i]->SetMinimum(th2dMinZaxisVector[i]);

    // if the maximum choosen by user is bigger than default, don't do anything, otherwise set the user value as the maximum   
    if (th2dMaxZaxisVector[i] < th2dVector[i]->GetBinContent(th2dVector[i]->GetMaximumBin())) {

      th2dVector[i]->SetMaximum(th2dMaxZaxisVector[i]);
      //profileIRVector[i]->SetMaximum(th2dMaxZaxisVector[i]);
      profileEtaRingVector[i]->SetMaximum(th2dMaxZaxisVector[i]);

    }

  }


  string wwwPath = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/" + dirName + "/" + iterNumber + "/2DMaps/Endcap/";
  wwwPath += whichEE + "/";

  vector<TCanvas *> c;
  vector<TCanvas *> cprof;

  for ( UInt_t i = 0; i < th2dVector.size(); i++ ) {

    c.push_back(NULL);
    string canvasName(th2dVector[i]->GetName());
    canvasName = "c_" + canvasName;
    c[i] = new TCanvas(canvasName.c_str(),canvasName.c_str());
    th2dVector[i]->Draw("COLZ");
    th2dVector[i]->GetXaxis()->SetTitle("iX");
    th2dVector[i]->GetXaxis()->SetTitleSize(0.06);
    th2dVector[i]->GetXaxis()->SetTitleOffset(0.7);
    th2dVector[i]->GetYaxis()->SetTitle("iY");
    th2dVector[i]->GetYaxis()->SetTitleSize(0.06);
    th2dVector[i]->GetYaxis()->SetTitleOffset(0.8);
    th2dVector[i]->SetStats(0);
    th2dVector[i]->Draw("COLZ");
    // after drawing, fix the palette                                                                                                                                  
    gPad->Update();
    TPaletteAxis *palette = (TPaletteAxis*)th2dVector[i]->GetListOfFunctions()->FindObject("palette");
    // the following lines move the palette. Choose the values you need for the position.                                                                              
    palette->SetX1NDC(0.91);
    palette->SetX2NDC(0.94);
    //palette->SetY1NDC(0.2);                                                                                                                                          
    //palette->SetY2NDC(0.8);                                                                                                                                          
    gPad->Modified();
    gPad->Update();
    // end of palette fixes      
    string name = wwwPath + th2dVector[i]->GetName() + "_" + whichEE;
    c[i]->SaveAs((name + ".pdf").c_str());
    c[i]->SaveAs((name + ".png").c_str());
    //cout << th2dVector[i]->GetName() << " :  mean , RMS -->  " << th2dVector[i]->GetMean(3) << " , " << th2dVector[i]->GetRMS(3) << endl;         

    cprof.push_back(NULL);
    canvasName = "c_" + string(profileEtaRingVector[i]->GetName());
    cprof[i] = new TCanvas(canvasName.c_str(),canvasName.c_str());
    profileEtaRingVector[i]->Draw("HE");
    profileEtaRingVector[i]->GetXaxis()->SetTitle("#eta Ring");
    profileEtaRingVector[i]->GetXaxis()->SetTitleSize(0.06);
    profileEtaRingVector[i]->GetXaxis()->SetTitleOffset(0.7);
    profileEtaRingVector[i]->GetYaxis()->SetTitle( profileYaxisTitle[i].c_str() );    
    profileEtaRingVector[i]->GetYaxis()->SetTitleSize(0.055);                                                                     
    profileEtaRingVector[i]->GetYaxis()->SetTitleOffset(0.8);                                                                                                          
    profileEtaRingVector[i]->SetStats(0);
    profileEtaRingVector[i]->Draw("HE");
    name = wwwPath + profileEtaRingVector[i]->GetName() + "_" + whichEE;
    cprof[i]->SaveAs((name + ".pdf").c_str());
    cprof[i]->SaveAs((name + ".png").c_str());

  }

}


Int_t main(int argc, char* argv[]) {

  TChain *chain = new TChain("calibEE");

  // string dirName = "";
  // string iter_number_as_string = "";

  string fileToChain = "";  
  // string eosPath = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/";
  // string dirName = "AlcaP0_fromRun273158_2016_v2";
  // string iter_num = "6";
  // string whichEE = "EEp";

  string whichEE(argv[1]);
  string path(argv[2]);
  string eosPath = "root://eoscms//eos/cms" + path;
  string dirName(argv[3]);
  string iter_num(argv[4]);
  string tagName(argv[5]);


  // if (argc > 1 ) {

  //   for (Int_t i = 1; i < argc; i++) {   // look at all possible options passed

  //     string thisArgument(argv[i]);

  //     if (thisArgument  == "EEp" ) {
  // 	whichEE = thisArgument;
  // 	cout << "Analysis of EE+" << endl; 
  //     }

  //     if (thisArgument  == "EEm" ) {
  // 	whichEE = thisArgument;
  // 	cout << "Analysis of EE-" << endl;
  //     }

      // if (thisArgument  == "-dir" ) {
      // 	dirName.assign(argv[i+1]);  // assign directory name after -dir option --> -dir <name_dir>
      // 	cout << "dirName = " << dirName << endl;
      // 	i += 1;
      // }

      // if (thisArgument  == "-iter" ) {
      //   iter_number_as_string.assign(argv[i+1]); 
      //   cout << "iter number = " << iter_number_as_string << endl;
      //   i += 1;
      // }


  //  }

  //}

  string Result;          // string which will contain the result                                                                                                    

  for (Int_t i = 0; i < 8; i++) {

    ostringstream convert;   // stream used for the conversion                                                                                                         
    convert << i;      // insert the textual representation of 'i' in the characters in the stream                                                                     
    Result = convert.str();
    fileToChain = eosPath + dirName + "/iter_" + iter_num + "/" + tagName + "Endcap_" + Result + "_calibMap.root";
    chain->Add(TString(fileToChain.c_str()));

  }

  //chain->Add(TString(fileToChain.c_str()));
  //cout << "chain = " << chain << endl;

  string iterNumber = "iter_" + iter_num;
  calibAnaEE *ana = new calibAnaEE(chain);
  ana->Loop(whichEE, dirName, iterNumber);

  return 0;

}
