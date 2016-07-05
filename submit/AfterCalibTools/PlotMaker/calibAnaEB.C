#define calibAnaEB_cxx
#include "calibAnaEB.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TPaletteAxis.h>

#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <math.h>
#include <sstream> //to use ostringstream to convert numbers to string in c++
#include <algorithm>  // to use max()

using namespace std;

void calibAnaEB::Loop(string dirName, string iterNumber)
{
//   In a ROOT session, you can do:
//      root> .L calibAnaEB.C
//      root> calibAnaEB t
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

  TH2D *hSignal = new TH2D("hSignal","Signal in EB",360,0,360,171,-85,86); 
  TH2D *hBackground = new TH2D("hBackground","Background in EB",360,0,360,171,-85,86); 
  TH2D *SoverB = new TH2D("SoverB","S/B in EB",360,0,360,171,-85,86); // up to 86 because crystal with ieta = 0 does not exists, it starts from 1
  TH2D *SoverSqrtSplusB = new TH2D("SoverSqrtSplusB","S/sqrt(S+B) in EB",360,0,360,171,-85,86);
  TH2D *SigmaMeanOverMean = new TH2D("SigmaMeanOverMean","sigma(fit_mean)/fit_mean in EB",360,0,360,171,-85,86); //mean is the fit mean (should be the pi0 peak mass)
  TH2D *mean = new TH2D("mean","fit_mean in EB",360,0,360,171,-85,86);
  TH2D *sigma = new TH2D("sigma","fit_sigma in EB",360,0,360,171,-85,86);

  // TProfile *hSignal_ietaProfile = new TProfile("hSignal_ietaProfile","S/B profile in EB",171,-85,86);
  // TProfile *hBackground_ietaProfile = new TProfile("hBackground_ietaProfile","S/B profile in EB",171,-85,86);
  // TProfile *SoverB_ietaProfile = new TProfile("SoverB_ietaProfile","S/B profile in EB",171,-85,86);
  // TProfile *SoverSqrtSplusB_ietaProfile = new TProfile("SoverSqrtSplusB_ietaProfile","S/sqrt(S+B) profile in EB",171,-85,86);
  // TProfile *SigmaMeanOverMean_ietaProfile = new TProfile("SigmaMeanOverMean_ietaProfile","sigma(fit_mean)/fit_mean profile in EB",171,-85,86);
  // TProfile *mean_ietaProfile = new TProfile("mean_ietaProfile","fit_mean profile in EB",171,-85,86);
  // TProfile *sigma_ietaProfile = new TProfile("sigma_ietaProfile","fit_sigma profile in EB",171,-85,86);

  TProfile *hSignal_etaProfile = new TProfile("hSignal_etaProfile","Signal profile in EB",170,-1.48,1.48);
  TProfile *hBackground_etaProfile = new TProfile("hBackground_etaProfile","Background profile in EB",170,-1.48,1.48);
  TProfile *SoverB_etaProfile = new TProfile("SoverB_etaProfile","S/B profile in EB",170,-1.48,1.48);
  TProfile *SoverSqrtSplusB_etaProfile = new TProfile("SoverSqrtSplusB_etaProfile","S/sqrt(S+B) profile in EB",170,-1.48,1.48);
  TProfile *SigmaMeanOverMean_etaProfile = new TProfile("SigmaMeanOverMean_etaProfile","sigma(fit_mean)/fit_mean profile in EB",170,-1.48,1.48);
  TProfile *mean_etaProfile = new TProfile("mean_etaProfile","fit_mean profile in EB",170,-1.48,1.48);
  TProfile *sigma_etaProfile = new TProfile("sigma_etaProfile","fit_sigma profile in EB",170,-1.48,1.48);

  TProfile *mean_iphiProfile = new TProfile("mean_iphiProfile","fit_mean profile in EB",360,0,360);

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

  // vector<TProfile*> profileIetaVector;
  // profileIetaVector.push_back(hSignal_ietaProfile);
  // profileIetaVector.push_back(hBackground_ietaProfile);
  // profileIetaVector.push_back(SoverB_ietaProfile);
  // profileIetaVector.push_back(SoverSqrtSplusB_ietaProfile);
  // profileIetaVector.push_back(SigmaMeanOverMean_ietaProfile);
  // profileIetaVector.push_back(mean_ietaProfile);
  // profileIetaVector.push_back(sigma_ietaProfile);

  vector<TProfile*> profileEtaVector;
  profileEtaVector.push_back(hSignal_etaProfile);
  profileEtaVector.push_back(hBackground_etaProfile);
  profileEtaVector.push_back(SoverB_etaProfile);
  profileEtaVector.push_back(SoverSqrtSplusB_etaProfile);
  profileEtaVector.push_back(SigmaMeanOverMean_etaProfile);
  profileEtaVector.push_back(mean_etaProfile);
  profileEtaVector.push_back(sigma_etaProfile);

  vector<string> profileYaxisTitle;
  profileYaxisTitle.push_back("events");
  profileYaxisTitle.push_back("events");
  profileYaxisTitle.push_back("S/B");
  profileYaxisTitle.push_back("S/#sqrt{S+B}");
  profileYaxisTitle.push_back("#sigma(mean)/mean");
  profileYaxisTitle.push_back("mean [GeV]");
  profileYaxisTitle.push_back("#sigma [GeV]");

  //  vector<TProfile*> profileEtaVector;

  // I was using these values starting for SoverB, SoverSqrtSplusB, SigmaMeanOverMean, mean, 
  // set only some of them with h->GetBinContent(h->GetMinimumBin()) and h->GetBinContent(h->GetMaximumBin())
  th2dMinZaxisVector.push_back(0.0);  // was set to hSignal->GetBinContent(hSignal->GetMinimumBin())
  th2dMinZaxisVector.push_back(0.0);  // hBackground->GetBinContent(hBackground->GetMinimumBin())
  th2dMinZaxisVector.push_back(0.0);
  th2dMinZaxisVector.push_back(0.0);
  th2dMinZaxisVector.push_back(0.001);//0.0
  th2dMinZaxisVector.push_back(0.130);
  th2dMinZaxisVector.push_back(0.005);

  Double_t normalizedS = 0.0;
  Double_t normalizedB = 0.0;

  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    
    if (jentry % 100000 == 0) cout << jentry << endl;

    if ((abs(Backgr) > 0.00001) && (abs(Signal) > 0.00001)) { // avoid empty crystals due to masked towers or whatever
     
      normalizedS = Signal * fit_Snorm;
      normalizedB = Backgr * fit_Bnorm;

      // to avoid that in 2D maps points below lower threshold in z axis are drawn white (as if they are empty), fill with the maximum between threshold and value

      hSignal->Fill((Double_t)iphi,(Double_t)ieta,max(th2dMinZaxisVector[0],(Double_t)normalizedS));
      hBackground->Fill((Double_t)iphi,(Double_t)ieta,max(th2dMinZaxisVector[1],(Double_t)normalizedB));
      SoverB->Fill((Double_t)iphi,(Double_t)ieta,max(th2dMinZaxisVector[2],(Double_t)normalizedS/normalizedB));
      SoverSqrtSplusB->Fill((Double_t)iphi,(Double_t)ieta, max(th2dMinZaxisVector[3],(Double_t)normalizedS/sqrt(normalizedS + normalizedB)));
      SigmaMeanOverMean->Fill((Double_t)iphi,(Double_t)ieta, max(th2dMinZaxisVector[4],(Double_t)fit_mean_err/fit_mean));
      mean->Fill((Double_t)iphi,(Double_t)ieta, max(th2dMinZaxisVector[5],(Double_t)fit_mean));
      sigma->Fill((Double_t)iphi,(Double_t)ieta, max(th2dMinZaxisVector[6],(Double_t)fit_sigma));

      // hSignal_ietaProfile->Fill((Double_t)ieta,normalizedS);
      // hBackground_ietaProfile->Fill((Double_t)ieta,normalizedB);
      // SoverB_ietaProfile->Fill((Double_t)ieta,normalizedS/normalizedB);
      // SoverSqrtSplusB_ietaProfile->Fill((Double_t)ieta, normalizedS/sqrt(normalizedS + normalizedB));
      // SigmaMeanOverMean_ietaProfile->Fill((Double_t)ieta, fit_mean_err/fit_mean);
      // mean_ietaProfile->Fill((Double_t)ieta, fit_mean);
      // sigma_ietaProfile->Fill((Double_t)ieta, fit_sigma);
 
      Double_t eta = ieta * 0.0174532925; // 0.0174532925 = pi/180; in the barrel eta is simply proportional to ieta 

      hSignal_etaProfile->Fill((Double_t)eta,normalizedS);
      hBackground_etaProfile->Fill((Double_t)eta,normalizedB);
      SoverB_etaProfile->Fill((Double_t)eta,normalizedS/normalizedB);
      SoverSqrtSplusB_etaProfile->Fill((Double_t)eta, normalizedS/sqrt(normalizedS + normalizedB));
      SigmaMeanOverMean_etaProfile->Fill((Double_t)eta, fit_mean_err/fit_mean);
      mean_etaProfile->Fill((Double_t)eta, fit_mean);
      sigma_etaProfile->Fill((Double_t)eta, fit_sigma);
 
      mean_iphiProfile->Fill((Double_t)iphi, fit_mean);

   }    

  }
 
  // for (UInt_t i = 0; i < th2dVector.size(); i++) {

  //   // set range to min and max, regardless the fact that values are overriden (for instance, in h goes from 0 to 10 and I tried to set the scale in 3-8, 3 and 8 would become the values returned by GetMinimum() and GetMaximum(). Here we want the real minimum and maximum).
  //   th2dMinZaxisVector.push_back(th2dVector[i]->GetBinContent(th2dVector[i]->GetMinimumBin()));
  //   th2dMaxZaxisVector.push_back(th2dVector[i]->GetBinContent(th2dVector[i]->GetMaximumBin()));

  // }


  th2dMaxZaxisVector.push_back(hSignal->GetBinContent(hSignal->GetMaximumBin()));
  th2dMaxZaxisVector.push_back(hBackground->GetBinContent(hBackground->GetMaximumBin()));
  th2dMaxZaxisVector.push_back(10e9);
  th2dMaxZaxisVector.push_back(10e9); // when this value is very large (bigger than the default) use the default to plot axis
  th2dMaxZaxisVector.push_back(0.0125);//0.02
  th2dMaxZaxisVector.push_back(0.140);
  th2dMaxZaxisVector.push_back(0.015);

  for (UInt_t i = 0; i < th2dVector.size(); i++) {

    if (th2dMinZaxisVector[i] < 0.0) th2dMinZaxisVector[i] = 0.0;
    th2dVector[i]->SetMinimum(th2dMinZaxisVector[i]);                                                                                                               
    //profileIetaVector[i]->SetMinimum(th2dMinZaxisVector[i]);
    profileEtaVector[i]->SetMinimum(th2dMinZaxisVector[i]);
 
    // if the maximum choosen by user is bigger than default, don't do anything, otherwise set the user value as the maximum                                        
    if (th2dMaxZaxisVector[i] < th2dVector[i]->GetBinContent(th2dVector[i]->GetMaximumBin())) {

      th2dVector[i]->SetMaximum(th2dMaxZaxisVector[i]);
      //profileIetaVector[i]->SetMaximum(th2dMaxZaxisVector[i]);
      profileEtaVector[i]->SetMaximum(th2dMaxZaxisVector[i]);

    }    

  }

  string wwwPath = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/" + dirName + "/" + iterNumber + "/2DMaps/Barrel/";
  //  string wwwPath = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/AlcaP0_fromRun273158_2016_v2/7thIteration/2DMaps/Barrel/test/";

  vector<TCanvas *> c;
  //vector<TCanvas *> cprof;
  vector<TCanvas *> cprof_eta;

  for ( UInt_t i = 0; i < th2dVector.size(); i++ ) {

    c.push_back(NULL);
    string canvasName(th2dVector[i]->GetName());
    canvasName = "c_" + canvasName;
    c[i] = new TCanvas(canvasName.c_str(),canvasName.c_str());
    th2dVector[i]->Draw("COLZ");
    th2dVector[i]->GetXaxis()->SetTitle("i #phi");
    th2dVector[i]->GetXaxis()->SetTitleSize(0.06);
    th2dVector[i]->GetXaxis()->SetTitleOffset(0.7);
    th2dVector[i]->GetYaxis()->SetTitle("i #eta");
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
    string name = wwwPath + th2dVector[i]->GetName() + "_EB";
    c[i]->SaveAs((name + ".pdf").c_str());
    c[i]->SaveAs((name + ".png").c_str());

    // NOW DRAWING PROFILE 
    
    // cprof.push_back(NULL);
    // canvasName = "c_" + string(profileIetaVector[i]->GetName()); 
    // cprof[i] = new TCanvas(canvasName.c_str(),canvasName.c_str());
    // profileIetaVector[i]->Draw("E");
    // profileIetaVector[i]->GetXaxis()->SetTitle("i #eta");
    // profileIetaVector[i]->GetXaxis()->SetTitleSize(0.06);
    // profileIetaVector[i]->GetXaxis()->SetTitleOffset(0.7);
    // // profileIetaVector[i]->GetYaxis()->SetTitle("");
    // // profileIetaVector[i]->GetYaxis()->SetTitleSize(0.06);
    // // profileIetaVector[i]->GetYaxis()->SetTitleOffset(0.8);
    // profileIetaVector[i]->SetStats(0);
    // profileIetaVector[i]->Draw("E");
    // name = wwwPath + profileIetaVector[i]->GetName() + "_EB";
    // cprof[i]->SaveAs((name + ".pdf").c_str());
    // cprof[i]->SaveAs((name + ".png").c_str());

    cprof_eta.push_back(NULL);
    canvasName = "c_" + string(profileEtaVector[i]->GetName()); 
    cprof_eta[i] = new TCanvas(canvasName.c_str(),canvasName.c_str());
    profileEtaVector[i]->Draw("HE");
    profileEtaVector[i]->GetXaxis()->SetTitle("#eta");
    profileEtaVector[i]->GetXaxis()->SetTitleSize(0.06);
    profileEtaVector[i]->GetXaxis()->SetTitleOffset(0.7);
    profileEtaVector[i]->GetYaxis()->SetTitle( profileYaxisTitle[i].c_str() );
    profileEtaVector[i]->GetYaxis()->SetTitleSize(0.055);
    profileEtaVector[i]->GetYaxis()->SetTitleOffset(0.8);
    profileEtaVector[i]->SetStats(0);
    profileEtaVector[i]->Draw("HE");
    name = wwwPath + profileEtaVector[i]->GetName() + "_EB";
    cprof_eta[i]->SaveAs((name + ".pdf").c_str());
    cprof_eta[i]->SaveAs((name + ".png").c_str());
    
  }  

  TCanvas *c_mean_iphi_prof = new TCanvas("c_mean_iphi_prof",("c_" + string(mean_iphiProfile->GetName())).c_str());
  mean_iphiProfile->Draw("HE");
  mean_iphiProfile->GetXaxis()->SetTitle("i #phi");
  mean_iphiProfile->GetXaxis()->SetTitleSize(0.06);
  mean_iphiProfile->GetXaxis()->SetTitleOffset(0.7);
  mean_iphiProfile->GetYaxis()->SetTitle("mean [GeV]");
  mean_iphiProfile->GetYaxis()->SetRangeUser(0.13,0.14);
  mean_iphiProfile->GetYaxis()->SetTitleSize(0.055);
  mean_iphiProfile->GetYaxis()->SetTitleOffset(0.8);  
  mean_iphiProfile->SetStats(0);
  mean_iphiProfile->Draw("HE");
  c_mean_iphi_prof->SaveAs( (wwwPath + mean_iphiProfile->GetName() + "_EB.pdf").c_str() );
  c_mean_iphi_prof->SaveAs( (wwwPath + mean_iphiProfile->GetName() + "_EB.png").c_str() );

}


Int_t main(int argc, char* argv[]) {

  TChain *chain = new TChain("calibEB");
 
  string fileToChain = "";
  string path(argv[1]);
  string eosPath = "root://eoscms//eos/cms" + path;
  string dirName(argv[2]);
  string iter_num(argv[3]);
  string tagName(argv[4]);

  // string eosPath = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/";
  // string dirName = "AlcaP0_fromRun273158_2016_v2";
  // string iter_num = "6"; 

  string Result;          // string which will contain the result  

  for (Int_t i = 0; i < 31; i++) {

    ostringstream convert;   // stream used for the conversion                                                                                                         
    convert << i;      // insert the textual representation of 'i' in the characters in the stream                                                                     
    Result = convert.str();
    fileToChain = eosPath + dirName + "/iter_" + iter_num + "/" + tagName + "Barrel_" + Result + "_calibMap.root";
    chain->Add(TString(fileToChain.c_str()));

  }

  
  // for some reason the already merged file doesn't work when running on the chain. So I use the TChain above --> reason is that variables end with _, they are different
  //  fileToChain = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2016/mciprian/AlcaP0_fromRun273158_2016_v2/iter_6/AlcaP0_fromRun273158_2016_v2_calibMap.root";
  //chain->Add(TString(fileToChain.c_str()));

  string iterNumber = "iter_" + iter_num;
  calibAnaEB *ana = new calibAnaEB(chain);
  ana->Loop(dirName, iterNumber);

  return 0;


}


