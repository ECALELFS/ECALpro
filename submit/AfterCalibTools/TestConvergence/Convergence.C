#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <set>
#include <cmath>
#include <cstdlib>
//Root Stuff
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::map;
using std::pair;
using std::stringstream;

static const int MAX_IETA = 85;
static const int MAX_IPHI = 360;
static const int MIN_IETA = 1;
static const int MIN_IPHI = 1;

//Usage: .x Convergence.C+("/store/caf/user/lpernie/","ALL_2010ArelevantFiles_NOXTALWRONG_01/",14,"")
//Usage: .x Convergence.C+("/store/group/alca_ecalcalib/lpernie/","ALL_2010_WithNEWSelection_01",6,"2012Cmerg_")
void Convergence( string Path_0, string Path, int nIter, string Tag ){

  string PathL = "root://eoscms//eos/cms" + Path_0 + Path;

  system( (string("mkdir -p plot_") + Path ).c_str());

  TCanvas* myc1 = new TCanvas("myc1", "CMS", 600, 600);
  for(int isEB=0; isEB<2; isEB++){  

    float *EB_RMS = NULL;
    EB_RMS = new float[nIter];
    float *iter= NULL;
    iter = new float[nIter];

    float hmean(0.), hrms(0.01);
    if(isEB==1){ hmean=0.02; hrms=0.015; }

    for(int i=0; i<nIter-1; i++){

      //Iter
      stringstream ss; ss<<i;
      stringstream ss1; ss1<<(i+1);
      string Iter = ss.str();
      string Iter1 = ss1.str();
      // Input
      string fileName = string(PathL) + "/iter_" + string(Iter) + "/" + string(Tag) + "calibMap.root";   
      cout<<"Opening: "<<fileName<<endl;
      TFile* fout = TFile::Open(fileName.c_str());
      fileName = string(PathL) + "/iter_" + string(Iter1) + "/" + string(Tag) +"calibMap.root";   
      cout<<"And: "<<fileName<<endl;
      TFile* fout1 = TFile::Open(fileName.c_str());

      TTree *Tree; TTree *Tree1;
      if(isEB==0){
           Tree  = (TTree*) fout->Get("calibEB");
           Tree1 = (TTree*) fout1->Get("calibEB");
      }
      if(isEB==1){
           Tree  = (TTree*) fout->Get("calibEE");
           Tree1 = (TTree*) fout1->Get("calibEE");
      }
      Float_t coeff, coeff1;
      Float_t Ndof, Ndof1;
      Tree->SetBranchAddress( "coeff_", &coeff);
      Tree1->SetBranchAddress( "coeff_", &coeff1);
      Tree->SetBranchAddress( "Ndof_", &Ndof);
      Tree1->SetBranchAddress( "Ndof_", &Ndof1);

      //Histo
      TH1F *h1; h1 =new TH1F("h1","",100,hmean-3*hrms,hmean+3*hrms);

      //Loop
      Long64_t nentries = Tree->GetEntriesFast();
      for(Long64_t iEntry=0; iEntry<nentries; iEntry++){
         Tree->GetEntry(iEntry);
         Tree1->GetEntry(iEntry);
         if( coeff1!=1. && coeff!=1. && coeff1!=coeff && coeff!=0 && Ndof>10 && Ndof1>10){
             if(isEB==0 )                              h1->Fill((coeff1-coeff));
             if(isEB==1 && coeff>0.97 && coeff1>0.97 ) h1->Fill((coeff1-coeff));
         }
      }
      h1->Draw();
      TString out;
      if(isEB==0) out = "plot_" + Path + "/EB_Iter_" + Iter + ".png";
      if(isEB==1) out = "plot_" + Path + "/EE_Iter_" + Iter + ".png";
      myc1->SaveAs(out.Data());
      EB_RMS[i]=h1->GetRMS();
      iter[i]=i+1;
      if(i<11){
         hmean = h1->GetMean();
         hrms  = h1->GetRMS();
      }
    }
    TGraph *Conv = new TGraph(nIter-1, iter, EB_RMS);
    Conv->SetLineColor(2);
    Conv->SetLineWidth(1);
    Conv->SetMarkerColor(2);
    Conv->SetMarkerStyle(20);
    Conv->SetMarkerSize(0.5);
    if(isEB==0) Conv->SetTitle("EB) IC Convergence");
    if(isEB==1) Conv->SetTitle("EE) IC Convergence");
    Conv->GetXaxis()->SetTitle("Iter");
    Conv->GetYaxis()->SetTitle("RMS(ICn+1 - IC)");
    Conv->Draw("ACP");
    myc1->cd();
    TString out;
    if(isEB==0) out = "plot_" + Path + "/EB_IC_Convergence.png";
    if(isEB==1) out = "plot_" + Path + "/EE_IC_Convergence.png";
    myc1->SaveAs(out.Data());
  }
}
