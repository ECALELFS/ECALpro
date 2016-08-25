#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooChebychev.h"

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
#include "TLatex.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"
#include "TFormula.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"

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

//================================================

// function used in Convergence()
// given the file endcap_ix_iy_zside_ietaRing.dat, it takes ix, iy and zside and returns the corresponding etaRing index 

// N.B.: this function is very slow !!! Will need to define a faster way, maybe save eta-Ring index in trees?

class Convergence {

public:
  Convergence( string Path_0, string Path, int nIter, string Tag, int nJump=1 );
  ~Convergence() {};

  void addExtension( string Path, int nIter, string Tag, int nJump=1 );
  void run();

private:
  Int_t getEtaRingInEE(Int_t &ix, Int_t &iy, Int_t &zside);
  
  string Path_0_;
  std::vector<string> Paths_, Tags_;
  std::vector<int> nIters_, nJumps_;

};
  
Convergence::Convergence( string Path_0, string Path, int nIter, string Tag, int nJump) {
  Path_0_ = Path_0;
  Paths_.push_back(Path);
  Tags_.push_back(Tag);
  nIters_.push_back(nIter);
  nJumps_.push_back(nJump);
}

void Convergence::addExtension( string Path, int nIter, string Tag, int nJump) {
  Paths_.push_back(Path);
  Tags_.push_back(Tag);
  nIters_.push_back(nIter);
  nJumps_.push_back(nJump);
}

Int_t Convergence::getEtaRingInEE(Int_t &ix, Int_t &iy, Int_t &zside) {

  Int_t etaRing = -1;

  string fileName = "endcap_ix_iy_zside_ietaRing.dat";
  ifstream inputFile(fileName.c_str());

  Int_t a,b,c,d;  // file format is --> a b c d                             
 
  if (inputFile.is_open()) {

    while ((etaRing == -1) && (inputFile >> a >> b >> c >> d)) {
      // first check zside, because it will remove half of the lines to check         
      if (c == zside && a == ix && b == iy) etaRing = d;
    }

  } else {
    std::cout << "Error: could not open file " << fileName << std::endl;
    exit(EXIT_FAILURE);
  }

  inputFile.close();

  return etaRing;

}

//=================================================

//5_3_6:  gROOT->ProcessLine(".include /afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms9/include/")
//Usage: .x Convergence.C+("/store/group/dpg_ecal/alca_ecalcalib/lpernie/","ALL_2015B_Multifit_01",13,"2015B_")
void Convergence::run() {

  system( (string("mkdir -p plot_") + Paths_[0] ).c_str());
  TCanvas* myc1 = new TCanvas("myc1", "CMS", 600, 600);
  TString outname = "plot_" + Paths_[0] + "/Differences.root";
  TFile* output = new TFile(outname.Data(),"RECREATE");
  TH2F* rms_EB  = new TH2F("rms_EB","IC(n)-IC(n-1) #phi on x #eta on y",MAX_IPHI, MIN_IPHI, MAX_IPHI, 2*MAX_IETA+1, -MAX_IETA-0.5, MAX_IETA+0.5 );
  TH2F* rms_EEp = new TH2F("rms_EEp","IC(n)-IC(n-1) iX on x iY on y (EEp)",100,0.5,100.5,100,0.5,100.5);
  TH2F* rms_EEm = new TH2F("rms_EEm","IC(n)-IC(n-1) iY on x iY on y (EEm)",100,0.5,100.5,100,0.5,100.5);

  vector<Int_t> etaRingEdges;  

  int nIter = -1; // lookig at differences: steps are N-1
  for(int c=0; c<(int)Paths_.size(); ++c) nIter += nIters_[c];

  std::cout << "====> Test of the convergence of the calibrations. <====" << std::endl
            << "Will run on " << Paths_.size() << " chunks of calibrations, for a total of "
            << nIter << " iterations" << std::endl;
  
  for(int isEB=0; isEB<2; isEB++){

    etaRingEdges.clear();  // erase all elements, as if it was created here (will be filled differently for EB and EE)
    
    float *EB_RMS = NULL;
    EB_RMS = new float[nIter];
    float *iter= NULL;
    iter = new float[nIter];
    
    float hmean(0.), hrms(0.01), sigma_plot(0);
    float hrange;
    int nbins;
    if(isEB==0){ // barrel
      // divide EB in bins of |ieta|, from  0 up to 86 (max(ieta) = 85, so ieta=86 is not considered in the range)
      //First range will be [1,20) and so on (upper value excluded)
      etaRingEdges.push_back(1);
      etaRingEdges.push_back(26);
      etaRingEdges.push_back(46);
      etaRingEdges.push_back(66);
      etaRingEdges.push_back(86);
      
      hmean=0.09; hrms=0.03; 
      hrange = 0.05;
      nbins = 200;
    } 
    else if (isEB==1){
      
      // divide EE in bins of ietaRing, from 0 up to 37 (max(ietaRing index) = 36, so etaRing=37 is not considered in the range). 
      //First range will be [1,9) and so on (upper value excluded) 
      etaRingEdges.push_back(1);
      etaRingEdges.push_back(10);
      etaRingEdges.push_back(19);
      etaRingEdges.push_back(28);
      etaRingEdges.push_back(37);
      
      hrange = 0.10; 
      nbins = 50;
    }

    Int_t n_hbinned;
    // Note: n edges --> n-1 bins
    n_hbinned = etaRingEdges.size() -1;		
    Float_t EB_RMS_etaRing[nIter][n_hbinned];
    
    int PrevChunkConsidered = false;
    int iterOffset = 0;
    for(int iChunk=0; iChunk<(int)Paths_.size(); ++iChunk) {

      std::cout << "Running chunk " << iChunk << " iteration offset = " << iterOffset << std::endl;

      string PathL = "root://eoscms//eos/cms" + Path_0_ + Paths_[iChunk];
      
      for(int i=0; i<(int)nIters_[iChunk]-1; ++i){

        //Iter
        stringstream ss; ss<<i;
        stringstream ss1; ss1<<(i+nJumps_[iChunk]);
        stringstream ssTot; ssTot<<(i+iterOffset);
        stringstream ssTot1; ssTot1<<(i+iterOffset+nJumps_[iChunk]);
        string Iter = ss.str();
        string Iter1 = ss1.str();
        string IterTot = ssTot.str();
        string IterTot1 = ssTot1.str();
        // Input
        TFile *fout=0;
        TFile *fout1=0;
        if(iChunk>0 && i==0 && !PrevChunkConsidered) {
          string PathLPrevChunk = "root://eoscms//eos/cms" + Path_0_ + Paths_[iChunk-1];
          stringstream ssNminus1; ssNminus1 << nIters_[iChunk-1]-1;
          string fileName = string(PathLPrevChunk) + "/iter_" + ssNminus1.str() + "/" + string(Tags_[iChunk-1]) + "calibMap.root";   
          cout<<"Opening: "<<fileName<<endl;
          fout = TFile::Open(fileName.c_str());
          fileName = string(PathL) + "/iter_" + string(Iter) + "/" + string(Tags_[iChunk]) +"calibMap.root";   
          cout<<"And: "<<fileName<<endl;
          fout1 = TFile::Open(fileName.c_str());
          PrevChunkConsidered = true;
          i=-1; iterOffset++;
        } else {
          string fileName = string(PathL) + "/iter_" + string(Iter) + "/" + string(Tags_[iChunk]) + "calibMap.root";   
          cout<<"Opening: "<<fileName<<endl;
          fout = TFile::Open(fileName.c_str());
          fileName = string(PathL) + "/iter_" + string(Iter1) + "/" + string(Tags_[iChunk]) +"calibMap.root";   
          cout<<"And: "<<fileName<<endl;
          fout1 = TFile::Open(fileName.c_str());
        }
        
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
        Int_t ieta, iphi, ix, iy, iz;
        Tree->SetBranchAddress( "coeff_", &coeff);
        Tree1->SetBranchAddress( "coeff_", &coeff1);
        Tree->SetBranchAddress( "Ndof_", &Ndof);
        Tree1->SetBranchAddress( "Ndof_", &Ndof1);
        Tree->SetBranchAddress( "Ndof_", &Ndof);
        Tree1->SetBranchAddress( "Ndof_", &Ndof1);
        if(isEB==0){
          Tree1->SetBranchAddress( "ieta_", &ieta);
          Tree1->SetBranchAddress( "iphi_", &iphi);
        }
        if(isEB==1){
          Tree1->SetBranchAddress( "ix_", &ix);
          Tree1->SetBranchAddress( "iy_", &iy);
          Tree1->SetBranchAddress( "zside_", &iz);
        }
        
        //Histo
        string hXaxisName = "IC_{" + IterTot1 + "}-IC_{" + IterTot + "}";  // e.g. IC_{1}-IC_{0} if nJump==1, or IC_{2}-IC_{0} if nJump==2...
        
        //TH1F *h1; h1 =new TH1F("h1","",1000,hmean-9*hrms,hmean+9*hrms);
        TH1F *h1; h1 =new TH1F("h1","",nbins,-1*hrange,hrange);
        h1->GetXaxis()->SetTitle(hXaxisName.c_str());
        
        TH1F *h_etaRing[n_hbinned];
        
        for (Int_t k = 0; k < n_hbinned; k++) {
          // create histograms and set title in printf() style using Form()
          h_etaRing[k] = new TH1F(Form("h_etaRing_%dTo%d",etaRingEdges[k],etaRingEdges[k+1]-1),Form("#eta-Ring %d To %d",etaRingEdges[k],etaRingEdges[k+1]-1),nbins,-1*hrange,hrange);
          h_etaRing[k]->GetXaxis()->SetTitle(hXaxisName.c_str());
        } 
        
        //Loop
        Long64_t nentries = Tree->GetEntriesFast();
        for(Long64_t iEntry=0; iEntry<nentries; iEntry++){
          
          Tree->GetEntry(iEntry);
          Tree1->GetEntry(iEntry);
          
          if( coeff1!=1. && coeff!=1. && coeff1!=coeff && coeff!=0 && coeff1!=0 /*&& Ndof>10 && Ndof1>10*/){
            
            h1->Fill((coeff1-coeff));
            // now fill histograms for different  ietaRing
            Int_t binIndex = 0;
            Int_t binFound = 0;
            
            while (!binFound && (binIndex < n_hbinned)) {
              if (isEB == 0) {
                if (abs(ieta) >= etaRingEdges[binIndex] && abs(ieta) < etaRingEdges[binIndex+1]) {
                  h_etaRing[binIndex]->Fill((coeff1-coeff));
                  binFound = 1; // get out of loop when bin is found
                }
              } else if (isEB == 1) {
                Int_t etaRing = getEtaRingInEE(ix,iy,iz);  // this function is very slow!!!
                if (etaRing >= etaRingEdges[binIndex] && etaRing < etaRingEdges[binIndex+1]) {
                  h_etaRing[binIndex]->Fill((coeff1-coeff));
                  binFound = 1; // get out of loop when bin is found
                }
              }
              binIndex++;
            }
            
          }
          
          if(i==nIters_[iChunk]-1){
            if(isEB==0 && coeff1!=1. && coeff!=1. && coeff1!=coeff && coeff!=0 /*&& Ndof>10 && Ndof1>10*/){ 
              rms_EB->SetBinContent(iphi, ieta+86., fabs(coeff1-coeff)/coeff1);
            }
            if(isEB==1){
              if(iz==1){    
                rms_EEp->SetBinContent(ix, iy, fabs(coeff1-coeff)/coeff1);
              }
              else if(iz==-1)rms_EEm->SetBinContent(ix, iy, fabs(coeff1-coeff)/coeff1);
              else cout<<"WARNING!!! zside_ not -1 or 1"<<endl;
            }
          }
          
        }
        gStyle->SetOptStat(111111);
        TString out;
        h1->Draw();
        if(isEB==0) out = "plot_" + Paths_[0] + "/EB_Iter_" + IterTot + ".png";
        else if(isEB==1) out = "plot_" + Paths_[0] + "/EE_Iter_" + IterTot + ".png";
        
        h1->Draw("hist");
        
        // TLatex lat;
        // char line[300];
        // lat.SetNDC();
        // lat.SetTextSize(0.030);
        // lat.SetTextColor(1);
        // sprintf(line,"SIGMA: %.6f ", h1->GetRMS() );
        // float xmin(0.55), yhi(0.80);// ypass(0.05);
        // lat.DrawLatex(xmin,yhi, line);
        
        myc1->SaveAs(out.Data());

        hmean=h1->GetMean();
        EB_RMS[i+iterOffset]=h1->GetRMS();
        iter[i+iterOffset]=i+iterOffset+1;

        for (Int_t k = 0; k < n_hbinned; k++) {
          
          EB_RMS_etaRing[i+iterOffset][k] = h_etaRing[k]->GetRMS();
          
          stringstream ssLow;
          stringstream ssUp;
          ssLow<<etaRingEdges[k];
          ssUp<<(etaRingEdges[k+1]-1);  // if edges are 1, 10, 19 ... we want first bin from 1 to 9 (included), then from 10 to 18 (included) and so on
          string etaRingLow = ssLow.str();
          string etaRingUp = ssUp.str();

          if(isEB==0) out = "plot_" + Paths_[0] + "/EB_Iter_" + IterTot + "_etaRing" + etaRingLow + "To" + etaRingUp + ".png";
          else if(isEB==1) out = "plot_" + Paths_[0] + "/EE_Iter_" + IterTot + "_etaRing" + etaRingLow + "To" + etaRingUp + ".png";
          h_etaRing[k]->Draw("hist");
          myc1->SaveAs(out.Data());
          
          delete h_etaRing[k];  // delete histogram before new iteration starts
          
        }        
      }
      iterOffset += (nIters_[iChunk]-1);
    }

    TGraph *Conv = new TGraph(nIter, iter, EB_RMS);
    Conv->SetLineColor(2);
    Conv->SetLineWidth(1);
    Conv->SetMarkerColor(2);
    Conv->SetMarkerStyle(20);
    Conv->SetMarkerSize(0.5);
    if(isEB==0) Conv->SetTitle("EB: IC Convergence");
    if(isEB==1) Conv->SetTitle("EE: IC Convergence");
    Conv->GetXaxis()->SetTitle("Iteration");
    //Conv->GetYaxis()->SetOffset(1.);
    //if(nJump==1) Conv->GetYaxis()->SetTitle("RMS(ICn+1 - IC)");
    //if(nJump==2) Conv->GetYaxis()->SetTitle("RMS(ICn+2 - IC)");
    if(nJumps_[0]==1) Conv->GetYaxis()->SetTitle("RMS(IC_n - IC_{n-1})");  // because X axis starts from 1, so we have RMS(IC_1 - IC_0) and so on 
    if(nJumps_[0]==2) Conv->GetYaxis()->SetTitle("RMS(IC_n - IC_{n-2})");  // because X axis starts from 2, so we have RMS(IC_2 - IC_0) and so on
    Conv->Draw("ACP");
    myc1->cd();
    TString out;
    if(isEB==0) out = "plot_" + Paths_[0] + "/EB_IC_Convergence.png";
    if(isEB==1) out = "plot_" + Paths_[0] + "/EE_IC_Convergence.png";
    myc1->SaveAs(out.Data());
    
    // now get proper y value from EB_RMS_etaRing[i][k] (need values at constant k), set points for TGraph and draw again
    for (Int_t k = 0; k < n_hbinned; k++) {
      
      stringstream ssLow;
      stringstream ssUp;
      ssLow<<etaRingEdges[k];
      ssUp<<(etaRingEdges[k+1]-1);  // if edges are 1, 10, 19 ... we want first bin from 1 to 9 (included), then from 10 to 18 (included) and so on
      string etaRingLow = ssLow.str();
      string etaRingUp = ssUp.str();
      
      for (Int_t iterIndex = 0; iterIndex < nIter; iterIndex++) {
        Conv->SetPoint(iterIndex, iter[iterIndex], EB_RMS_etaRing[iterIndex][k]);
      }	    
      Conv->Draw("ACP");
      if(isEB==0) out = "plot_" + Paths_[0] + "/EB_IC_Convergence_etaRing" + etaRingLow + "To" + etaRingUp + ".png";
      if(isEB==1) out = "plot_" + Paths_[0] + "/EE_IC_Convergence_etaRing" + etaRingLow + "To" + etaRingUp + ".png";
      myc1->SaveAs(out.Data());
    }
    
    
  }
  output->cd();
  rms_EB->Write();
  rms_EEp->Write();
  rms_EEm->Write();
  output->Close();
    
}
