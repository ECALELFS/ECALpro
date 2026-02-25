#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
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

using namespace RooFit;

//Usage: .x CheckFit.C+("root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/lpernie/ALL_2015A_RAW_RECHIT_SMIC_estimTime_01/iter_0/2015A_calibMap.root",
//                       "PlotFIT2015A_pi0", 0, true)
//EBorEE: 0(EB+EE) 1(EB) 2(EE);
void CheckFit( TString File, TString folder , int EBorEE=0, bool Are_pi0_=true ){

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit(111110);
  gStyle->SetOptFile(1); 
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.3);
  gStyle->SetMarkerColor(1);
  TCanvas* myc1 = new TCanvas("myc1"," ",500,500);
  myc1->SetLeftMargin(0.16);
  TString Comm = "mkdir -p " + folder;
  system( Comm.Data() );
  //File, trees and variables
  ofstream myfileEB, myfileEE, myfileEE_range;
  myfileEB.open (folder + "/NotFittedEB.txt");
  myfileEE.open (folder + "/NotFittedEE.txt");
  myfileEE_range.open (folder + "/BadRangeEE.txt");
  cout<<"Opening: "<<File.Data()<<endl;
  TFile* fin = TFile::Open(File.Data());
  TTree *treeEB = (TTree*) fin->Get("calibEB");
  TTree *treeEE = (TTree*) fin->Get("calibEE");
  Float_t EBfit_b0_, EBfit_b1_, EBfit_b2_, EBfit_b3_, EBfit_Bnorm_, EBfit_mean_, EBfit_sigma_, EBfit_Snorm_, EBChisqu_, EBNdof_, EBcoeff_;
  Int_t EBieta_, EBiphi_, EBhashedIndex_;
  treeEB->SetBranchAddress("fit_b0_",&EBfit_b0_);
  treeEB->SetBranchAddress("fit_b1_",&EBfit_b1_);
  treeEB->SetBranchAddress("fit_b2_",&EBfit_b2_);
  treeEB->SetBranchAddress("fit_b3_",&EBfit_b3_);
  treeEB->SetBranchAddress("fit_Bnorm_",&EBfit_Bnorm_);
  treeEB->SetBranchAddress("fit_mean_",&EBfit_mean_);
  treeEB->SetBranchAddress("fit_sigma_",&EBfit_sigma_);
  treeEB->SetBranchAddress("fit_Snorm_",&EBfit_Snorm_);
  treeEB->SetBranchAddress("ieta_",&EBieta_);
  treeEB->SetBranchAddress("iphi_",&EBiphi_);
  treeEB->SetBranchAddress("Chisqu_",&EBChisqu_);
  treeEB->SetBranchAddress("Ndof_",&EBNdof_);
  treeEB->SetBranchAddress("hashedIndex_",&EBhashedIndex_);
  treeEB->SetBranchAddress("coeff_",&EBcoeff_);
  Float_t EEfit_b0_, EEfit_b1_, EEfit_b2_, EEfit_b3_, EEfit_Bnorm_, EEfit_mean_, EEfit_sigma_, EEfit_Snorm_, EEChisqu_, EENdof_, EEcoeff_;
  Int_t EEix_, EEiy_, EEzside_, EEhashedIndex_;
  treeEE->SetBranchAddress("fit_b0_",&EEfit_b0_);
  treeEE->SetBranchAddress("fit_b1_",&EEfit_b1_);
  treeEE->SetBranchAddress("fit_b2_",&EEfit_b2_);
  treeEE->SetBranchAddress("fit_b3_",&EEfit_b3_);
  treeEE->SetBranchAddress("fit_Bnorm_",&EEfit_Bnorm_);
  treeEE->SetBranchAddress("fit_mean_",&EEfit_mean_);
  treeEE->SetBranchAddress("fit_sigma_",&EEfit_sigma_);
  treeEE->SetBranchAddress("fit_Snorm_",&EEfit_Snorm_);
  treeEE->SetBranchAddress("ix_",&EEix_);
  treeEE->SetBranchAddress("iy_",&EEiy_);
  treeEE->SetBranchAddress("zside_",&EEzside_);
  treeEE->SetBranchAddress("Chisqu_",&EEChisqu_);
  treeEE->SetBranchAddress("Ndof_",&EENdof_);
  treeEE->SetBranchAddress("hashedIndex_",&EEhashedIndex_);
  treeEE->SetBranchAddress("coeff_",&EEcoeff_);

  //Loop for EB and EE
  std::vector<int> IndexBadEB, Coo1BadEB, Coo2BadEB, Coo3BadEB; IndexBadEB.clear(); Coo1BadEB.clear(); Coo2BadEB.clear(); Coo3BadEB.clear();
  std::vector<int> IndexBadEE, Coo1BadEE, Coo2BadEE, Coo3BadEE; IndexBadEE.clear(); Coo1BadEE.clear(); Coo2BadEE.clear(); Coo3BadEE.clear();
  std::vector<int> IndexBadRangeEE, Coo1BadRangeEE, Coo2BadRangeEE, Coo3BadRangeEE; IndexBadRangeEE.clear(); Coo1BadRangeEE.clear(); Coo2BadRangeEE.clear(); Coo3BadRangeEE.clear();
  TH2F *NoFit_EB    = new TH2F("NoFit_EB","NoFit_EB",2*MAX_IETA+1,-MAX_IETA-0.5,MAX_IETA+0.5,MAX_IPHI, MIN_IPHI-0.5, MAX_IPHI+0.5 );
  TH2F *NoFit_EEp   = new TH2F("NoFit_EEp","NoFit_EEp",101,-0.5,100.5,101,-0.5,100.5);
  TH2F *NoFit_EEm   = new TH2F("NoFit_EEm","NoFit_EEm",101,-0.5,100.5,101,-0.5,100.5);
  TH2F *BadRangeEEp = new TH2F("BadRangeEEp","BadRangeEEp",101,-0.5,100.5,101,-0.5,100.5);
  TH2F *BadRangeEEm = new TH2F("BadRangeEEm","BadRangeEEm",101,-0.5,100.5,101,-0.5,100.5);
  for( int EBEE=0; EBEE<2; EBEE++ ){
    if( (EBEE==0 && EBorEE==2) || (EBEE==1 && EBorEE==1) ) continue;
    //Change things to run on EB and EE
    TString nameEBEE = "";
    TTree *tree=0;
    if(EBEE==0){
	nameEBEE = "EB";
	tree = treeEB;
    }
    if(EBEE==1){
	nameEBEE = "EE";
	tree = treeEE;
    }
    Int_t eventN = tree->GetEntries();
    //fastCheck Chi2
    tree->Draw("Chisqu_");
    TString Hname = folder + "/" + nameEBEE + "_";
    if(Are_pi0_)  Hname += "pi0_Chi2.png";
    else          Hname += "eta_Chi2.png";
    myc1->SaveAs( Hname.Data() );
    tree->Draw("Chisqu_","Chisqu_<5");
    myc1->Update();
    Hname = folder + "/" + nameEBEE + "_";
    if(Are_pi0_)  Hname += "pi0_Chi2_Good.png";
    else          Hname += "eta_Chi2_Good.png";
    myc1->SaveAs( Hname.Data() );
    tree->Draw("Chisqu_","Chisqu_>5");
    myc1->Update();
    Hname = folder + "/" + nameEBEE + "_";
    if(Are_pi0_)  Hname += "pi0_Chi2_Bad.png";
    else          Hname += "eta_Chi2_Bad.png";
    myc1->SaveAs( Hname.Data() );
    //Long Check
    for(int i=0; i<eventN; i++){
	//Get Variables
	tree->GetEntry(i);
	float fit_b0_, fit_b1_, fit_b2_, fit_b3_, fit_mean_, fit_sigma_, fit_Snorm_, fit_Bnorm_, Chisqu_, Ndof_, ietaX_, iphiY_, iZ_, Index_, Coeff_;
	if(EBEE==0){
	  fit_b0_=EBfit_b0_;fit_b1_=EBfit_b1_;fit_b2_=EBfit_b2_;fit_b3_=EBfit_b3_;fit_mean_=EBfit_mean_;fit_sigma_=EBfit_sigma_;fit_Snorm_=EBfit_Snorm_;fit_Bnorm_=EBfit_Bnorm_;Chisqu_=EBChisqu_;Ndof_=EBNdof_;
	  ietaX_=EBieta_; iphiY_=EBiphi_; iZ_=0, Index_=EBhashedIndex_; Coeff_=EBcoeff_;
	}
	if(EBEE==1){
	  fit_b0_=EEfit_b0_;fit_b1_=EEfit_b1_;fit_b2_=EEfit_b2_;fit_b3_=EEfit_b3_;fit_mean_=EEfit_mean_;fit_sigma_=EEfit_sigma_;fit_Snorm_=EEfit_Snorm_;fit_Bnorm_=EEfit_Bnorm_;Chisqu_=EEChisqu_;Ndof_=EENdof_;
	  ietaX_=EEix_; iphiY_=EEiy_; iZ_=EEzside_; Index_=EEhashedIndex_; Coeff_=EEcoeff_;
	}
	if(Chisqu_==0 && EBEE==0){ IndexBadEB.push_back(Index_); Coo1BadEB.push_back(ietaX_); Coo2BadEB.push_back(iphiY_); Coo3BadEB.push_back(iZ_); }
	if(Chisqu_==0 && EBEE==1){ IndexBadEE.push_back(Index_); Coo1BadEE.push_back(ietaX_); Coo2BadEE.push_back(iphiY_); Coo3BadEE.push_back(iZ_); }
	if( EBEE==1 && Coeff_<0.64 ){ IndexBadRangeEE.push_back(Index_); Coo1BadRangeEE.push_back(ietaX_); Coo2BadRangeEE.push_back(iphiY_); Coo3BadRangeEE.push_back(iZ_);}
	//if(Chisqu_<5) continue;
	//Background
	float xlo = Are_pi0_? 0.08:0.4, xhi = Are_pi0_? 0.22:0.65;
	RooRealVar x("x","#gamma#gamma invariant mass",xlo, xhi, "GeV/c^{2}");
	RooRealVar cb0("cb0","cb0", fit_b0_, fit_b0_, fit_b0_); cb0.setVal(fit_b0_); cb0.setRange(fit_b0_,fit_b0_);
	RooRealVar cb1("cb1","cb1", fit_b1_, fit_b1_, fit_b1_); cb1.setVal(fit_b1_); cb1.setRange(fit_b1_,fit_b1_);
	RooRealVar cb2("cb2","cb2", fit_b2_, fit_b2_, fit_b2_); cb2.setVal(fit_b2_); cb2.setRange(fit_b2_,fit_b2_);
	RooRealVar cb3("cb3","cb3", fit_b3_, fit_b3_, fit_b3_); cb3.setVal(fit_b3_); cb3.setRange(fit_b3_,fit_b3_);
	RooArgList cbpars(cb0,cb1,cb2);
	if(Are_pi0_) cbpars.add( cb3);
	RooChebychev bkg("bkg","bkg model", x, cbpars );
	RooRealVar Nbkg("Nbkg","background yield",fit_Bnorm_,fit_Bnorm_,fit_Bnorm_);
	//Signal
	RooRealVar mean("mean","#pi^{0} peak position", fit_mean_, fit_mean_, fit_mean_, "GeV/c^{2}"); mean.setVal(fit_mean_); mean.setRange(fit_mean_,fit_mean_);
	RooRealVar sigma("sigma","#pi^{0} core #sigma", fit_sigma_, fit_sigma_, fit_sigma_, "GeV/c^{2}"); sigma.setVal(fit_sigma_); sigma.setRange(fit_sigma_,fit_sigma_);
	RooGaussian gaus("gaus","Core Gaussian",x, mean,sigma);
	RooRealVar Nsig("Nsig","#pi^{0} yield",fit_Snorm_,fit_Snorm_,fit_Snorm_);
	RooAddPdf model("model","sig+bkg",RooArgList(gaus,bkg),RooArgList(Nsig,Nbkg));
	//Plot
	RooPlot* xframe = x.frame();
	model.plotOn(xframe);
	model.plotOn(xframe,Components(RooArgSet(bkg)),LineStyle(kDashed), LineColor(kRed));
	xframe->Draw();
	std::stringstream ind_ietaX; ind_ietaX << (int) ietaX_;
	std::stringstream ind_iphiY; ind_iphiY << (int) iphiY_;
	std::stringstream ind_iZ; ind_iZ << (int) iZ_;
	TString nameH = folder + "/" + nameEBEE + "_";
	if( Chisqu_ > 5.) nameH += "HIGH_CHI2_";
	if(Are_pi0_)  nameH += "pi0_" + TString(ind_iZ.str()) + "_" + TString(ind_ietaX.str()) + "_" + TString(ind_iphiY.str()) + ".png";
	else          nameH += "eta_" + TString(ind_iZ.str()) + "_" + TString(ind_ietaX.str()) + "_" + TString(ind_iphiY.str()) + ".png";
	myc1->SaveAs( nameH.Data() );
    }// for event
  }//For EB and EE
  //Fill and write histos
  for( unsigned int j=0; j<Coo1BadEB.size(); j++){
     int BinEta = Coo1BadEB[j] + 86;
     int BinPhi = Coo2BadEB[j] + 1;
     NoFit_EB->SetBinContent(BinEta, BinPhi, 1);
     myfileEB<<"Not Fitted channels: Index "<<IndexBadEB[j]<<" coo1 "<<Coo1BadEB[j]<<" coo2 "<<Coo2BadEB[j]<<" coo3 "<<Coo3BadEB[j]<<"\n";
  }
  for( unsigned int j=0; j<Coo1BadEE.size(); j++){
    if(Coo3BadEE[j]<0) NoFit_EEm->SetBinContent(Coo1BadEE[j]+1, Coo2BadEE[j]+1, 1);
    if(Coo3BadEE[j]>0) NoFit_EEp->SetBinContent(Coo1BadEE[j]+1, Coo2BadEE[j]+1, 1);
    myfileEE<<"Not Fitted channels: Index "<<IndexBadEE[j]<<" coo1 "<<Coo1BadEE[j]<<" coo2 "<<Coo2BadEE[j]<<" coo3 "<<Coo3BadEE[j]<<"\n";
  }
  for( unsigned int j=0; j<Coo1BadRangeEE.size(); j++){
    if( Coo3BadRangeEE[j]<0 ) BadRangeEEm->SetBinContent(Coo1BadRangeEE[j]+1, Coo2BadRangeEE[j]+1, 1);
    if( Coo3BadRangeEE[j]>0 ) BadRangeEEp->SetBinContent(Coo1BadRangeEE[j]+1, Coo2BadRangeEE[j]+1, 1);
    myfileEE_range<<"Bad Range channels: Index "<<IndexBadRangeEE[j]<<" coo1 "<<Coo1BadRangeEE[j]<<" coo2 "<<Coo2BadRangeEE[j]<<" coo3 "<<Coo3BadRangeEE[j]<<"\n";
  }

  gStyle->SetOptStat(0); 
  TString nameH = folder + "/ZeroFitMapEB.png";
  NoFit_EB->Draw("colz");
  myc1->SaveAs( nameH.Data() );
  nameH = folder + "/ZeroFitMapEEm.png";
  NoFit_EEm->Draw("colz");
  myc1->SaveAs( nameH.Data() );
  nameH = folder + "/ZeroFitMapEEp.png";
  NoFit_EEp->Draw("colz");
  myc1->SaveAs( nameH.Data() );
  nameH = folder + "/BadRangeEEm.png";
  BadRangeEEp->Draw("colz");
  myc1->SaveAs( nameH.Data() );
  nameH = folder + "/BadRangeEEp.png";
  BadRangeEEm->Draw("colz");
  myc1->SaveAs( nameH.Data() );

  myfileEB.close();
  myfileEE.close();
  myfileEE_range.close();
}
