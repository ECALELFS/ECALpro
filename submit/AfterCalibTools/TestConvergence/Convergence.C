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

//5_3_6:  gROOT->ProcessLine(".include /afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms9/include/")
//Usage: .x Convergence.C+("/store/group/dpg_ecal/alca_ecalcalib/lpernie/","ALL_2015B_Multifit_01",13,"2015B_")
void Convergence( string Path_0, string Path, int nIter, string Tag, int nJump=1 ){

    string PathL = "root://eoscms//eos/cms" + Path_0 + Path;
    system( (string("mkdir -p plot_") + Path ).c_str());
    TCanvas* myc1 = new TCanvas("myc1", "CMS", 600, 600);
    TString outname = "plot_" + Path + "/Differences.root";
    TFile* output = new TFile(outname.Data(),"RECREATE");
    TH2F* rms_EB  = new TH2F("rms_EB","IC(n)-IC(n-1) #phi on x #eta on y",MAX_IPHI, MIN_IPHI, MAX_IPHI, 2*MAX_IETA+1, -MAX_IETA-0.5, MAX_IETA+0.5 );
    TH2F* rms_EEp = new TH2F("rms_EEp","IC(n)-IC(n-1) iX on x iY on y (EEp)",100,0.5,100.5,100,0.5,100.5);
    TH2F* rms_EEm = new TH2F("rms_EEm","IC(n)-IC(n-1) iY on x iY on y (EEm)",100,0.5,100.5,100,0.5,100.5);

    for(int isEB=0; isEB<2; isEB++){

	  float *EB_RMS = NULL;
	  EB_RMS = new float[nIter];
	  float *iter= NULL;
	  iter = new float[nIter];

	  float hmean(0.), hrms(0.01), sigma_plot(0);
	  if(isEB==1){ hmean=0.09; hrms=0.03; } //rms 0.015

	  for(int i=0; i<nIter; i++){

		//Iter
		stringstream ss; ss<<i;
		stringstream ss1; ss1<<(i+nJump);
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
		//TH1F *h1; h1 =new TH1F("h1","",1000,hmean-9*hrms,hmean+9*hrms);
		TH1F *h1; h1 =new TH1F("h1","",2000,-0.1,0.1);

		//Loop
		Long64_t nentries = Tree->GetEntriesFast();
		for(Long64_t iEntry=0; iEntry<nentries; iEntry++){
		    Tree->GetEntry(iEntry);
		    Tree1->GetEntry(iEntry);
		    if( coeff1!=1. && coeff!=1. && coeff1!=coeff && coeff!=0 && coeff1!=0 /*&& Ndof>10 && Ndof1>10*/){
			  if(isEB==0 ) h1->Fill((coeff1-coeff));
			  if(isEB==1 ) h1->Fill((coeff1-coeff));
		    }

		    if(i==nIter-1){
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
		if(isEB==0) out = "plot_" + Path + "/EB_Iter_" + Iter + ".png";
		if(isEB==1) out = "plot_" + Path + "/EE_Iter_" + Iter + ".png";

		hmean = h1->GetMean();
		hrms  = h1->GetRMS();
		//Fit Method
		RooRealVar x("x","IC distribution",hmean-2.3*hrms, hmean+2.3*hrms,"");
		RooDataHist dh("dh","#gamma#gamma invariant mass",RooArgList(x),h1);
		RooRealVar mean("mean","mean",hmean, hmean-1.5*hrms,hmean+1.5*hrms,"");
		RooRealVar sigma("sigma","#sigma",hrms, hrms-hrms/40.,hrms+hrms/40.,"");
		mean.setRange(hmean-hmean/2.,hmean+hmean/2.);
		sigma.setRange(0.,hrms+hrms/2.);
		RooGaussian gaus("gaus","Core Gaussian",x, mean,sigma);
		RooRealVar Nsig("Nsig","#pi^{0} yield",1000.,0.,1.e7);
		Nsig.setVal( h1->Integral()-h1->Integral()/1000. );
		Nsig.setRange(h1->Integral()/50.,h1->Integral());

		RooRealVar cb0("cb0","cb0", 0.2, -1.,1.);
		RooRealVar cb1("cb1","cb1",-0.1, -1.,1.);
		RooRealVar cb2("cb2","cb2", 0.1, -1.,1.);
		RooArgList cbpars(cb0,cb1,cb2);
		RooChebychev bkg("bkg","bkg model", x, cbpars );
		RooRealVar Nbkg("Nbkg","background yield",h1->Integral()/100.,0.,h1->Integral());
		Nbkg.setVal( h1->Integral()/100. );


		RooAddPdf model1("model","only_gaus",RooArgList(gaus),RooArgList(Nsig));
		RooAddPdf model2("model","sig+bkg",RooArgList(gaus,bkg),RooArgList(Nsig,Nbkg));
		RooAbsPdf* model=0; model = &model2;
		RooNLLVar nll("nll","log likelihood var",*model,dh);//,Extended());
		RooMinuit m(nll);
		m.setVerbose(kFALSE);
		m.migrad();
		RooFitResult* res = m.save();

		RooPlot*  xframe = x.frame(h1->GetNbinsX());
		xframe->SetTitle("");
		dh.plotOn(xframe);
		model->plotOn(xframe);
		//model->plotOn(xframe,Components(bkg),LineStyle(kDashed), LineColor(kRed));
		h1->Draw();
		xframe->Draw("same");
		myc1->SaveAs(out.Data());
		sigma_plot  = sigma.getVal();

		TLatex lat;
		char line[300];
		lat.SetNDC();
		lat.SetTextSize(0.030);
		lat.SetTextColor(1);
		sprintf(line,"SIGMA: %.6f ", sigma_plot );
		float xmin(0.55), yhi(0.80);// ypass(0.05);
		lat.DrawLatex(xmin,yhi, line);
//
//		hmean=mean.getVal();
//		EB_RMS[i]=sigma_plot;
		hmean=h1->GetMean();
		EB_RMS[i]=h1->GetRMS();
		iter[i]=i+1;
	  }

	  TGraph *Conv = new TGraph(nIter, iter, EB_RMS);
	  Conv->SetLineColor(2);
	  Conv->SetLineWidth(1);
	  Conv->SetMarkerColor(2);
	  Conv->SetMarkerStyle(20);
	  Conv->SetMarkerSize(0.5);
	  if(isEB==0) Conv->SetTitle("EB) IC Convergence");
	  if(isEB==1) Conv->SetTitle("EE) IC Convergence");
	  Conv->GetXaxis()->SetTitle("Iter");
	  //Conv->GetYaxis()->SetOffset(1.);
	  if(nJump==1) Conv->GetYaxis()->SetTitle("RMS(ICn+1 - IC)");
	  if(nJump==2) Conv->GetYaxis()->SetTitle("RMS(ICn+2 - IC)");
	  Conv->Draw("ACP");
	  myc1->cd();
	  TString out;
	  if(isEB==0) out = "plot_" + Path + "/EB_IC_Convergence.png";
	  if(isEB==1) out = "plot_" + Path + "/EE_IC_Convergence.png";
	  myc1->SaveAs(out.Data());
    }
    output->cd();
    rms_EB->Write();
    rms_EEp->Write();
    rms_EEm->Write();
    output->Close();

//    delete output;
//    delete rms_EB;
//    delete rms_EEp;
//    delete rms_EEm;
}
