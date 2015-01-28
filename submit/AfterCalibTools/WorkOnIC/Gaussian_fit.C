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
//7_1_0:  gROOT->ProcessLine(".include /afs/cern.ch/cms/slc6_amd64_gcc481/lcg/roofit/5.34.22-cms2/include/")
//Usage: .x Gaussian_fit.C+("Compare2012DMINECCnot","2012D_EtaRingMY_Vs_MYnoCC.root")
void Gaussian_fit( string Dir, string File ){

    string PathL = Dir + "/" + File;
    TCanvas* myc1 = new TCanvas("myc1", "CMS", 600, 600);
    TString outname = Dir + "/Gaus_Fit.root";

    //Files
    cout<<"Opening: "<<PathL<<endl;
    TFile* fin = TFile::Open(PathL.c_str());
    //Histos
    string List[4]={"Ratio1D_EB","Ratio1D_EEm","Ratio1D_EEp","Ratio1D_EB_eta1"};

    for(int nH=0; nH<4; nH++){

	  TH1F *h_toFit = (TH1F*) fin->Get( List[nH].c_str() );
	  float hmean = h_toFit->GetMean();
	  float hrms = h_toFit->GetRMS();
	  TString outName = Dir + "/" + "Ratio1D_EB.png";
	  if(nH==1) outName = Dir + "/" + "Ratio1D_EEm.png";
	  if(nH==2) outName = Dir + "/" + "Ratio1D_EEp.png";
	  if(nH==3) outName = Dir + "/" + "Ratio1D_EB_eta1.png";

	  //Fit Method
	  RooRealVar x("x","IC distribution",hmean-1.6*hrms, hmean+1.6*hrms, "");
	  RooDataHist dh("dh","#gamma#gamma invariant mass",RooArgList(x),h_toFit);
	  RooRealVar mean("mean","mean",hmean, hmean-1.5*hrms,hmean+1.5*hrms,"");
	  RooRealVar sigma("sigma","#sigma",hrms, hrms-hrms/40.,hrms+hrms/40.,"");
	  mean.setRange(hmean-hmean/2.,hmean+hmean/2.);
	  sigma.setRange(0.,hrms+hrms/2.);
	  RooGaussian gaus("gaus","Core Gaussian",x, mean,sigma);
	  RooRealVar Nsig("Nsig","#pi^{0} yield",1000.,0.,1.e7);
	  Nsig.setVal( h_toFit->Integral()-h_toFit->Integral()/1000. );
	  Nsig.setRange(h_toFit->Integral()/50.,h_toFit->Integral());

	  RooRealVar cb0("cb0","cb0", 0.2, -1.,1.);
	  RooRealVar cb1("cb1","cb1",-0.1, -1.,1.);
	  RooRealVar cb2("cb2","cb2", 0.1, -1.,1.);
	  RooArgList cbpars(cb0,cb1,cb2);
	  RooChebychev bkg("bkg","bkg model", x, cbpars );
	  RooRealVar Nbkg("Nbkg","background yield",h_toFit->Integral()/100.,0.,h_toFit->Integral());
	  Nbkg.setVal( h_toFit->Integral()/100. );


	  RooAddPdf model1("model","only_gaus",RooArgList(gaus),RooArgList(Nsig));
	  RooAddPdf model2("model","sig+bkg",RooArgList(gaus,bkg),RooArgList(Nsig,Nbkg));
	  RooAbsPdf* model=0; model = &model2;
	  RooNLLVar nll("nll","log likelihood var",*model,dh);//,Extended());
	  RooMinuit m(nll);
	  m.setVerbose(kFALSE);
	  m.migrad();
	  //RooFitResult* res = m.save();

	  RooPlot*  xframe = x.frame(h_toFit->GetNbinsX());
	  xframe->SetTitle("");
	  dh.plotOn(xframe);
	  model->plotOn(xframe);
	  //model->plotOn(xframe,Components(bkg),LineStyle(kDashed), LineColor(kRed));
	  h_toFit->Draw();
	  xframe->Draw("same");
	  //myc1->SaveAs(outName.Data());
	  float sigma_plot  = sigma.getVal();

	  TLatex lat;
	  char line[300];
	  lat.SetNDC();
	  lat.SetTextSize(0.030);
	  lat.SetTextColor(1);
	  sprintf(line,"SIGMA: %.6f ", sigma_plot );
	  float xmin(0.55), yhi(0.80);// ypass(0.05);
	  lat.DrawLatex(xmin,yhi, line);
	  myc1->SaveAs(outName.Data());
    }
    //rms_EB->Write();
    //rms_EEp->Write();
    //rms_EEm->Write();
    //output->Close();
}
