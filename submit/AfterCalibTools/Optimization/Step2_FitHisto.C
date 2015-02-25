#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TF2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TFormula.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMarker.h>
#include <TMarker.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TFitResultPtr.h>
#include <TVirtualFitter.h>
#include <TMinuit.h>
#include <TChain.h>
#include <TObject.h>
#include <iostream>
#include <memory>
#include <string>

#include "TTree.h"
#include "TLatex.h"
#include "TMath.h"
#include "TBranch.h"
#include "TFile.h"
#include "TStyle.h"
#include "TPad.h"
#include "TProfile.h"
#include "TPave.h"
#include "TPaveText.h"
#include "TString.h"
#include "TPaveStats.h"

using namespace std;
using namespace RooFit;

struct FFitResult {

    float S;     
    float B;  
    float SoB;      

    float chi2;
    int   dof;
    double neve; 
    double errore;

    double taglio;
    double divisione;
    RooPlot*  fram;
};

FFitResult FitFit(TH1D* h, double xmin, double xmax, int events, int Name, FILE *file_txt, float eff , bool Are_pi0_);

//.x Step2_FitHisto.C("ALL_MINBIAS_UNCAL_L1_NOL1FILTER", "Fstep1_EB_eta.root", true, false)
void Step2_FitHisto(TString folder="ALL_MINBIAS_UNCAL_L1_NOL1FILTER", TString file="Fstep1_EE_eta.root", bool isEB=false, bool Are_pi0_=false){

    //Open Input File
    TFile *f1 = new TFile( (folder + "/" + file).Data(),"r" );
    TH2F *histo=(TH2F*)f1->Get("hmass");
    TH2F *histo_tot=(TH2F*)f1->Get("hmass_tot");
    //Create Output files
    TString out_f;
    if(isEB){
	  out_f="Fstep2_EB_pi0.root";
	  if(!Are_pi0_) out_f="Fstep2_EB_eta.root";
    }
    else{
	  out_f="Fstep2_EE_pi0.root";
	  if(!Are_pi0_) out_f="Fstep2_EE_eta.root";
    }
    TFile* outPut = new TFile( (folder + "/" + out_f).Data(),"RECREATE" );
    outPut->cd();
    FILE *file_txt;
    if(isEB){
	  out_f="Fstep2_notSorted_EB_pi0.txt";
	  if(!Are_pi0_) out_f="Fstep2_notSorted_EB_eta.txt";
    }
    else    {
	  out_f="Fstep2_notSorted_EE_pi0.txt";
	  if(!Are_pi0_) out_f="Fstep2_notSorted_EE_eta.txt";
    }
    file_txt=fopen( (folder + "/" + out_f).Data(),"w" );
    //FIT
    cout << "Now starting the Fit procedure ..." << endl;
    double                 xmin=0.08, xmax=0.2;   //Pi0 EB
    if(Are_pi0_ && !isEB){ xmin=0.07; xmax=0.21; }//Pi0 EE
    if(!Are_pi0_){         xmin=0.23; xmax=0.6;  }//Eta
    cout<<"You have "<<histo->GetNbinsX()<<" bins..."<<endl;  
    for(int i=0; i<histo->GetNbinsX(); i++){
	  TH1D *h1;
	  TH1D *h1_tot;
	  h1 = histo->ProjectionY("_py",i+1,i+1);
	  h1_tot = histo_tot->ProjectionY("_py",i+1,i+1);
	  int Nentr = h1->GetEntries();
	  int Nentr_tot = h1_tot->GetEntries();
	  int iMin(0), iMax(0);
	  if(Are_pi0_){ 
		iMin = h1->GetXaxis()->FindBin(0.08);
		iMax = h1->GetXaxis()->FindBin(0.18);
	  }
	  else{
		iMin = h1->GetXaxis()->FindBin(0.35);
		iMax = h1->GetXaxis()->FindBin(0.6);
	  }
	  double integral = h1->Integral(iMin, iMax);
	  float eff =  (float)Nentr/(float)Nentr_tot;
	  if(eff>0.002){
		FFitResult res;
		res = FitFit(h1, xmin, xmax, Nentr, i, file_txt, eff, Are_pi0_);
	  }
    }
    outPut->Write();
    outPut->Close();
}

FFitResult FitFit(TH1D* h, double xmin, double xmax, int events, int Name, FILE *file_txt, float eff , bool Are_pi0_){

    stringstream ss; ss << Name; 
    TString NameTrue = "Bin_"+ss.str();
    TCanvas* myc1 = new TCanvas(NameTrue.Data(), "CMS", 600, 600); 
    float effic = ( h->GetEntries()/events)*100;

    RooRealVar x("x","#pi_{0} invariant mass",xmin, xmax, "GeV/c^{2}");
    RooDataHist dh("dh","#gamma#gamma invariant mass",RooArgList(x),h);
    RooRealVar mean("mean","#pi^{0} peak position", Are_pi0_?0.116:0.57,  Are_pi0_?0.105:0.45, Are_pi0_?0.150:0.65,"GeV/c^{2}");
    RooRealVar sigma("sigma","#pi^{0} core #sigma",Are_pi0_?0.010:0.08, 0.005,Are_pi0_?0.025:0.15,"GeV/c^{2}");

    RooRealVar Nsig("Nsig","#pi^{0} yield",1000,0.,1.e7);
    Nsig.setVal( h->GetSum()*0.1);
    RooGaussian gaus("gaus","Core Gaussian",x, mean,sigma);

    RooRealVar cb0("cb0","cb0", 0.2, -1.,1.);
    RooRealVar cb1("cb1","cb1",-0.1, -1.,1.);
    RooRealVar cb2("cb2","cb2", 0.1, -1.,1.);
    RooRealVar cb3("cb3","cb3",-0.1, -1.,1.);
    RooRealVar cb4("cb4","cb4", 0.1, -1.,1.);
    RooRealVar cb5("cb5","cb5", 0.1, -1.,1.);
    RooRealVar cb6("cb6","cb6", 0.3, -1.,1.);

    RooArgList cbpars(cb0,cb1,cb2,cb3);
    RooChebychev bkg("bkg","bkg model", x, cbpars );

    RooRealVar Nbkg("Nbkg","background yield",1.e3,0.,h->GetSum());
    Nbkg.setVal( h->GetSum()*0.9);

    RooAbsPdf* model=0;

    RooAddPdf model1("model","sig+bkg",RooArgList(gaus,bkg),RooArgList(Nsig,Nbkg));

    model = &model1;

    RooNLLVar nll("nll","log likelihood var",*model,dh,Extended());
    RooMinuit m(nll);
    m.migrad();

    RooFitResult* res = m.save() ;

    RooChi2Var chi2("chi2","chi2 var",*model,dh, Extended());
    int ndof = h->GetNbinsX() - res->floatParsFinal().getSize();

    x.setRange("sobRange",mean.getVal()-2.*sigma.getVal(), mean.getVal()+2.*sigma.getVal());
    RooAbsReal* integralSig = gaus.createIntegral(x,NormSet(x),Range("sobRange"));
    RooAbsReal* integralBkg = bkg.createIntegral(x,NormSet(x),Range("sobRange"));

    float normSig = integralSig->getVal();
    float normBkg = integralBkg->getVal();

    RooPlot*  xframe = x.frame();
    xframe->SetTitle("frame");
    dh.plotOn(xframe);
    model->plotOn(xframe,Components(bkg),LineStyle(kDashed), LineColor(3));
    model->plotOn(xframe,Components(gaus),LineColor(2));

    model->plotOn(xframe);

    xframe->GetXaxis()->SetTitleOffset(1.5);

    xframe->Draw();

    FFitResult result;
    result.S=normSig*Nsig.getVal();
    result.B=normBkg*Nbkg.getVal();
    result.SoB=(result.S/result.B);     
    result.chi2=xframe->chiSquare();                                                                                                                    
    result.dof=ndof;                                                                                                                               
    result.neve=effic;                                                                                         
    result.errore = result.SoB*sqrt( pow(normSig*Nsig.getError()/result.S,2) +pow(normBkg*Nbkg.getError()/result.B,2) ) ;

    result.fram=xframe;  

    TLatex lat;
    char line[300];
    lat.SetNDC();
    lat.SetTextSize(0.03);
    lat.SetTextColor(1);
    float Xmin(0.55), Yhi(0.70), Ypass(0.05);
    sprintf(line,"S/B: %.5f", result.S/result.B);
    lat.DrawLatex(Xmin,Yhi, line);
    sprintf(line,"S(B): %.5f(%.5f) ", result.S, result.B );
    lat.DrawLatex(Xmin,Yhi-Ypass, line);
    sprintf(line,"Chi^2: %.5f", result.chi2/result.dof  );
    lat.DrawLatex(Xmin,Yhi-2*Ypass, line);
    sprintf(line,"#mu: %.5f , #sigma(#mu): %.5f", mean.getVal(), sigma.getVal() );
    lat.DrawLatex(Xmin,Yhi-3.*Ypass, line);
    sprintf(line,"#sigma(#mu)/#mu: %.5f ", mean.getError()/mean.getVal() );
    lat.DrawLatex(Xmin,Yhi-4.*Ypass, line);
    sprintf(line,"Efficiency: %.5f", eff );
    lat.DrawLatex(Xmin,Yhi-5.*Ypass, line);
    if( xframe->chiSquare()/result.dof<0.02 && result.S/result.B>0.01 && mean.getVal()<Are_pi0_?0.2:0.6 && mean.getVal()>Are_pi0_?0.01:0.45 ) 
	  fprintf(file_txt,"BIN %i  SB %.5f  MuSi %.5f  CHI %.5f  Eff %.5f \n", Name, result.S/result.B, mean.getError()/mean.getVal(), result.chi2/result.dof, eff );
    myc1->Write();
    delete myc1;
}
