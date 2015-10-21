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

FFitResult FitFit(TH1D* h, double xmin, double xmax, int events, int Name, FILE *file_txt, float eff , bool Are_pi0_, int LH, bool isEB);

//.x Step2_FitHisto.C("ALL_MINBIAS_UNCAL_L1_NOL1FILTER", "Fstep1_EB_eta.root", true, false)
void Step2_FitHisto(TString folder="ALL_MINBIAS_UNCAL_L1_NOL1FILTER_20bx25_7e33_noCC", TString file="Fstep1_EE_pi0.root", bool isEB=false, bool Are_pi0_=true){

  //Open Input File
  TFile *f1 = new TFile( (folder + "/" + file).Data(),"r" );
  TH2F *histoL, *histoH, *histo_totL, *histo_totH;
  //Create Output files
  TString out_f;
  if(isEB){
    out_f="Fstep2_EB_pi0.root";
    histoL = (TH2F*)f1->Get("hmassEBL");
    histoH = (TH2F*)f1->Get("hmassEBH");
    histo_totL = (TH2F*)f1->Get("hmass_totEBL");
    histo_totH = (TH2F*)f1->Get("hmass_totEBH");
    if(!Are_pi0_) out_f="Fstep2_EB_eta.root";
  }
  else{
    histoL = (TH2F*)f1->Get("hmassEEL");
    histoH = (TH2F*)f1->Get("hmassEEH");
    histo_totL = (TH2F*)f1->Get("hmass_totEEL");
    histo_totH = (TH2F*)f1->Get("hmass_totEEH");
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
  double                  xmin=0.1, xmax=0.18;   //Pi0 EB L
  //Low and High eta region
  TH2F *histo, *histo_tot;
  for(int LH=0; LH<2; LH++){
    if(Are_pi0_ && isEB && LH==1 ){ xmin=0.1; xmax=0.19; } //Pi0 EB H
    if(Are_pi0_ && !isEB){          xmin=0.07; xmax=0.2; } //Pi0 EE
    if(!Are_pi0_ && isEB && LH==0){ xmin=0.39; xmax=0.68; }//Eta EB L
    if(!Are_pi0_ && isEB && LH==1){ xmin=0.40; xmax=0.65; } //Eta EB H
    if(!Are_pi0_ && !isEB){         xmin=0.32; xmax=0.9; } //Eta EE
    if( LH==0 ){
	histo = histoL;
	histo_tot = histo_totL;
    }
    if( LH==1 ){
	histo = histoH;
	histo_tot = histo_totH;
    }
    cout<<"You have "<<histo->GetNbinsX()<<" bins..."<<endl;  
    //Loop on the histo
    for(int i=0; i<histo->GetNbinsX(); i++){
//if( LH==0 ) continue;
//if( i!=5604 && i!=3392 && i!=6256 && i!=304 && i!=4160 && i!=2512 ) continue;
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
	  iMin = h1->GetXaxis()->FindBin(0.47);
	  iMax = h1->GetXaxis()->FindBin(0.57);
	}
	double integral = h1->Integral(iMin, iMax);
	float eff =  (float)Nentr/(float)Nentr_tot;
	if(eff>0.002){
	  vector<float> chi2s_sb; chi2s_sb.clear();
	  int Try=1;
	  chi2s_sb = FitFit(h1, xmin, xmax, Nentr, i, file_txt, eff, Are_pi0_, LH, isEB);
	  while( (chi2s_sb[0]>0.05 || chi2s_sb[2]>2.) && Try<10 ){
	    chi2s_sb = FitFit(h1, xmin+Try*0.0005, xmax-Try*0.001, Nentr, i, file_txt, eff, Are_pi0_, LH, isEB);
	    Try++;
	  }
	}
    }
  }
  outPut->Write();
  outPut->Close();
  delete outPut;
  delete f1;
}

vector<float> FitFit(TH1D* h, double xmin, double xmax, int events, int Name, FILE *file_txt, float eff , bool Are_pi0_, int LH, bool isEB){

  stringstream ss; ss << Name;
  TString preName = "BinL_";
  if( LH==1 ) preName = "BinH_";
  TString NameTrue = preName.Data() + ss.str();
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
  RooRealVar cb1("cb1","cb1", 0., (Are_pi0_ && LH==1) ? -1.:-0.5, (Are_pi0_ && LH==1) ? 1.:0.5);
  RooRealVar cb2("cb2","cb2", 0., (Are_pi0_ && LH==1) ? -1.:-0.5, (Are_pi0_ && LH==1) ? 1.:0.5);

  RooArgList cbpars(cb0,cb1);
  if( Are_pi0_  && !isEB )         cbpars.add(cb2);
  if( Are_pi0_  && isEB && LH==0 ) cbpars.add(cb2);
  if( !Are_pi0_ &&         LH==0 ) cbpars.add(cb2);
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
  delete res;
  x.setRange("sobRange",mean.getVal()-2.*sigma.getVal(), mean.getVal()+2.*sigma.getVal());
  RooAbsReal* integralSig = gaus.createIntegral(x,NormSet(x),Range("sobRange"));
  RooAbsReal* integralBkg = bkg.createIntegral(x,NormSet(x),Range("sobRange"));

  float normSig = integralSig->getVal();
  float normBkg = integralBkg->getVal();
  delete integralSig;
  delete integralBkg;
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
  sprintf(line,"S/B: %.3f", result.S/result.B);
  lat.DrawLatex(Xmin,Yhi, line);
  sprintf(line,"S(B): %.3f(%.3f) ", result.S, result.B );
  lat.DrawLatex(Xmin,Yhi-Ypass, line);
  sprintf(line,"Chi^2: %.3f", result.chi2/result.dof  );
  lat.DrawLatex(Xmin,Yhi-2*Ypass, line);
  sprintf(line,"#mu: %.3f , #sigma(#mu): %.3f", mean.getVal(), sigma.getVal() );
  lat.DrawLatex(Xmin,Yhi-3.*Ypass, line);
  sprintf(line,"#sigma(#mu)/#mu: %.3f ", mean.getError()/mean.getVal() );
  lat.DrawLatex(Xmin,Yhi-4.*Ypass, line);
  sprintf(line,"Efficiency: %.3f", eff );
  lat.DrawLatex(Xmin,Yhi-5.*Ypass, line);
  if( xframe->chiSquare()/result.dof<0.05 && result.S/result.B<2. ){
    if(LH==0) fprintf(file_txt,"L_BIN %i  SB %.5f  MuSi %.5f  CHI %.5f  Eff %.5f \n", Name, result.S/result.B, mean.getError()/mean.getVal(), result.chi2/result.dof, eff );
    if(LH==1) fprintf(file_txt,"H_BIN %i  SB %.5f  MuSi %.5f  CHI %.5f  Eff %.5f \n", Name, result.S/result.B, mean.getError()/mean.getVal(), result.chi2/result.dof, eff );
    myc1->Write();
  }
  delete myc1;
  delete xframe;
  vector<float> chi2s_sb; //chi2s_sb.push_back(result.chi2/result.dof); chi2s_sb.push_back(result.S/result.B);
  chi2s_sb.push_back(0.); chi2s_sb.push_back(0.);
  return chi2s_sb;
}
