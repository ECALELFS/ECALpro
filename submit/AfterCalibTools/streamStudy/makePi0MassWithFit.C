#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2F.h>
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
#include <iomanip>
#include <cmath>
#include <sstream> //to use ostringstream to convert numbers to string in c++    

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMinimizer.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "./CMS_lumi.h"
//#include "/afs/cern.ch/work/m/mciprian/w_mass_analysis/CMSSW_8_0_25/src/WmassAnalysis/macros/utility.h"
#include "./utility.h"

using namespace std;
using namespace RooFit;

int Xtal_ID[170][360]={0};
int Xtal_Ieta[61200]={0};
int Xtal_Iphi[61200]={0};

int Xtal_Ix[14648]={0};
int Xtal_Iy[14648]={0};
int Xtal_Iz[14648]={0};

static string endcap_ix_iy_zside_ietaRing = "/afs/cern.ch/user/m/mciprian/public/ECALproTools/EE_xyzToEtaRing/eerings_modified.root";
static string deadXtalFileName = "/afs/cern.ch/user/m/mciprian/public/ECALproTools/test_DeadXtal_AlCaP0_Run2017B_3July_upToRun297723/h_DeadXtal.root";
static bool drawAllMassPlot = true;

//=================================================

class vectorManager {

public:
  vectorManager() { };

  vectorManager(const vector<TH1*> & histPtrs,
		const vector<string> & histNames,
		const vector<string> & histLegends
		)
  {
    histPtrs_    = vector<TH1*  >(histPtrs);
    histNames_   = vector<string>(histNames);
    histLegends_ = vector<string>(histLegends);
  };

  ~vectorManager() {};

  vector<TH1*>   getHistPtrs()    const { return histPtrs_;    };
  vector<string> getHistNames()   const { return histNames_;   };
  vector<string> getHistLegends() const { return histLegends_; };

  void addComponent(TH1* histPtr = NULL, const string& histName = "name", const string& histLegend = "leg")  { 
    histPtrs_.push_back(histPtr);
    histNames_.push_back(histName);
    histLegends_.push_back(histLegend);
  };

private:

  vector<TH1*> histPtrs_;
  vector<string> histNames_;
  vector<string> histLegends_;

};


//=================================================

bool noDeadXtalIn3x3matrixSeededByThisXtal(const TH2F* hDeadXtals = NULL, const int x = 1, const int y = 1) {

  // WARNING: it is assumed that the seed is already not adjacent to a gap, therefore we won't have abs(eta)=0 or abs(ieta)=85 or iphi=1 or iphi=360 for the seed 

  int nDeadXtals = 0;

  for (int xspan = x-1; xspan <= x+1 && nDeadXtals == 0; xspan++) {
    for (int yspan = y-1; yspan <= y+1 && nDeadXtals == 0; yspan++) {
      nDeadXtals += (int) (0.5 + hDeadXtals->GetBinContent(xspan,yspan)); // histogram returns float, to avoid bad truncation sum 0.5 and then round to integer 
    }
  }

  return (nDeadXtals == 0) ? true : false;

}

//=================================================

// void drawHisto(TH1* hSum = NULL, 
// 	       const bool isEB = true, 
// 	       const string& outDir = "./", 
// 	       const string& hName = "", 
// 	       const double lumi = 8.6) {

//   TGaxis::SetMaxDigits(3); 

//   TCanvas* canvas = new TCanvas("canvas","",600,600);
//   canvas->cd();
//   canvas->SetTickx(1);
//   canvas->SetTicky(1);
//   canvas->cd();
//   canvas->SetRightMargin(0.06);

//   hSum->SetStats(0);
//   hSum->SetLineColor(kBlack);
//   hSum->SetMarkerColor(kBlack);
//   hSum->SetMarkerStyle(20);
//   hSum->SetMarkerSize(1);

//   hSum->SetTitle(0);
  
//   hSum->GetXaxis()->SetLabelSize(0.04);
//   hSum->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");
//   hSum->GetXaxis()->SetTitleSize(0.05);  
//   hSum->GetXaxis()->SetTitleOffset(0.9);
//   hSum->GetXaxis()->SetRangeUser(0.05,0.25);

//   double maxY = hSum->GetBinContent(hSum->GetMaximumBin());
//   hSum->GetYaxis()->SetRangeUser(0.0, 1.2*maxY);
//   hSum->GetYaxis()->SetTitle("#gamma#gamma pairs / 0.004 GeV/c^{2}");
//   hSum->GetYaxis()->SetTitleOffset(1.1);
//   hSum->GetYaxis()->SetTitleSize(0.05);
//   hSum->Draw("EP");

//   /////////////////////////
//   /////////////////////////
//   // fit and draw result

//   // better not to use extreme values of histogram as the fit range, because the actual range is shorter
//   // see example here --> http://mciprian.web.cern.ch/mciprian/test_plot/pi0Mass_EB_h_xtal_iter0.png
//   // I suggest using 0.080 and 0.21 for pi0
//   RooRealVar x("x","#gamma#gamma invariant mass", Are_pi0_? 0.07:0.4, Are_pi0_? 0.21:0.65, "GeV/c^2");
//   if (Are_pi0_ && not isEB) x.setRange(0.075, 0.24);

//   RooDataHist dh("dh","#gamma#gamma invariant mass",RooArgList(x),hSum);

//   RooRealVar mean("mean","#pi^{0} peak position", Are_pi0_? 0.13:0.52,  Are_pi0_? 0.105:0.5, Are_pi0_? upper_bound_pi0mass_EB:upper_bound_etamass_EB,"GeV/c^{2}");
//   RooRealVar sigma("sigma","#pi^{0} core #sigma",0.011, 0.005,0.015,"GeV/c^{2}");
//   if(not isEB)  {
//     mean.setRange( Are_pi0_? 0.1:0.45, Are_pi0_? upper_bound_pi0mass_EE:upper_bound_etamass_EE);
//     mean.setVal(Are_pi0_? 0.13:0.55);
//     sigma.setRange(0.005, 0.020);
//   }

//   RooRealVar Nsig("Nsig","#pi^{0} yield", hSum->Integral()*0.15,0.,hSum->Integral()*10.0);
//   //Nsig.setVal( hSum->Integral()*0.1);

//   //sig model
//   RooGaussian gaus("gaus","Core Gaussian",x, mean,sigma);

//   // bkg model
//   RooRealVar cb0("cb0","cb0", 0.2, -1.,1.);
//   RooRealVar cb1("cb1","cb1",-0.1, -1.,1.);
//   RooRealVar cb2("cb2","cb2", 0.1,  -1.,1.);
//   RooRealVar cb3("cb3","cb3", -0.1, -0.5,0.5);
//   RooArgList cbpars(cb0,cb1,cb2,cb3);
//   RooChebychev bkg("bkg","bkg model", x, cbpars );

//   RooRealVar Nbkg("Nbkg","background yield",hSum->Integral()*0.85,0.,hSum->Integral()*10.0);
//   //Nbkg.setVal( hSum->Integral()*0.8 );

//   RooAbsPdf* model = 0;
//   // can use many models
//   RooAddPdf model1("model","sig+bkg",RooArgList(gaus,bkg),RooArgList(Nsig,Nbkg));
//   // modelXXX ...
//   model = &model1;

//   RooNLLVar nll("nll","log likelihood var",*model,dh, RooFit::Extended(true));
//   //RooAbsReal * nll = model->createNLL(dh); //suggetsed way, taht should be the same                                                                                      

//   // FIT 1
//   // copied from ECALpro
//   // RooMinimizer m(nll);
//   // m.setVerbose(kFALSE);
//   // //m.setVerbose(kTRUE);                                                                                                       
//   // m.migrad();
//   // //m.hesse();                                                                                                                          
//   // RooFitResult* res = m.save() ;

//   // FIT2
//   // copied from Raffaele Gerosa
//   RooMinimizer mfit(nll);
//   mfit.setVerbose(kFALSE);
//   mfit.setPrintLevel(-1);
//   cout << "######### Minimize" << endl;
//   mfit.minimize("Minuit2","minimize");
//   cout << "######### Minimize hesse " << endl;
//   mfit.minimize("Minuit2","hesse");
//   cout<<"######### Estimate minos errors for all parameters"<<endl;
//   mfit.minos(RooArgSet(Nsig,Nbkg));
//   RooFitResult* res = mfit.save("res") ;

//   // FIT 1 and FIT 2 yields practically the same result, using the second

//   cout << "print fit result" << endl;
//   res->Print();

//   RooChi2Var chi2("chi2","chi2 var",*model,dh, true);

//   int ndof = hSum->GetNbinsX() - res->floatParsFinal().getSize();

//   //compute S/B and chi2                 
//   // use 3 sigma range around mean to get S/B                                                                                       
//   x.setRange("sobRange",mean.getVal()-3.*sigma.getVal(), mean.getVal()+3.*sigma.getVal());
//   RooAbsReal* integralSig = gaus.createIntegral(x,NormSet(x),Range("sobRange"));

//   RooAbsReal* integralBkg = bkg.createIntegral(x,NormSet(x),Range("sobRange"));

//   float normSig = integralSig->getVal();
//   float normBkg = integralBkg->getVal();

//   Pi0FitResult pi0res; // this is the output value of this method                                                                                                          
//   pi0res.res = res;

//   pi0res.S = normSig*Nsig.getVal();
//   pi0res.Serr = normSig*Nsig.getError();

//   pi0res.B = normBkg*Nbkg.getVal();
//   pi0res.Berr = normBkg*Nbkg.getError();

//   pi0res.SoB =  pi0res.S/pi0res.B;
//   pi0res.SoBerr =  pi0res.SoB*sqrt( pow(pi0res.Serr/pi0res.S,2) +
// 				    pow(pi0res.Berr/pi0res.B,2) ) ;
//   pi0res.dof = ndof;

//   RooPlot*  xframe = x.frame(hSum->GetNbinsX());
//   xframe->SetTitle(0);
//   dh.plotOn(xframe);  // already drawn
//   model->plotOn(xframe,Components(bkg),LineStyle(kDashed), LineColor(kRed));
//   model->plotOn(xframe);

//   xframe->Draw("same");

//   TLatex lat;
//   char line[300];
//   lat.SetNDC();
//   lat.SetTextSize(0.035);
//   lat.SetTextColor(1);

//   float xmin(0.6), yhi(0.85), ypass(0.05);
//   if(isEB and not Are_pi0_) yhi=0.30;
//   sprintf(line,"Nsig: %.0f #pm %.0f", Nsig.getVal(), Nsig.getError() );
//   lat.DrawLatex(xmin,yhi, line);

//   sprintf(line,"m_{#gamma#gamma}: %.2f #pm %.2f", mean.getVal()*1000., mean.getError()*1000. );
//   lat.DrawLatex(xmin,yhi-ypass, line);

//   sprintf(line,"#sigma: %.2f #pm %.2f (%.2f%s)", sigma.getVal()*1000., sigma.getError()*1000., sigma.getVal()*100./mean.getVal(), "%" );
//   lat.DrawLatex(xmin,yhi-2.*ypass, line);

//   sprintf(line,"S/B(3#sigma): %.2f #pm %.2f", pi0res.SoB, pi0res.SoBerr );                                                                                               
//   lat.DrawLatex(xmin,yhi-3.*ypass, line);

//   sprintf(line,"#Chi^{2}: %.2f", xframe->chiSquare()/pi0res.dof );
//   lat.DrawLatex(xmin,yhi-4.*ypass, line);

//   // if using a function to do the fit, can return the output
//   // Pi0FitResult fitres = pi0res;
//   // return fitres;

//   // end of fit part
//   /////////////////////////
//   /////////////////////////

//   canvas->RedrawAxis("sameaxis");
//   if (lumi < 1.0) CMS_lumi(canvas,Form("%.2f",lumi),false,false);
//   else CMS_lumi(canvas,Form("%.1f",lumi),false,false);
//   setTDRStyle();

//   string title = "pi0mass";
//   title += isEB ? "_EB_" : "_EE_";
//   string canvasTitle = outDir + title + hName;
//   canvas->SaveAs((canvasTitle + ".png").c_str());

//   // save fit parameters in file named as the canvas but with txt extension
  
//   string fitParameterFileName = canvasTitle + ".txt";
//   ofstream fitParameterFile(fitParameterFileName.c_str(),ios::out);
//   if ( !fitParameterFile.is_open() ) {  
//     cout<<"Error: unable to open file " << fitParameterFileName <<" !"<<endl;
//     exit(EXIT_FAILURE);
//   } else {
//     pi0res.probchi2 = TMath::Prob(xframe->chiSquare(), ndof);
//     pi0res.chi2 = xframe->chiSquare();
//     fitParameterFile << "FIT PARAMETERS:" << endl;
//     fitParameterFile << "-------------------" << endl;
//     fitParameterFile << setw(15) <<  " Nsig: " << Nsig.getVal() << " +/- " << Nsig.getError() << endl;
//     fitParameterFile << setw(15) <<  " Nbkg: " << Nbkg.getVal() << " +/- " << Nbkg.getError() << endl;
//     fitParameterFile << setw(15) <<  " Nsig(3sigma): " << pi0res.S << " +/- " << pi0res.Serr << endl;
//     fitParameterFile << setw(15) <<  " Nbkg(3sigma): " << pi0res.B << " +/- " << pi0res.Berr << endl;
//     fitParameterFile << setw(15) <<  " S/B(3sigma) : " << pi0res.SoB << " +/- " << pi0res.SoBerr << endl;
//     fitParameterFile << setw(15) <<  " mean : " << mean.getVal()*1000. << " +/- " << mean.getError()*1000. << endl;
//     fitParameterFile << setw(15) <<  " sigma : " << sigma.getVal()*1000. << " +/- " << sigma.getError()*1000. << endl;
//     fitParameterFile << setw(15) <<  " chi2: " << xframe->chiSquare() << endl;
//     fitParameterFile << setw(15) <<  " DOF: " << pi0res.dof << endl;
//     fitParameterFile << setw(15) <<  " prob(chi2): " << pi0res.probchi2 << endl;
//     fitParameterFile.close();
//   }


// }

//=================================================   

void doPi0MassWithFit(TH1* h,
		      const string& hName = "h",
		      const bool isEB = true,
		      const string& plotType = "etaring",
		      const string& outDir = "./",
		      const double lumi = 8.6,
		      const string& filename = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/AlCaP0_IC2017_upTo31July2017_noCC/iter_0/AlCaP0_IC2017_upTo31July2017_noCC_epsilonPlots.root"
		      ) 
{

  // force digits on axis (I wanted it only on y axis to trigger exponential notation, but it looks like it is not implemented for a single axis) 
  TGaxis::SetMaxDigits(3); 

  TH1::SetDefaultSumw2();
  gROOT->SetBatch(kTRUE);

  if (isEB) {

    for(int i = 0; i < 61200; i++)
      {
	int det_ID = EBDetId::detIdFromDenseIndex(i);

	EBDetId ebseed(det_ID);
        int ieta = ebseed.ieta();
        int iphi = ebseed.iphi();		
	Xtal_Ieta[i] = ieta;
	Xtal_Iphi[i] = iphi;
	//	cout<<i<<"   "<<ieta<<"  "<<iphi<<endl;
      }

  } else {

    // EE
    for(int i = 0; i < 14648; i++)
      {
	int det_ID = EEDetId::detIdFromDenseIndex(i);

	// TO BE TESTED
	EEDetId eeseed(det_ID);
        int ix = eeseed.ix();
        int iy = eeseed.iy();		       
	int iz = eeseed.zside();		
	Xtal_Ix[i] = ix;
	Xtal_Iy[i] = iy;
	Xtal_Iz[i] = iz;
	//	cout<<i<<"   "<<ix<<"  "<<iy<<endl;
      }

  }
  
  TH2F* hDeadXtalEB = NULL;
  TH2F* hDeadXtalEEm = NULL;
  TH2F* hDeadXtalEEp = NULL;

  TFile* deadXtalFile = TFile::Open(deadXtalFileName.c_str(),"READ");
  if (!deadXtalFile || !deadXtalFile->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<< deadXtalFileName <<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  hDeadXtalEB = (TH2F*) deadXtalFile->Get("rms_EB_r"); // #eta on x #phi on y
  hDeadXtalEEm = (TH2F*) deadXtalFile->Get("rms_EEm");
  hDeadXtalEEp = (TH2F*) deadXtalFile->Get("rms_EEp");

  if (!hDeadXtalEB || hDeadXtalEB == NULL) {
    cout << "Error: histogram rms_EB not found in file " << deadXtalFileName << ". End of programme." << endl;
    exit(EXIT_FAILURE);
  }

  if (!hDeadXtalEEm || hDeadXtalEEm == NULL) {
    cout << "Error: histogram rms_EEm not found in file " << deadXtalFileName << ". End of programme." << endl;
    exit(EXIT_FAILURE);
  }

  if (!hDeadXtalEEp || hDeadXtalEEp == NULL) {
    cout << "Error: histogram rms_EEp not found in file " << deadXtalFileName << ". End of programme." << endl;
    exit(EXIT_FAILURE);
  }

  // needed for EE to convert xyz to ietaRing
  TFile *EEetaRingFile = NULL;
  TH2F *hEE_etaRing = NULL;
  TH2F *hEEp_etaRing = NULL;
  TH2F *hEEm_etaRing = NULL;

  if (not isEB) {

    EEetaRingFile = new TFile((endcap_ix_iy_zside_ietaRing).c_str(),"READ");
    if (!EEetaRingFile || !EEetaRingFile->IsOpen()) {
      cout << "Error: file \"" << endcap_ix_iy_zside_ietaRing << "\" was not opened." << endl;
      exit(EXIT_FAILURE);
    }
    hEEp_etaRing = (TH2F*) EEetaRingFile->Get("hEEp");
    hEEm_etaRing = (TH2F*) EEetaRingFile->Get("hEEm");
    if (!hEEp_etaRing || hEEp_etaRing == NULL) {
      cout << "Error: histogram 'hEEp' not found in file ' " << endcap_ix_iy_zside_ietaRing << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    if (!hEEm_etaRing || hEEm_etaRing == NULL) {
      cout << "Error: histogram 'hEEm' not found in file ' " << endcap_ix_iy_zside_ietaRing << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }

  }

  TFile* f = TFile::Open(filename.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<filename<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

  f->cd();

  string directoryName = isEB ? "Barrel/" : "Endcap/";
  // bool IamInDirectory = f->cd(directoryName.c_str());
  // if (not IamInDirectory) {
  //   cout<<"Error: I couldn't change directory in file.\nApplication will be terminated."<<endl;
  //   exit(EXIT_FAILURE);
  // }

  TDirectory *dir = NULL;
  dir = f->GetDirectory(directoryName.c_str());
  if (!dir || dir == NULL) {
    cout<<"Error: I couldn't get directory in file.\nApplication will be terminated."<<endl;
    exit(EXIT_FAILURE);
  }

  dir->cd();

  UInt_t nObjectNotFound = 0;

  TH1F* hist = NULL;
  TH1F* hSum = NULL;

  bool isFirstHistogram = true;
  int nHistogramAdded = 0;
  int nTotalObjects = isEB ? 61200 : 14648;
  int nEvents = 0;

  if (plotType == "xtal") {
    
    string histName = isEB ? "epsilon_EB_iR_30003" : "epsilon_EE_iR_6397";
    hist = (TH1F*) dir->Get(histName.c_str());
    if (!hist) {
      cout << "Warning: TH1F object not found in file and plotType = " << plotType << ". Please check. Abort" <<endl;
      exit(EXIT_FAILURE);
    }
    hSum = new TH1F(*((TH1F*) hist->Clone("hSum")));

  } else {

    TIter next(dir->GetListOfKeys());
    TKey *key = NULL;

    while ( (key = (TKey*)next()) ) {

      cout.flush();
      if(nEvents % 50 == 0) cout<<"\r"<<"Crystals processed: "<<double(nEvents)/nTotalObjects*100<<" % ";
      //cout << "entry : " << nEvents << endl;                                                                                                                                
      nEvents++;

      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1F")) continue;

      hist = (TH1F*) key->ReadObj();
      if (!hist) {
	cout << "Warning: TH1F object not found in file. Skipping and going on with next object" <<endl;
	nObjectNotFound++;
	continue;
      }

      //cout << "check" << endl;

      string hname(hist->GetName());
      string hnameTag = isEB ? "epsilon_EB_iR_" : "epsilon_EE_iR_";
      if (hname.find(hnameTag.c_str()) == string::npos) continue;
    
      string fitIndexStr = ""; 
      fitIndexStr.assign(hname, hnameTag.size(), string::npos);
      //cout << "hname " << hname << "     fitIndexStr " << fitIndexStr << endl;
      Int_t fitIndex = std::stoi(fitIndexStr);

      // here we select some crystals with some algorithm
    
      bool conditionFulfilled = false;

      if (isEB) {

	int ieta = Xtal_Ieta[fitIndex];
	int iphi = Xtal_Iphi[fitIndex];

	// hardcoded, implement a flag to choose selection algorithm
	if (plotType == "etaring") {

	  // ieta == -2 and removing crystals near gaps in iphi and removing dead xtals or xtals adjacent to a dead xtals
	  conditionFulfilled = (ieta == -2 && iphi%20 != 0 && (iphi-1)%20 != 0 && noDeadXtalIn3x3matrixSeededByThisXtal(hDeadXtalEB,ieta+86,iphi));

	} else if (plotType == "xtal") {

	  if (ieta == -2 && iphi == 124) conditionFulfilled = true;

	} else {

	  conditionFulfilled = fabs(ieta) < 15 && iphi > 350 && iphi < 360;

	}

      } else {

	int ix = Xtal_Ix[fitIndex];
	int iy = Xtal_Iy[fitIndex];
	int iz = Xtal_Iz[fitIndex];


	if (plotType == "etaring") {
	
	  // etaring value is hardcoded (we used a crystal with ix,iy,iz = 19,83,-1 which is at etaRing=5, so we stick to that ring)
	  if (iz > 0) conditionFulfilled = hEEp_etaRing->GetBinContent(ix,iy) == 5 && noDeadXtalIn3x3matrixSeededByThisXtal(hDeadXtalEEp,ix,iy);
	  else        conditionFulfilled = hEEm_etaRing->GetBinContent(ix,iy) == 5 && noDeadXtalIn3x3matrixSeededByThisXtal(hDeadXtalEEm,ix,iy);

	} else if (plotType == "xtal") {

	  if (ix == 19 && iy == 83 && iz == -1) conditionFulfilled = true;

	} else {

	  int ixFromCenter = (ix < 51) ? (ix - 51) : (ix - 50); // avoid ixFromCenter = 0
	  int iyFromCenter = (iy < 51) ? (iy - 51) : (iy - 50); // avoid ixFromCenter = 0
	  double radius = sqrt(ixFromCenter*ixFromCenter + iyFromCenter*iyFromCenter);
	  //if (sqrt(radius) > 40) histToSum.push_back((TH1F*) hist->Clone());
	  if (radius > 42.0 && radius < 47.0) conditionFulfilled = true;

	}

      }

      if (conditionFulfilled) {
      
	if (isFirstHistogram) {
	  hSum = new TH1F(*((TH1F*) hist->Clone("hSum")));
	  isFirstHistogram = false;
	} else {
	  hSum->Add((TH1F*)hist->Clone());
	}
	nHistogramAdded++;
	//if (nHistogramAdded >= 50) break;
	//if (plotType == "xtal") break;
      
      }

    }

  }
  
  if (plotType != "xtal") {
    if (nHistogramAdded == 0) {
      cout << "Warning: no histogram used. End of programme." <<endl;
      exit(EXIT_FAILURE);
    }
    cout << "Selecting " << nHistogramAdded << " crystals in " << ((isEB) ? "EB" : "EE") << endl; 
  }
  if (nObjectNotFound > 0) cout << nObjectNotFound << " crystals were not found in file" << endl; 
  cout << "hSum->Integral() " << hSum->Integral() << endl;

  // it seems that the first time CMS_lumi is used the settings are screwed up 
  // produce a dummy plot (either do not save it or remove it) 
  //double lumi = 0.18; //in fb-1 
  TCanvas*ctmp = new TCanvas("ctmp","");
  ctmp->cd();
  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  htmp1->Fill(0.5);
  htmp1->Draw("H");
  CMS_lumi(ctmp,Form("%.2f",lumi),false,false);
  setTDRStyle();
  delete htmp1;
  delete ctmp;

  //h = new TH1F(hName.c_str(),"",hSum->GetNbinsX(),hSum->GetBinLowEdge(1),hSum->GetBinLowEdge(hSum->GetNbinsX()+1));
  Double_t minRange = isEB ? 0.06: 0.06;
  Double_t maxRange = isEB ? 0.22: 0.25;
  // get bin window slightly smaller than that defined by the range
  Int_t minBin = hSum->FindBin(minRange) + 1; 
  Int_t maxBin = hSum->FindBin(maxRange);
  h->SetBins(maxBin - minBin ,hSum->GetBinLowEdge(minBin),hSum->GetBinLowEdge(maxBin));
  //h->SetBins(,hSum->GetBinLowEdge(1),hSum->GetBinLowEdge(hSum->GetNbinsX()+1));
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) { 
    h->SetBinContent(i, hSum->GetBinContent(i + minBin -1));
    h->SetBinError(i, hSum->GetBinError(i + minBin -1));
  }

  if (drawAllMassPlot) drawHisto((TH1F*) h->Clone(), isEB, outDir, hName, lumi);

  f->Close();
  delete f;
  deadXtalFile->Close();
  delete deadXtalFile;
  if (not isEB) {  
    EEetaRingFile->Close();
    delete EEetaRingFile;
  }

}

//===================================

void makePi0MassWithFit(const bool isEB = true,
			const bool singleXtalOnly = false,
			const string& outDir = "/afs/cern.ch/user/m/mciprian/www/test_plot/",
			const double lumi = 8.6,
			const string& dirName = "AlCaP0_IC2017_upTo31July2017_noCC",
			const string& filePath = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/"
			) 
{

  createPlotDirAndCopyPhp(outDir);

  TH1F* h_xtal_iter0 = new TH1F();
  TH1F* h_xtal_iter2 = new TH1F();
  TH1F* h_xtal_iter6 = new TH1F();
  TH1F* h_etaring_iter0 = new TH1F();
  TH1F* h_etaring_iter2 = new TH1F();
  TH1F* h_etaring_iter6 = new TH1F();

  vectorManager* histData = new vectorManager();
  // (to be fixed) name format must be h_<type>_iter<n>, where type = "xtal" or "etaring", and n = 0,1,...
  histData->addComponent(h_xtal_iter0, "h_xtal_iter0", "1 xtal, 1st iter");
  histData->addComponent(h_xtal_iter2, "h_xtal_iter2", "1 xtal, 3rd iter");
  histData->addComponent(h_xtal_iter6, "h_xtal_iter6", "1 xtal, 7th iter");
  if (not singleXtalOnly) {
    histData->addComponent(h_etaring_iter0, "h_etaring_iter0", "#eta-ring, 1st iter");
    histData->addComponent(h_etaring_iter2, "h_etaring_iter2", "#eta-ring, 3rd iter");
    histData->addComponent(h_etaring_iter6, "h_etaring_iter6", "#eta-ring, 7th iter");
  }

  vector <TH1*> hlist(histData->getHistPtrs());
  vector<string> hname(histData->getHistNames());
  vector<string> legentry(histData->getHistLegends());


  for (uint i = 0; i < hlist.size(); i++) {

    hlist[i]->SetNameTitle(hname[i].c_str(),"");
    string niter = hname[i].substr(hname[i].find("iter")+4,string::npos);
    string filename = Form("%s%s/iter_%s/%s_epsilonPlots.root",
			   filePath.c_str(),
			   dirName.c_str(),
			   niter.c_str(),
			   dirName.c_str()); 
    string plotType = hname[i].substr(hname[i].find("h_")+2,hname[i].find("_iter")-2);
    cout << "plotType = " << plotType << endl;
    doPi0MassWithFit(hlist[i],hname[i],isEB, plotType, outDir, lumi, filename);
    if (!hlist[i] || hlist[i]==NULL) {
      cout << "Error: hlist[" << i << "] is NULL, please check. Abort" << endl;
      exit(EXIT_FAILURE);
      
    }
    
  }

  //draw_nTH1(hlist,"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.24","a.u.", canvasname, outDir, legentry, "", lumi, 1, false, false);
  if (histData->getHistPtrs().size() > 1) {

    cout << "Now going to plot all histograms together" << endl;
    string canvasname = isEB ? "pi0mass_comparison_EB" : "pi0mass_comparison_EE";
    if (singleXtalOnly) canvasname += "_singleXtal";

    if (isEB) draw_nTH1(histData->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.06,0.22","a.u.", canvasname, outDir, histData->getHistLegends(), "", lumi, 1, false, false);
    else draw_nTH1(histData->getHistPtrs(),"#gamma#gamma invariant mass (GeV/c^{2})::0.05,0.25","a.u.", canvasname, outDir, histData->getHistLegends(), "", lumi, 1, false, false);
  } else {
    cout << " Only 1 histogram was used. Not plotting more histograms together" << endl;
  }

}
