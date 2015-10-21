// -*- C++ -*-
//
// Package:    FitEpsilonPlot
// Class:      FitEpsilonPlot
// 
/**\class FitEpsilonPlot FitEpsilonPlot.cc CalibCode/FitEpsilonPlot/src/FitEpsilonPlot.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Marco Grassi, CMS
//         Created:  Tue Nov  8 17:18:54 CET 2011
// $Id: FitEpsilonPlot.cc,v 1.9 2013/06/17 13:40:42 lpernie Exp $
//
//

// system include files
#include <memory>
#include <iostream>
#include <string>

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
#include "TMath.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

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
#include "RooMinuit.h"

#include "CalibCode/FitEpsilonPlot/interface/FitEpsilonPlot.h"

using std::cout;
using std::endl;
using std::string;
using namespace RooFit;


FitEpsilonPlot::FitEpsilonPlot(const edm::ParameterSet& iConfig)

{
    /// to be moved in parameters.py
    useMassInsteadOfEpsilon_ = 1;

    //now do what ever initialization is needed
    currentIteration_ =  iConfig.getUntrackedParameter<int>("CurrentIteration");
    epsilonPlotFileName_ = iConfig.getUntrackedParameter<std::string>("EpsilonPlotFileName");
    outputDir_ = iConfig.getUntrackedParameter<std::string>("OutputDir");
    outfilename_          = iConfig.getUntrackedParameter<std::string>("OutputFile");
    calibMapPath_ = iConfig.getUntrackedParameter<std::string>("calibMapPath");
    inRangeFit_ = iConfig.getUntrackedParameter<int>("NInFit");
    finRangeFit_ = iConfig.getUntrackedParameter<int>("NFinFit");    
    EEoEB_ = iConfig.getUntrackedParameter<std::string>("EEorEB");
    isNot_2010_ = iConfig.getUntrackedParameter<bool>("isNot_2010");
    Are_pi0_ = iConfig.getUntrackedParameter<bool>("Are_pi0");
    StoreForTest_ = iConfig.getUntrackedParameter<bool>("StoreForTest","false");
    Barrel_orEndcap_ = iConfig.getUntrackedParameter<std::string>("Barrel_orEndcap");

    /// setting calibration type
    calibTypeString_ = iConfig.getUntrackedParameter<std::string>("CalibType");
    if(     calibTypeString_.compare("xtal")    == 0 ) { calibTypeNumber_ = xtal;    regionalCalibration_ = &xtalCalib; } 
    else if(calibTypeString_.compare("tt")      == 0 ) { calibTypeNumber_ = tt;      regionalCalibration_ = &TTCalib;   }
    else if(calibTypeString_.compare("etaring") == 0 ) { calibTypeNumber_ = etaring; regionalCalibration_ = &etaCalib;  }
    else throw cms::Exception("CalibType") << "Calib type not recognized\n";
    cout << "FIT_EPSILON: crosscheck: selected type: " << regionalCalibration_->printType() << endl;

    /// retrieving calibration coefficients of the previous iteration
    char fileName[200];
    if(currentIteration_ < 0) throw cms::Exception("IterationNumber") << "Invalid negative iteration number\n";
    else if(currentIteration_ > 0)
    {
	  //sprintf(fileName,"%s/iter_%d/calibMap.root", outputDir_.c_str(), currentIteration_-1);
	  sprintf(fileName,"%s", calibMapPath_.c_str());
	  regionalCalibration_->getCalibMap()->loadCalibMapFromFile(fileName);
    }

    // load epsilon from current iter
    epsilon_EB_h = new TH1F*[regionalCalibration_->getCalibMap()->getNRegionsEB()];
    epsilon_EE_h = new TH1F*[regionalCalibration_->getCalibMap()->getNRegionsEE()];
    //sprintf(fileName,"%s/iter_%d/EcalNtp.root", outputDir_.c_str(), currentIteration_);
    //sprintf(fileName,"%s/iter_%d/%s", outputDir_.c_str(), currentIteration_, epsilonPlotFileName_.c_str());
    sprintf(fileName,"%s", epsilonPlotFileName_.c_str());
    cout << "FIT_EPSILON: FitEpsilonPlot:: loading epsilon plots from file: " << epsilonPlotFileName_ << endl;
    loadEpsilonPlot(fileName);


}


FitEpsilonPlot::~FitEpsilonPlot()
{
    delete epsilon_EB_h;
    delete epsilon_EE_h;

    if(inputEpsilonFile_->IsOpen())
	  inputEpsilonFile_->Close();
}


//
// member functions
//


void FitEpsilonPlot::loadEpsilonPlot(char *filename)
{
    char line[100];

    inputEpsilonFile_ = TFile::Open(filename);
    if(!inputEpsilonFile_) 
	  throw cms::Exception("loadEpsilonPlot") << "Cannot open file " << string(filename) << "\n"; 
    if( EEoEB_ == "Barrel" && (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ){
	  for(int iR=inRangeFit_; iR <= finRangeFit_ && iR < regionalCalibration_->getCalibMap()->getNRegionsEB(); iR++)
	  {
		sprintf(line,"Barrel/epsilon_EB_iR_%d",iR);
		epsilon_EB_h[iR] = (TH1F*)inputEpsilonFile_->Get(line);

		if(!epsilon_EB_h[iR])
		    throw cms::Exception("loadEpsilonPlot") << "Cannot load histogram " << string(line) << "\n";
		else if(!(iR%1000))
		    cout << "FIT_EPSILON: Epsilon distribution for EB region " << iR << " loaded" << endl;
	  }
    }
    else if( EEoEB_ == "Endcap" && (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ){
	  for(int jR=inRangeFit_; jR <= finRangeFit_ && jR<EEDetId::kSizeForDenseIndexing; jR++)
	  {
		sprintf(line,"Endcap/epsilon_EE_iR_%d",jR);
		epsilon_EE_h[jR] = (TH1F*)inputEpsilonFile_->Get(line);
		if(!epsilon_EE_h[jR])
		    throw cms::Exception("loadEpsilonPlot") << "Cannot load histogram " << string(line) << "\n";
		else if(!(jR%1000))
		    cout << "FIT_EPSILON: Epsilon distribution for EE region " << jR << " loaded" << endl;
	  }
    }
}



void  FitEpsilonPlot::deleteEpsilonPlot(TH1F **h, int size)
{
    for(int jR=0; jR<size; jR++)
	  delete h[jR];

    delete h;
}


void FitEpsilonPlot::saveCoefficients() 
{
    /// output file
    char fileName[200];
    sprintf(fileName,"%s/%s", outputDir_.c_str(), outfilename_.c_str());
    outfile_ = new TFile(fileName,"RECREATE");
    cout << "FIT_EPSILON: Saving Calibration Coefficients in " << endl;
    //@cout << string(fileName) << " ... ";
    cout <<"FIT_EPSILON: "<< string(fileName) << " ... ";
    if(!outfile_) throw cms::Exception("WritingOutputFile") << "It was no possible to create output file " << string(fileName) << "\n";
    outfile_->cd();

    // 2D calib map in the barrel
    TH2F* hmap_EB = new TH2F("calibMap_EB","EB calib coefficients: #eta on x, #phi on y",
		2*EBDetId::MAX_IETA+1,-EBDetId::MAX_IETA-0.5,EBDetId::MAX_IETA+0.5,
		EBDetId::MAX_IPHI, EBDetId::MIN_IPHI-0.5, EBDetId::MAX_IPHI+0.5 );
    hmap_EB->GetXaxis()->SetTitle("i#eta");
    hmap_EB->GetYaxis()->SetTitle("i#phi");
    TH2F* hmap_EEp = new TH2F("calibMap_EEp","EE+ calib coefficients",100,0.5,100.5,100,0.5,100.5);
    hmap_EEp->GetXaxis()->SetTitle("ix");
    hmap_EEp->GetYaxis()->SetTitle("iy");
    TH2F* hmap_EEm = new TH2F("calibMap_EEm","EE- calib coefficients",100,0.5,100.5,100,0.5,100.5);
    hmap_EEm->GetXaxis()->SetTitle("ix");
    hmap_EEm->GetYaxis()->SetTitle("iy");
    TH1F* hint = new TH1F("hint","Bin1: inRangeFit_ Bin2: finRangeFit_ Bin3: Barrel(0)/Endcap(1)",3,0.,3.);
    hint->SetBinContent(1,inRangeFit_);
    hint->SetBinContent(2,finRangeFit_);
    if( EEoEB_ == "Barrel" ) hint->SetBinContent(3,0);
    else                     hint->SetBinContent(3,1);
    hint->Write();

    /// filling Barrel Map
    for(int j=0; j<regionalCalibration_->getCalibMap()->getNRegionsEB(); ++j)  
    {
	  std::vector<DetId> ids = regionalCalibration_->allDetIdsInEBRegion(j);
	  for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) {
		EBDetId ebid(*iid);
		int ix = ebid.ieta()+EBDetId::MAX_IETA+1;

		float coeffValue = regionalCalibration_->getCalibMap()->coeff(*iid) > 0. ? regionalCalibration_->getCalibMap()->coeff(*iid) : 1.;
		hmap_EB->SetBinContent( ix, ebid.iphi(), coeffValue );
	  } // loop over DetId in regions
    }
    hmap_EB->SetMinimum(0.9);
    hmap_EB->SetStats(false);
    hmap_EB->Write();

    for(int jR=0; jR < regionalCalibration_->getCalibMap()->getNRegionsEE(); jR++)
    {
	  std::vector<DetId> ids =  regionalCalibration_->allDetIdsInEERegion(jR);
	  for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) 
	  { 
		EEDetId eeid(*iid);
		float coeffValue =  regionalCalibration_->getCalibMap()->coeff(*iid) > 0. ?  regionalCalibration_->getCalibMap()->coeff(*iid) : 1.;

		if(eeid.positiveZ())
		    hmap_EEp->Fill(eeid.ix(), eeid.iy(), coeffValue); 
		else 
		    hmap_EEm->Fill(eeid.ix(), eeid.iy(), coeffValue);
	  }
    }

    hmap_EEp->SetMinimum(0.9);
    hmap_EEp->SetStats(false);
    hmap_EEp->Write();

    hmap_EEm->SetMinimum(0.9);
    hmap_EEm->SetStats(false);
    hmap_EEm->Write();

    /*------------- TTREE --------------*/

    uint32_t   rawId;
    int        hashedIndex;
    int        ieta;
    int        iphi;
    int        iSM;
    int        iMod;
    int        iTT;
    int        iTTeta;
    int        iTTphi;
    int        iter = currentIteration_;
    float      regCoeff;
    float      Signal;//#
    float      Backgr; 
    float      Chisqu; 
    float      Ndof; 
    float      fit_mean;
    float      fit_mean_err;
    float      fit_sigma;
    float      fit_Snorm;
    float      fit_b0;
    float      fit_b1;    
    float      fit_b2;    
    float      fit_b3;    
    float      fit_Bnorm; 
    /// endcap variables
    int ix;
    int iy;
    int zside;
    int sc; 
    int isc;
    int ic;
    int iquadrant;

    TTree* treeEB = new TTree("calibEB","Tree of EB Inter-calibration constants");
    TTree* treeEE = new TTree("calibEE","Tree of EE Inter-calibration constants");


    /// barrel
    treeEB->Branch("rawId",&rawId,"rawId/i");
    treeEB->Branch("hashedIndex",&hashedIndex,"hashedIndex/I");
    treeEB->Branch("ieta",&ieta,"ieta/I");
    treeEB->Branch("iphi",&iphi,"iphi/I");
    treeEB->Branch("iSM",&iSM,"iSM/I");
    treeEB->Branch("iMod",&iMod,"iMod/I");
    treeEB->Branch("iTT",&iTT,"iTT/I");
    treeEB->Branch("iTTeta",&iTTeta,"iTTeta/I");
    treeEB->Branch("iTTphi",&iTTphi,"iTTphi/I");
    treeEB->Branch("iter",&iter,"iter/I");
    treeEB->Branch("coeff",&regCoeff,"coeff/F");
    treeEB->Branch("Signal",&Signal,"Signal/F");//#
    treeEB->Branch("Backgr",&Backgr,"Backgr/F");
    treeEB->Branch("Chisqu",&Chisqu,"Chisqu/F");
    treeEB->Branch("Ndof",&Ndof,"Ndof/F");
    treeEB->Branch("fit_mean",&fit_mean,"fit_mean/F");
    treeEB->Branch("fit_mean_err",&fit_mean_err,"fit_mean_err/F");
    treeEB->Branch("fit_sigma",&fit_sigma,"fit_sigma/F");
    treeEB->Branch("fit_Snorm",&fit_Snorm,"fit_Snorm/F");
    treeEB->Branch("fit_b0",&fit_b0,"fit_b0/F");
    treeEB->Branch("fit_b1",&fit_b1,"fit_b1/F");
    treeEB->Branch("fit_b2",&fit_b2,"fit_b2/F");
    treeEB->Branch("fit_b3",&fit_b3,"fit_b3/F");
    treeEB->Branch("fit_Bnorm",&fit_Bnorm,"fit_Bnorm/F");

    /// endcap
    treeEE->Branch("ix",&ix,"ix/I");
    treeEE->Branch("iy",&iy,"iy/I");
    treeEE->Branch("zside",&zside,"zside/I");
    treeEE->Branch("sc",&sc,"sc/I");
    treeEE->Branch("isc",&isc,"isc/I");
    treeEE->Branch("ic",&ic,"ic/I");
    treeEE->Branch("iquadrant",&iquadrant,"iquadrant/I");
    treeEE->Branch("hashedIndex",&hashedIndex,"hashedIndex/I");
    treeEE->Branch("iter",&iter,"iter/I");
    treeEE->Branch("coeff",&regCoeff,"coeff/F");
    treeEE->Branch("Signal",&Signal,"Signal/F");//#
    treeEE->Branch("Backgr",&Backgr,"Backgr/F");
    treeEE->Branch("Chisqu",&Chisqu,"Chisqu/F");
    treeEE->Branch("Ndof",&Ndof,"Ndof/F");
    treeEE->Branch("fit_mean",&fit_mean,"fit_mean/F");
    treeEE->Branch("fit_mean_err",&fit_mean_err,"fit_mean_err/F");
    treeEE->Branch("fit_sigma",&fit_sigma,"fit_sigma/F");
    treeEE->Branch("fit_Snorm",&fit_Snorm,"fit_Snorm/F");
    treeEE->Branch("fit_b0",&fit_b0,"fit_b0/F");
    treeEE->Branch("fit_b1",&fit_b1,"fit_b1/F");
    treeEE->Branch("fit_b2",&fit_b2,"fit_b2/F");
    treeEE->Branch("fit_b3",&fit_b3,"fit_b3/F");
    treeEE->Branch("fit_Bnorm",&fit_Bnorm,"fit_Bnorm/F");


    for(int iR=0; iR < regionalCalibration_->getCalibMap()->getNRegionsEB(); ++iR)  {
	  std::vector<DetId> ids = regionalCalibration_->allDetIdsInEBRegion(iR);
	  for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) {
		EBDetId ebid(*iid);
		hashedIndex = ebid.hashedIndex();
		ieta = ebid.ieta();
		iphi = ebid.iphi();
		iSM = ebid.ism();
		iMod = ebid.im();
		iTT  = ebid.tower().hashedIndex();
		iTTeta = ebid.tower_ieta();
		iTTphi = ebid.tower_iphi();
		Signal = EBmap_Signal[ebid.hashedIndex()];//#
		Backgr = EBmap_Backgr[ebid.hashedIndex()];
		Chisqu = EBmap_Chisqu[ebid.hashedIndex()];
		Ndof = EBmap_ndof[ebid.hashedIndex()];
		fit_mean     = EBmap_mean[ebid.hashedIndex()];
		fit_mean_err = EBmap_mean_err[ebid.hashedIndex()];
		fit_sigma  = EBmap_sigma[ebid.hashedIndex()];
		fit_Snorm  = EBmap_Snorm[ebid.hashedIndex()];
		fit_b0     = EBmap_b0[ebid.hashedIndex()];
		fit_b1     = EBmap_b1[ebid.hashedIndex()];
		fit_b2     = EBmap_b2[ebid.hashedIndex()];
		fit_b3     = EBmap_b3[ebid.hashedIndex()];
		fit_Bnorm  = EBmap_Bnorm[ebid.hashedIndex()];

		regCoeff = regionalCalibration_->getCalibMap()->coeff(*iid);

		treeEB->Fill();
	  } // loop over DetId in regions
    } // loop over regions

    for(int jR=0; jR < regionalCalibration_->getCalibMap()->getNRegionsEE() ; jR++)
    {
	  std::vector<DetId> ids = regionalCalibration_->allDetIdsInEERegion(jR);
	  for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) 
	  { 
		EEDetId eeid(*iid);
		ix = eeid.ix();
		iy = eeid.iy();
		zside = eeid.zside();
		sc = eeid.sc();
		isc = eeid.isc();
		ic = eeid.ic();
		iquadrant = eeid.iquadrant();
		hashedIndex = eeid.hashedIndex();
		regCoeff = regionalCalibration_->getCalibMap()->coeff(*iid);
		Signal = EEmap_Signal[eeid.hashedIndex()];//#
		Backgr = EEmap_Backgr[eeid.hashedIndex()];
		Chisqu = EEmap_Chisqu[eeid.hashedIndex()];            
		Ndof = EEmap_ndof[eeid.hashedIndex()];            
		fit_mean     = EEmap_mean[eeid.hashedIndex()];
		fit_mean_err = EEmap_mean_err[eeid.hashedIndex()];
		fit_sigma  = EEmap_sigma[eeid.hashedIndex()];
		fit_Snorm  = EEmap_Snorm[eeid.hashedIndex()];
		fit_b0     = EEmap_b0[eeid.hashedIndex()];
		fit_b1     = EEmap_b1[eeid.hashedIndex()];
		fit_b2     = EEmap_b2[eeid.hashedIndex()];
		fit_b3     = EEmap_b3[eeid.hashedIndex()];
		fit_Bnorm  = EEmap_Bnorm[eeid.hashedIndex()];

		treeEE->Fill();
	  }
    }

    treeEB->Write();
    treeEE->Write();

    outfile_->Write();
    outfile_->Close();
    cout << "FIT_EPSILON:  done" << endl;
}

// ------------ method called for each event  ------------
    void
FitEpsilonPlot::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    TF1 ffit("gausa","gaus(0)+[3]*x+[4]",-0.5,0.5);
    ffit.SetParameters(100,0,0.1);
    ffit.SetParNames("Constant","Mean_value","Sigma","a","b");

    ffit.SetParLimits(3,-500,500);
    ffit.SetParLimits(2,0.05,0.22);

    cout << "FIT_EPSILON: About to fit epsilon distributions" << endl; 

    /// compute average weight, eps, and update calib constant
    if( (EEoEB_ == "Barrel") && (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ){
	  for(uint32_t j= (uint32_t)inRangeFit_; j <= (uint32_t)finRangeFit_ && j < (uint32_t)regionalCalibration_->getCalibMap()->getNRegionsEB(); ++j)  
	  {
		cout<<"FIT_EPSILON: Fitting EB Cristal--> "<<j<<endl;

		if(!(j%1000)) cout << "FIT_EPSILON: fitting EB region " << j << endl;

		float mean = 0.;
		if(!useMassInsteadOfEpsilon_ && epsilon_EB_h[j]->Integral(epsilon_EB_h[j]->GetNbinsX()*(1./6.),epsilon_EB_h[j]->GetNbinsX()*0.5) > 20) 
		{

		    double Max = 0.;
		    double Min = -0.5, bin = 0.0125;
		    Max = Min+(bin*(double)epsilon_EB_h[j]->GetMaximumBin());
		    double Bound1 = -0.15, Bound2 = 0.25;
		    if ( fabs(Max+Bound1) > 0.24  ){ Bound1 = -0.1;}
		    if ( Max+Bound2 > 0.34  ){ Bound2 = 0.15;}
		    if ( fabs(Max+Bound1) > 0.24  ){ Bound1 = -0.075;}
		    if ( Max+Bound2 > 0.34  ){ Bound2 = 0.1;}
		    if ( fabs(Max+Bound1) > 0.24  ){ Bound1 = -0.03;}
		    if ( Max+Bound2 > 0.34  ){ Bound2 = 0.05;}
		    if ( fabs(Max+Bound1) > 0.24  ){ Bound1 = -0.009;}
		    if ( Max+Bound2 > 0.34  ){ Bound2 = 0.01;}

		    epsilon_EB_h[j]->Fit(&ffit,"qB","", Max+Bound1,Max+Bound2);
		    if(ffit.GetNDF() != 0) {
			  double chi2 = ( ffit.GetChisquare()/ffit.GetNDF() );

			  if ( chi2  > 11 ){
				ffit.SetParLimits(2,0.05,0.15);
				ffit.SetParameters(100,0,0.1);
				epsilon_EB_h[j]->Fit(&ffit,"qB","", Max+Bound1,Max+Bound2);
				chi2 = (ffit.GetChisquare()/ffit.GetNDF());
				if ( chi2  < 11 ){   cout<<"Saved 1 Level!!"<<endl;  }
				else{
				    ffit.SetParameters(100,0,0.1);
				    ffit.SetParLimits(2,0.05,0.1);
				    epsilon_EB_h[j]->Fit(&ffit,"qB","",  Max+Bound1,Max+Bound2);
				    chi2 = (ffit.GetChisquare()/ffit.GetNDF());
				    if ( chi2  < 11 ){ cout<<"Saved 2 Level!!"<<endl; }
				    else{ cout<<"DAMN: High Chi square..."<<endl; }
				}
			  }
		    }
		    else cout<<"DAMN: NDF == 0"<<endl;
		    mean = ffit.GetParameter(1);
		}
		else if(useMassInsteadOfEpsilon_)
		{
		    int iMin = epsilon_EB_h[j]->GetXaxis()->FindBin(Are_pi0_? 0.08:0.4 ); 
		    int iMax = epsilon_EB_h[j]->GetXaxis()->FindBin(Are_pi0_? 0.18:0.65 );
		    double integral = epsilon_EB_h[j]->Integral(iMin, iMax);  
		    if(integral>60.)
		    {
			  Pi0FitResult fitres = FitMassPeakRooFit( epsilon_EB_h[j], Are_pi0_? 0.08:0.4, Are_pi0_? 0.21:0.65, j, 1, Pi0EB, 0, isNot_2010_); //0.05-0.3
			  RooRealVar* mean_fitresult = (RooRealVar*)(((fitres.res)->floatParsFinal()).find("mean"));
			  mean = mean_fitresult->getVal();
			  float r2 = mean/(Are_pi0_? PI0MASS:ETAMASS);
			  r2 = r2*r2;
			  //cout<<"EBMEAN::"<<j<<":"<<mean<<" Saved if: "<<fitres.SoB<<">(isNot_2010_ ? 0.04:0.1) "<<(fitres.chi2/fitres.dof)<<" < 0.2 "<<fabs(mean-0.15)<<" >0.0000001) "<<endl;
			  //if( fitres.SoB>(isNot_2010_ ? 0.04:0.1) && (fitres.chi2/fitres.dof)< 0.5 && fabs(mean-0.15)>0.0000001) mean = 0.5 * ( r2 - 1. );
			  if( fitres.chi2 < 5 && fabs(mean-0.15)>0.0000001) mean = 0.5 * ( r2 - 1. );
			  else                                              mean = 0.;
		    }
		    else{
			  mean = 0.;
		    }
		}


		std::vector<DetId> ids = regionalCalibration_->allDetIdsInEBRegion(j);
		for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) 
		{
		    regionalCalibration_->getCalibMap()->coeff(*iid) *= (mean==0.) ? 1. : 1./(1.+mean);
		} // loop over DetId in regions
	  } // loop over regions
    }// if you have to fit barrel

    /// loop over EE crystals
    if( (EEoEB_ == "Endcap") && (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ){
	  for(int jR = inRangeFit_; jR <=finRangeFit_ && jR < regionalCalibration_->getCalibMap()->getNRegionsEE(); jR++)
	  {
		cout << "FIT_EPSILON: Fitting EE Cristal--> " << jR << endl;
		if(!(jR%1000))
		    cout << "FIT_EPSILON: fitting EE region " << jR << endl;

		float mean = 0.;
		if(!useMassInsteadOfEpsilon_ && epsilon_EE_h[jR]->Integral(epsilon_EE_h[jR]->GetNbinsX()*(1./6.),epsilon_EE_h[jR]->GetNbinsX()*0.5) > 20) 
		{
		    TF1 *ffit = new TF1("gausa","gaus(0)+[3]*x+[4]",-0.5,0.5);
		    ffit->SetParameters(100,0,0.1);
		    ffit->SetParNames("Constant","Mean_value","Sigma","a","b");

		    ffit->SetParLimits(0,0.,epsilon_EE_h[jR]->GetEntries()*1.1);
		    ffit->SetParLimits(3,-500,500);
		    ffit->SetParLimits(2,0.05,0.3);

		    double Max = 0.;
		    double Min = -0.5, bin = 0.0125;
		    Max = Min+(bin*(double)epsilon_EE_h[jR]->GetMaximumBin());
		    double Bound1 = -0.35, Bound2 = 0.35;
		    if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.3;}
		    if ( Max+Bound2 > 0.48  ){ Bound2 = 0.3;}
		    if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.25;}
		    if ( Max+Bound2 > 0.48  ){ Bound2 = 0.2;}
		    if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.2;}
		    if ( Max+Bound2 > 0.48  ){ Bound2 = 0.15;}
		    if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.15;}
		    if ( Max+Bound2 > 0.48  ){ Bound2 = 0.1;}
		    if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.1;}
		    if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.05;}
		    //@@IterativeFit(epsilon_EE_h[jR], *ffit);
		    //@@mean = ffit.GetParameter(1); 
		    epsilon_EE_h[jR]->Fit(ffit,"qB","", Max+Bound1,Max+Bound2);

		    if(ffit->GetNDF() != 0) {
			  double chi2 = ( ffit->GetChisquare()/ffit->GetNDF() );
			  if(chi2 > 11  ) { cout<<"DAMN:(EE) High Chi square..."<<endl; }
		    }
		    else cout<<"DAMN: NDF == 0"<<endl;
		    mean = ffit->GetParameter(1);
		}
		else if(useMassInsteadOfEpsilon_)
		{
		    int iMin = epsilon_EE_h[jR]->GetXaxis()->FindBin(Are_pi0_? 0.08:0.4 ); 
		    int iMax = epsilon_EE_h[jR]->GetXaxis()->FindBin(Are_pi0_? 0.18:0.65 );
		    double integral = epsilon_EE_h[jR]->Integral(iMin, iMax);  

		    if(integral>70.)
		    {
			  Pi0FitResult fitres = FitMassPeakRooFit( epsilon_EE_h[jR], Are_pi0_? 0.08:0.4, Are_pi0_? 0.25:0.65, jR, 1, Pi0EE, 0, isNot_2010_);//0.05-0.3
			  RooRealVar* mean_fitresult = (RooRealVar*)(((fitres.res)->floatParsFinal()).find("mean"));
			  mean = mean_fitresult->getVal();
			  float r2 = mean/(Are_pi0_? PI0MASS:ETAMASS);
			  r2 = r2*r2;
			  //cout<<"EEMEAN::"<<jR<<":"<<mean<<" Saved if: "<<fitres.SoB<<">0.3 "<<(fitres.chi2/fitres.dof)<<" < (isNot_2010_? 0.07:0.35) "<<fabs(mean-0.14)<<" >0.0000001) "<<endl;
			  //if( (fitres.chi2/fitres.dof)<0.3 && fitres.SoB>(isNot_2010_? 0.07:0.35) && fabs(mean-0.14)>0.0000001 ) mean = 0.5 * ( r2 - 1. );
			  if( fitres.chi2 < 5 && fabs(mean-0.16)>0.0000001 ) mean = 0.5 * ( r2 - 1. );
			  else                                              mean = 0.;
		    }
		    else
		    {
			  mean = 0.; 
		    }
		}

		std::vector<DetId> ids = regionalCalibration_->allDetIdsInEERegion(jR);
		for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) 
		{
		    regionalCalibration_->getCalibMap()->coeff(*iid) *= (mean==0.) ? 1. : 1./(1.+mean);
		}
	  }//for EE
    }// if you have to fit Endcap

}



    void 
FitEpsilonPlot::IterativeFit(TH1F* h, TF1 & ffit) 
{
    float chi2 = 1000.;
    const int iterMax = 10;

    h->Fit(&ffit,"q","",-0.4,0.4);

    float mean = (ffit.GetParameters())[1];
    float sigma = (ffit.GetParameters())[2];
    float xmin = mean-2.*sigma;
    float xmax = mean+2.*sigma;

    ffit.SetRange(xmin,xmax);

    double par[3] = { ffit.GetParameters()[0], mean, sigma };

    for(int iter=0; iter< iterMax && chi2>5.; ++iter) 
    {
	  ffit.SetParameters(par[0],par[1], par[2]);

	  h->Fit(&ffit,"q","",xmin,xmax);
	  par[0] = (ffit.GetParameters())[0];
	  par[1] = (ffit.GetParameters())[1];
	  par[2] = (ffit.GetParameters())[2];

	  if(ffit.GetNDF()!=0) {
		chi2 = ffit.GetChisquare()/ffit.GetNDF();
	  }

    }
    return;
}



//-----------------------------------------------------------------------------------
Pi0FitResult FitEpsilonPlot::FitMassPeakRooFit(TH1F* h, double xlo, double xhi,  uint32_t HistoIndex, int ngaus, FitMode mode, int niter, bool isNot_2010_) 
{
    //-----------------------------------------------------------------------------------


    RooRealVar x("x","#gamma#gamma invariant mass",xlo, xhi, "GeV/c^2");

    RooDataHist dh("dh","#gamma#gamma invariant mass",RooArgList(x),h);

    RooRealVar mean("mean","#pi^{0} peak position", Are_pi0_? 0.13:0.52,  Are_pi0_? 0.105:0.5, Are_pi0_? 0.15:0.62,"GeV/c^{2}");
    RooRealVar sigma("sigma","#pi^{0} core #sigma",0.013, 0.005,0.020,"GeV/c^{2}");


    if(mode==Pi0EE)  {
	  mean.setRange( Are_pi0_? 0.1:0.45, Are_pi0_? 0.16:0.62);
	  mean.setVal(Are_pi0_? 0.13:0.55);
	  sigma.setRange(0.005, 0.060);
    }
    if(mode==Pi0EB && niter==1){
	  mean.setRange(Are_pi0_? 0.105:0.47, Are_pi0_? 0.15:0.62);
	  sigma.setRange(0.003, 0.030);
    }

    RooRealVar Nsig("Nsig","#pi^{0} yield",1000.,0.,1.e7);
    Nsig.setVal( h->GetSum()*0.1);

    RooGaussian gaus("gaus","Core Gaussian",x, mean,sigma);

    RooRealVar sigmaTail("sigmaTail","#pi^{0} tail #sigma",0.040, 0.020,0.065,"GeV/c^{2}");

    RooGaussian gaus2("gaus2","Tail Gaussian",x, mean,sigmaTail);

    RooRealVar fcore("fcore","f_{core}",0.9,0.,1.);
    RooAddPdf  signal("signal","signal model",RooArgList(gaus,gaus2),fcore);

    RooRealVar p0("p0","p0", 1000.,-1.e5,1.e5);
    RooRealVar p1("p1","p1", -3000.,-1.e5,1.e5);
    RooRealVar p2("p2","p2", 10000.,-1.e5,1.e5);
    RooRealVar p3("p3","p3", -10000.,-1.e5,1.e5);
    RooRealVar p4("p4","p4",-4000.,-1.e5,1.e5);
    RooRealVar p5("p5","p5", 5.,-1.e5,1.e5);
    RooRealVar p6("p6","p6", 6.,-1.e5,1.e5);

    RooRealVar cb0("cb0","cb0", 0.2, -1.,1.);
    RooRealVar cb1("cb1","cb1",-0.1, -1.,1.);
    RooRealVar cb2("cb2","cb2", 0.1, -1.,1.);
    RooRealVar cb3("cb3","cb3",-0.1, -0.5,0.5);
    RooRealVar cb4("cb4","cb4", 0.1, -1.,1.);
    RooRealVar cb5("cb5","cb5", 0.1, -1.,1.);
    RooRealVar cb6("cb6","cb6", 0.3, -1.,1.);


    //RooChebychev bkg("bkg","bkg model", x, RooArgList(cb0,cb1,cb2) );
    //RooChebychev bkg("bkg","bkg model", x, RooArgList(cb0,cb1,cb2,cb3) );

    RooArgList cbpars(cb0,cb1,cb2);
    if(mode==Pi0EB || mode==Pi0EE) cbpars.add( cb3);
    //if(mode==Pi0EE) cbpars.add( cb4);
    //if(mode==Pi0EE) cbpars.add( cb5);

    if(mode==Pi0EB && niter==1){
	  cb3.setRange(-1,1.);
	  cb4.setRange(-0.3,0.3);
	  cbpars.add( cb4 );
    }
    if(mode==Pi0EB && niter==2){
	  cb3.setRange(-1,1.);
	  cb4.setRange(-1,1);
	  cbpars.add( cb4 );
    }
    if(mode==Pi0EB && niter==3){
	  cb3.setRange(-1,1.);
	  cb4.setRange(-1,1);
	  cb5.setRange(-0.5, 0.5);
	  cbpars.add( cb4 );
	  cbpars.add( cb5 );
    }

    RooChebychev bkg("bkg","bkg model", x, cbpars );

    //RooPolynomial bkg("bkg","background model",x,RooArgList(p0,p1,p2,p3,p4,p5,p6) );
    //RooPolynomial bkg("bkg","background model",x,RooArgList(p0,p1,p2,p3) );

    RooRealVar Nbkg("Nbkg","background yield",1.e3,0.,1.e8);
    Nbkg.setVal( h->GetSum()*0.8 );

    RooAbsPdf* model=0;

    RooAddPdf model1("model","sig+bkg",RooArgList(gaus,bkg),RooArgList(Nsig,Nbkg));
    RooAddPdf model2("model","sig+bkg",RooArgList(signal,bkg),RooArgList(Nsig,Nbkg));

    if(ngaus==1)      model = &model1;
    else if(ngaus==2) model = &model2;


    RooNLLVar nll("nll","log likelihood var",*model,dh, RooFit::Extended(true));
    //RooAbsReal * nll = model->createNLL(dh); //suggetsed way, taht should be the same
    RooMinuit m(nll);
    m.setVerbose(kFALSE);
    //m.setVerbose(kTRUE);
    m.migrad();
    //m.hesse();
    RooFitResult* res = m.save() ;

    RooChi2Var chi2("chi2","chi2 var",*model,dh, true);
    int ndof = h->GetNbinsX() - res->floatParsFinal().getSize();

    // compute S/B

    //compute S/B and chi2
    x.setRange("sobRange",mean.getVal()-3.*sigma.getVal(), mean.getVal()+3.*sigma.getVal());
    RooAbsReal* integralSig = gaus.createIntegral(x,NormSet(x),Range("sobRange"));

    RooAbsReal* integralBkg = bkg.createIntegral(x,NormSet(x),Range("sobRange"));

    float normSig = integralSig->getVal();
    float normBkg = integralBkg->getVal();

    Pi0FitResult pi0res; // this is the output value of this method
    pi0res.res = res;

    pi0res.S = normSig*Nsig.getVal();
    pi0res.Serr = normSig*Nsig.getError();

    pi0res.B = normBkg*Nbkg.getVal();
    pi0res.Berr = normBkg*Nbkg.getError();

    pi0res.SoB =  pi0res.S/pi0res.B;
    pi0res.SoBerr =  pi0res.SoB*sqrt( pow(pi0res.Serr/pi0res.S,2) + 
		pow(pi0res.Berr/pi0res.B,2) ) ;
    pi0res.dof = ndof;

    RooPlot*  xframe = x.frame(h->GetNbinsX());
    xframe->SetTitle(h->GetTitle());
    dh.plotOn(xframe);
    model->plotOn(xframe,Components(bkg),LineStyle(kDashed), LineColor(kRed));
    model->plotOn(xframe);

    xframe->Draw();
    pi0res.probchi2 = TMath::Prob(xframe->chiSquare(), ndof);
    pi0res.chi2 = xframe->chiSquare();
    cout << "FIT_EPSILON: Nsig: " << Nsig.getVal() 
	  << " nsig 3sig: " << normSig*Nsig.getVal()
	  << " nbkg 3sig: " << normBkg*Nbkg.getVal()
	  << " S/B: " << pi0res.SoB << " +/- " << pi0res.SoBerr
	  << " chi2: " << xframe->chiSquare()
	  << " DOF: " << pi0res.dof
	  << " prob(chi2): " << pi0res.probchi2
					<< endl;

    if(mode==Pi0EB){
	  EBmap_Signal[HistoIndex]=pi0res.S;
	  EBmap_Backgr[HistoIndex]=pi0res.B;
	  EBmap_Chisqu[HistoIndex]=xframe->chiSquare();
	  EBmap_ndof[HistoIndex]=ndof;
	  EBmap_mean[HistoIndex]=mean.getVal();
	  EBmap_mean_err[HistoIndex]=mean.getError();
	  EBmap_sigma[HistoIndex]=sigma.getVal();
	  EBmap_Snorm[HistoIndex]=normSig;
	  EBmap_b0[HistoIndex]=cb0.getVal();
	  EBmap_b1[HistoIndex]=cb1.getVal();
	  EBmap_b2[HistoIndex]=cb2.getVal();
	  EBmap_b3[HistoIndex]=cb3.getVal();
	  EBmap_Bnorm[HistoIndex]=normBkg;
    }
    if(mode==Pi0EE){
	  EEmap_Signal[HistoIndex]=pi0res.S;
	  EEmap_Backgr[HistoIndex]=pi0res.B;
	  EEmap_Chisqu[HistoIndex]=xframe->chiSquare();
	  EEmap_ndof[HistoIndex]=ndof;
	  EEmap_mean[HistoIndex]=mean.getVal();
	  EEmap_mean_err[HistoIndex]=mean.getError();
	  EEmap_sigma[HistoIndex]=sigma.getVal();
	  EEmap_Snorm[HistoIndex]=normSig;
	  EEmap_b0[HistoIndex]=cb0.getVal();
	  EEmap_b1[HistoIndex]=cb1.getVal();
	  EEmap_b2[HistoIndex]=cb2.getVal();
	  EEmap_b3[HistoIndex]=cb3.getVal();
	  EEmap_Bnorm[HistoIndex]=normBkg;
    }

    TLatex lat;
    char line[300];
    lat.SetNDC();
    lat.SetTextSize(0.040);
    lat.SetTextColor(1);

    float xmin(0.55), yhi(0.80), ypass(0.05);
    if(mode==EtaEB) yhi=0.30;
    sprintf(line,"Yield: %.0f #pm %.0f", Nsig.getVal(), Nsig.getError() );
    lat.DrawLatex(xmin,yhi, line);

    sprintf(line,"m_{#gamma#gamma}: %.2f #pm %.2f", mean.getVal()*1000., mean.getError()*1000. );
    lat.DrawLatex(xmin,yhi-ypass, line);

    sprintf(line,"#sigma: %.2f #pm %.2f (%.2f%s)", sigma.getVal()*1000., sigma.getError()*1000., sigma.getVal()*100./mean.getVal(), "%" );
    lat.DrawLatex(xmin,yhi-2.*ypass, line);

    //sprintf(line,"S/B(3#sigma): %.2f #pm %.2f", pi0res.SoB, pi0res.SoBerr );
    sprintf(line,"S/B(3#sigma): %.2f", pi0res.SoB );
    lat.DrawLatex(xmin,yhi-3.*ypass, line);

    sprintf(line,"#Chi^2: %.2f", xframe->chiSquare()/pi0res.dof );
    lat.DrawLatex(xmin,yhi-4.*ypass, line);

    sprintf(line,"Attempt: %d", niter );
    lat.DrawLatex(xmin,yhi-5.*ypass, line);

    Pi0FitResult fitres = pi0res;
    //if(mode==Pi0EB && ( xframe->chiSquare()/pi0res.dof>0.35 || pi0res.SoB<0.6 || fabs(mean.getVal()-(Are_pi0_? 0.150:0.62))<0.0000001 ) ){
    if(mode==Pi0EB && ( xframe->chiSquare()>5 || fabs(mean.getVal()-(Are_pi0_? 0.150:0.62))<0.0000001 ) ){
	  if(niter==0) fitres = FitMassPeakRooFit( h, xlo, xhi, HistoIndex, ngaus, mode, 1, isNot_2010_);
	  if(niter==1) fitres = FitMassPeakRooFit( h, xlo, xhi, HistoIndex, ngaus, mode, 2, isNot_2010_);
	  if(niter==2) fitres = FitMassPeakRooFit( h, xlo, xhi, HistoIndex, ngaus, mode, 3, isNot_2010_);
    }
    if(StoreForTest_ && niter==0){
	  std::stringstream ind;
	  ind << (int) HistoIndex;
	  TString nameHistofit = "Fit_n_" + ind.str();
	  xframe->SetName(nameHistofit.Data());
	  outfileTEST_->cd();
	  xframe->Write();
    }

    return fitres;
}

// ------------ method called once each job just before starting event loop  ------------
    void 
FitEpsilonPlot::beginJob()
{
    if(StoreForTest_){
	  outfileTEST_ = new TFile("/tmp/Fit_Stored.root","RECREATE");
	  if(!outfileTEST_) cout<<"WARNING: file with fit not created."<<endl;
    }
}

// ------------ method called once each job just after ending the event loop  ------------
    void 
FitEpsilonPlot::endJob() 
{
    saveCoefficients();
    if(StoreForTest_){
	  cout<<"Fit stored in /tmp/Fit_Stored.root"<<endl;
	  outfileTEST_->Write();
	  outfileTEST_->Close();
    }
}

// ------------ method called when starting to processes a run  ------------
    void 
FitEpsilonPlot::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
    void 
FitEpsilonPlot::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
    void 
FitEpsilonPlot::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
    void 
FitEpsilonPlot::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FitEpsilonPlot::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FitEpsilonPlot);
