#include <TROOT.h>   
#include <TAttFill.h>
#include <TAxis.h>  
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>  
#include <TFile.h>
#include <TFitResult.h>  
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TVirtualFitter.h>

#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h                                
#include <cstdio>
#include <cmath>
#include <sstream> //to use ostringstream to convert numbers to string in c++          

#include "/afs/cern.ch/work/m/mciprian/w_mass_analysis/CMSSW_8_0_25/src/WmassAnalysis/macros/utility.h"

using namespace std;

void drawControlHistogramsEoverEtrue(const string& inputFile = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/pi0Gun_MC_EoverEtrue_foldSM/iter_0/pi0Gun_MC_EoverEtrue_foldSM_epsilonPlots.root", const string& outDir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/CC_EoverEtrue/") {

  system(Form("mkdir -p %s",outDir.c_str()));
  system(Form("cp /afs/cern.ch/user/m/mciprian/public/index.php %s",outDir.c_str()));

  TFile* f = TFile::Open(inputFile.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error opening file \"" << inputFile << "\".\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }

  vector<string> histoName;  
  vector<string> histoName2D;

  vector<string> regId;
  regId.push_back("region1EB");
  regId.push_back("region2EB");
  regId.push_back("region1EE");
  regId.push_back("region2EE");
  vector<string> etaRange;
  etaRange.push_back("EB |#eta| < 1.0");
  etaRange.push_back("EB |#eta| > 1.0");
  etaRange.push_back("EE |#eta| < 1.8");
  etaRange.push_back("EE |#eta| > 1.8");


  vector<plotManager> myPlot;
  for (UInt_t i = 0; i  < regId.size(); ++i) {

    myPlot.push_back(plotManager(Form("pi0pt_afterCuts_%s",regId[i].c_str()),"p_{T}(#pi^{0}) [GeV]",Form("pi0pt_afterCuts_%s",regId[i].c_str()),1,1));
    myPlot.push_back(plotManager(Form("g1pt_afterCuts_%s",regId[i].c_str()),"p_{T}(#gamma_{1}) [GeV]",Form("g1pt_afterCuts_%s",regId[i].c_str()),1,1));
    myPlot.push_back(plotManager(Form("g2pt_afterCuts_%s",regId[i].c_str()),"p_{T}(#gamma_{2}) [GeV]",Form("g2pt_afterCuts_%s",regId[i].c_str()),1,1));
    myPlot.push_back(plotManager(Form("g1Nxtal_afterCuts_%s",regId[i].c_str()),"Number of crystals (#gamma_{1})",Form("g1Nxtal_afterCuts_%s",regId[i].c_str()),1,1));
    myPlot.push_back(plotManager(Form("g2Nxtal_afterCuts_%s",regId[i].c_str()),"Number of crystals (#gamma_{2})",Form("g2Nxtal_afterCuts_%s",regId[i].c_str()),1,1));
    myPlot.push_back(plotManager(Form("pi0PhotonsNoverlappingXtals_afterCuts_%s",regId[i].c_str()),"Number of overlapping crystals",Form("pi0PhotonsNoverlappingXtals_afterCuts_%s",regId[i].c_str()),1,1));

    histoName2D.push_back(Form("pi0MassVsPU_%s",regId[i].c_str()));

  }  

  histoName.push_back("h_numberUnmergedGenPhotonPairs_EB");
  histoName.push_back("h_numberUnmergedGenPhotonPairs_EE");
  histoName.push_back("h_numberUnmergedGenPhotonPairs");
  histoName.push_back("h_numberMatchedGenPhotonPairs_EB");
  histoName.push_back("h_numberMatchedGenPhotonPairs_EE");
  histoName.push_back("h_numberMatchedGenPhotonPairs");

  ///////////////////////////////////
  // tmp plot to be removed to adjust settings in CMS_lumi                                           
  ///////////////////////
  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  TH1D* htmp2 = new TH1D("htmp2","",1,0,1);
  htmp1->Fill(0.5);
  htmp2->Fill(0.5);
  vector<TH1*> htmpVec; htmpVec.push_back(htmp2);
  drawTH1dataMCstack(htmp1, htmpVec, "variable", "Events", "tmpToBeRemoved", outDir);
  system(("rm " + outDir + "*tmpToBeRemoved*").c_str());
  ///////////////////////////////



  for (UInt_t i = 0; i < myPlot.size(); ++i) {

    TH1F* htmp = nullptr;
    
    htmp  = (TH1F*) f->Get(myPlot[i].getHistName().c_str());
    if (!htmp) {
      cout << "Error: could not get histogram named " << myPlot[i].getHistName() << ". End of programme" << endl;
      exit(EXIT_FAILURE);
    }
    htmp->SetDirectory(0);

    Int_t etaRangeIndex = i/(myPlot.size()/etaRange.size()); 
    drawSingleTH1(htmp,myPlot[i].getXaxisName(),"Events",myPlot[i].getCanvasName(),outDir,etaRange[etaRangeIndex],-1.0,myPlot[i].getRebinFactor(),false,myPlot[i].getLogy());

  }

  for (UInt_t i = 0; i < histoName.size(); ++i) {

    TH1F* htmp = nullptr;
    
    htmp  = (TH1F*) f->Get(histoName[i].c_str());
    if (!htmp) {
      cout << "Error: could not get histogram named " << histoName[i] << ". End of programme" << endl;
      exit(EXIT_FAILURE);
    }
    htmp->SetDirectory(0);

    TCanvas* canvas = new TCanvas("canvas","");
    canvas->cd();
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    canvas->cd();
    canvas->SetRightMargin(0.06);
    htmp->Draw("HIST");
    htmp->SetStats(0);
    htmp->GetXaxis()->SetTitleSize(0.05);
    htmp->GetXaxis()->SetTitleOffset(0.7);
    htmp->GetXaxis()->SetTitle(htmp->GetTitle());
    htmp->GetYaxis()->SetTitleSize(0.06);
    htmp->GetYaxis()->SetTitleOffset(0.7);
    htmp->GetYaxis()->SetTitle("Events");
    canvas->SaveAs(Form("%s/%s.png",outDir.c_str(),htmp->GetName()));
    canvas->SaveAs(Form("%s/%s.pdf",outDir.c_str(),htmp->GetName()));
    delete canvas;

  }

  for (UInt_t i = 0; i < histoName2D.size(); ++i) {

    TH2F* htmp = nullptr;
    
    htmp  = (TH2F*) f->Get(histoName2D[i].c_str());
    if (!htmp) {
      cout << "Error: could not get histogram named " << histoName2D[i] << ". End of programme" << endl;
      exit(EXIT_FAILURE);
    }
    htmp->SetDirectory(0);

    vector<Int_t> PUrangeForProfile = {1,15,20,25,30,35,51}; // range includes lower value but not upper value (1->14, 15->19,20->24,...) 
    vector<TH1*> hproj_PU;
    vector<string> legEntries;
    vector<Int_t> PUbinToDraw = {1,0,1,0,0,1}; 

    for (UInt_t jPU = 0; jPU < (PUrangeForProfile.size() - 1); ++jPU) {

      if (PUbinToDraw[jPU]) {
	hproj_PU.push_back( htmp->ProjectionX(Form("%s_projX_PU_%dto%d",htmp->GetName(),PUrangeForProfile[jPU],PUrangeForProfile[jPU+1]-1),
					      htmp->GetYaxis()->FindFixBin(PUrangeForProfile[jPU]),
					      htmp->GetYaxis()->FindFixBin(PUrangeForProfile[jPU+1]-1)
					      ) 
			    );
	legEntries.push_back(Form("%d < PU < %d",PUrangeForProfile[jPU],PUrangeForProfile[jPU+1]-1));
      }      

    }

    draw_nTH1(hproj_PU, "#pi^{0} mass [GeV]::0.08,0.21", "a.u.", histoName2D[i]+"_compareProjection", outDir, legEntries, "", -1.0, 1, false, false);

    drawCorrelationPlot(htmp,"#pi^{0} mass [GeV]","number of true PU events","Events",histoName2D[i],"",outDir,
			2,1,false,false,false,1);

  }

  f->Close();


} 
