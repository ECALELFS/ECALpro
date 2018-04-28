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

#include "CMS_lumi.h"

using namespace std;

class plotManager {

public:
  plotManager(const string & histName,
              const string & xAxisName,
              const string & canvasName,
              const int    & logy = 0 // 0 for both, 1 for no log, 2 for only log                                                                                            
              )
  {
    histName_    = histName;
    xAxisName_   = xAxisName;
    canvasName_  = canvasName;
  };
  
  ~plotManager() {};
  
  string getHistName()    const { return histName_;    };
  string getXaxisName()   const { return xAxisName_;   };
  string getCanvasName()  const { return canvasName_;  };

  void setHistName   (const string & histName  )  { histName_    = histName;    };
  void setXaxisName  (const string & xAxisName )  { xAxisName_   = xAxisName;   };
  void setCanvasName (const string & canvasName)  { canvasName_  = canvasName;  };

private:

  string histName_;
  string xAxisName_;
  string canvasName_;

};



// macro to plot compare pairs of plots in the calibration epsilon files
// pass:
// 1) the name of the folder where to store the output plots
// 2) the input file name where the maps are stored (the file is probably on EOS, use root://eoscms//eos/cms/...)

//=============================================================

void drawTH1pair(TH1* h1, TH1* h2, 
		 const string& xAxisNameTmp = "", const string& yAxisName = "Events", const string& canvasName = "default", 
		 const string& outputDIR = "./", 
		 const string& legEntry1 = "data", const string& legEntry2 = "MC", const string& ratioPadYaxisName = "data/MC", 
		 const Double_t lumi = -1.0, 
		 const Bool_t drawPlotLogY = true) 
{

  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2() 

  string xAxisName = xAxisNameTmp;
  // string xAxisName = "";
  // Double_t xmin = 0;
  // Double_t xmax = 0;
  // Bool_t setXAxisRangeFromUser = getAxisRangeFromUser(xAxisName, xmin, xmax, xAxisNameTmp);

  // cout << "xAxisName = " << xAxisName << "   xmin = " << xmin << "  xmax = " << xmax << endl;

  Double_t intNum, intDen, errNum, errDen;
  intNum = h1->IntegralAndError(1,h1->GetNbinsX(),errNum);
  intDen = h2->IntegralAndError(1,h2->GetNbinsX(),errDen);
  Double_t IntegralRatio = intNum/intDen;
  //Double_t ratioError = IntegralRatio * sqrt(errNum*errNum/(intNum*intNum) + errDen*errDen/(intDen*intDen));

  if (yAxisName == "a.u.") {
    h1->Scale(1./h1->Integral());
    h2->Scale(1./h2->Integral());
  }

  h1->SetStats(0);
  h2->SetStats(0);
  h1->SetTitle("");
  h2->SetTitle("");

  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);

  TH1* frame =  (TH1*) h1->Clone("frame");
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->SetStats(0);

  h1->SetLineColor(kBlack);
  h1->SetMarkerColor(kBlack);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);

  h1->GetXaxis()->SetLabelSize(0);
  h1->GetYaxis()->SetTitle(yAxisName.c_str());
  h1->GetYaxis()->SetTitleOffset(1.1);
  // h1->GetYaxis()->SetTitleOffset(0.8);  // was 1.03 without setting also the size
  h1->GetYaxis()->SetTitleSize(0.05);
  //h1->GetYaxis()->SetRangeUser(0.0, max(h1->GetMaximum(),h2->GetMaximum()) * 1.2);
  h1->GetYaxis()->SetRangeUser(0.0, max(h1->GetBinContent(h1->GetMaximumBin()),h2->GetBinContent(h2->GetMaximumBin())) * 1.2);
  //if (setXAxisRangeFromUser) h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->Draw("EP");

  h2->SetLineColor(kRed);
  h2->SetLineWidth(2);
  h2->Draw("hist same");

  TLegend leg (0.5,0.7,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(h1,legEntry1.c_str(),"PLE");
  leg.AddEntry(h2,legEntry2.c_str(),"L");
  leg.Draw("same");
  canvas->RedrawAxis("sameaxis");

  TPaveText *pvtxt = NULL;
  if (yAxisName == "a.u.") {
    pvtxt = new TPaveText(0.5,0.6,0.90,0.7, "BR NDC");
    pvtxt->SetFillColor(0);
    pvtxt->SetFillStyle(0);
    pvtxt->SetBorderSize(0);
    pvtxt->AddText(Form("norm num/den = %.2f",IntegralRatio));
    //pvtxt->AddText(Form("norm num/den = %.2f +/- %.2f",IntegralRatio,ratioError));
    pvtxt->Draw();
  }

  //  CMS_lumi(canvas,Form("%.1f",lumi));
  if (lumi < 0) CMS_lumi(canvas,"",false,false);
  else CMS_lumi(canvas,Form("%.1f",lumi),false,false);
  setTDRStyle();

  pad2->Draw();
  pad2->cd();

  frame->Reset("ICES");
  frame->GetYaxis()->SetRangeUser(0.5,1.5);
  frame->GetYaxis()->SetNdivisions(5);
  frame->GetYaxis()->SetTitle(ratioPadYaxisName.c_str());
  frame->GetYaxis()->SetTitleOffset(1.2);
  // frame->GetYaxis()->SetTitleSize(0.15);
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle(xAxisName.c_str());
  // if (setXAxisRangeFromUser) frame->GetXaxis()->SetRangeUser(xmin,xmax);
  // frame->GetXaxis()->SetTitleOffset(0.8);
  frame->GetXaxis()->SetTitleSize(0.05);

  TH1D* ratio = (TH1D*) h1->Clone("ratio");
  TH1D* den_noerr = (TH1D*) h2->Clone("den_noerr");
  TH1D* den = (TH1D*) h2->Clone("den");
  for(int iBin = 1; iBin < den->GetNbinsX()+1; iBin++)
    den_noerr->SetBinError(iBin,0.);

  ratio->Divide(den_noerr);
  den->Divide(den_noerr);
  den->SetFillColor(kGray);
  frame->Draw();
  ratio->SetMarkerSize(0.85);
  ratio->Draw("EPsame");
  den->Draw("E2same");

  TF1* line = new TF1("horiz_line","1",ratio->GetXaxis()->GetBinLowEdge(1),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1));
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->Draw("Lsame");
  ratio->Draw("EPsame");
  pad2->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR + canvasName + ".png").c_str());
  canvas->SaveAs((outputDIR + canvasName + ".pdf").c_str());

  if (drawPlotLogY) {

    if (yAxisName == "a.u.") h1->GetYaxis()->SetRangeUser(max(0.0001,min(h1->GetMinimum(),h2->GetMinimum())*0.8),max(h1->GetMaximum(),h2->GetMaximum())*100);
    else h1->GetYaxis()->SetRangeUser(max(0.001,min(h1->GetMinimum(),h2->GetMinimum())*0.8),max(h1->GetMaximum(),h2->GetMaximum())*100);
    canvas->SetLogy();
    /* if (lumi < 0) CMS_lumi(canvas,"",true,false); */
    /* else CMS_lumi(canvas,Form("%.1f",lumi),true,false); */
    canvas->SaveAs((outputDIR + canvasName + "_logY.png").c_str());
    canvas->SaveAs((outputDIR + canvasName + "_logY.pdf").c_str());
    canvas->SetLogy(0);

  }    

  delete canvas;

}

//=============================================================                                                                                                              
  
void adjustSettings_CMS_lumi(const string& outputDir = "./") {

  // tmp plot to be removed to adjust settings in CMS_lumi                                                                                                                  
  TH1D* htmp1 = new TH1D("htmp1","",1,0,1);
  TH1D* htmp2 = new TH1D("htmp2","",1,0,1);
  htmp1->Fill(0.5);
  htmp2->Fill(0.5);
  drawTH1pair(htmp1, htmp2, "variable", "Events", "tmpToBeRemoved", outputDir);
  system(("rm " + outputDir + "*tmpToBeRemoved*").c_str());
  delete htmp1;
  delete htmp2;

}


//================================================
 

void compareDistributionsTwoFile(const string& outDir = "",			  
				 const string& inputFile = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC/iter_5/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_epsilonPlots.root",
				 const string& inputFile2 = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_nxtal9both_ext1_fromIter3/iter_1/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_nxtal9both_ext1_fromIter3_epsilonPlots.root"
			  ) 
{

  TH1::SetDefaultSumw2();

  system(Form("mkdir -p %s",outDir.c_str()));
  system(Form("cp /afs/cern.ch/user/m/mciprian/public/index.php %s",outDir.c_str()));

  vector<plotManager> myPlot;
  myPlot.push_back(plotManager("pi0pt_afterCuts_region1EB","#pi^{0} p_{T} (EB: |#eta| < 1.0)","pi0pt_region1EB"));
  myPlot.push_back(plotManager("g1pt_afterCuts_region1EB","leading (seed) #gamma p_{T} (EB: |#eta| < 1.0)","g1pt_region1EB"));
  myPlot.push_back(plotManager("g2pt_afterCuts_region1EB","trailing (seed) #gamma p_{T} (EB: |#eta| < 1.0)","g2pt_region1EB"));
  myPlot.push_back(plotManager("g1Nxtal_afterCuts_region1EB","leading (seed) #gamma Nxtal (EB: |#eta| < 1.0)","g1Nxtal_region1EB"));
  myPlot.push_back(plotManager("g2Nxtal_afterCuts_region1EB","trailing (seed) #gamma Nxtal (EB: |#eta| < 1.0)","g2Nxtal_region1EB"));

  myPlot.push_back(plotManager("pi0pt_afterCuts_region2EB","#pi^{0} p_{T} (EB: |#eta| > 1.0)","pi0pt_region2EB"));
  myPlot.push_back(plotManager("g1pt_afterCuts_region2EB","leading (seed) #gamma p_{T} (EB: |#eta| > 1.0)","g1pt_region2EB"));
  myPlot.push_back(plotManager("g2pt_afterCuts_region2EB","trailing (seed) #gamma p_{T} (EB: |#eta| > 1.0)","g2pt_region2EB"));
  myPlot.push_back(plotManager("g1Nxtal_afterCuts_region2EB","leading (seed) #gamma Nxtal(EB: |#eta| > 1.0)","g1Nxtal_region2EB"));
  myPlot.push_back(plotManager("g2Nxtal_afterCuts_region2EB","trailing (seed) #gamma Nxtal (EB: |#eta| > 1.0)","g2Nxtal_region2EB"));

  vector<TH1*> h1;
  vector<TH1*> h2;

  TFile* f1 = TFile::Open(inputFile.c_str(),"READ");
  if (!f1 || !f1->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error opening file \"" << inputFile << "\".\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }

  for (UInt_t i = 0; i < myPlot.size(); ++i) {
    h1.push_back( (TH1F*) f1->Get(myPlot[i].getHistName().c_str()) );
    if (!h1.back()) {
      cout << "Error: could not get histogram '" << myPlot[i].getHistName() << "'. End of programme" << endl;
      exit(EXIT_FAILURE);
    }
    h1.back()->SetDirectory(0);
  }
  f1->Close();


  TFile* f2 = TFile::Open(inputFile2.c_str(),"READ");
  if (!f2 || !f2->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error opening file \"" << inputFile2 << "\".\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }

  for (UInt_t i = 0; i < myPlot.size(); ++i) {
    h2.push_back( (TH1F*) f2->Get(myPlot[i].getHistName().c_str()) );
    if (!h2.back()) {
      cout << "Error: could not get histogram '" << myPlot[i].getHistName() << "'. End of programme" << endl;
      exit(EXIT_FAILURE);
    }
    h2.back()->SetDirectory(0);
  }
  f2->Close();
  

  adjustSettings_CMS_lumi(outDir);

  for (UInt_t i = 0; i < myPlot.size(); ++i) {
    drawTH1pair(h1[i], h2[i], myPlot[i].getXaxisName(), "a.u.", myPlot[i].getCanvasName(), outDir, 
		"Nxtal >= 7 (both #gamma)", "Nxtal == 9 (both #gamma)", "(>=7)/(==9)", -1.0, false); 
  }    

}

//=========================================

void compareDistributionsSameFile(const string& outDir = "",			  
				  const string& inputFile = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC/iter_3/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_epsilonPlots.root"
				  ) 
{

  TH1::SetDefaultSumw2();

  system(Form("mkdir -p %s",outDir.c_str()));
  system(Form("cp /afs/cern.ch/user/m/mciprian/public/index.php %s",outDir.c_str()));

  vector<plotManager> myPlot;
  myPlot.push_back(plotManager("pi0pt_afterCuts","#pi^{0} p_{T}","pi0pt_compRegion"));
  myPlot.push_back(plotManager("g1pt_afterCuts","leading (seed) #gamma p_{T}","g1pt_compRegion"));
  myPlot.push_back(plotManager("g2pt_afterCuts","trailing (seed) #gamma p_{T}","g2pt_compRegion"));
  myPlot.push_back(plotManager("g1Nxtal_afterCuts","leading (seed) #gamma Nxtal","g1Nxtal_compRegion"));
  myPlot.push_back(plotManager("g2Nxtal_afterCuts","trailing (seed) #gamma Nxtal","g2Nxtal_compRegion"));

  vector<TH1*> h1;
  vector<TH1*> h2;

  TFile* f1 = TFile::Open(inputFile.c_str(),"READ");
  if (!f1 || !f1->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error opening file \"" << inputFile << "\".\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }

  for (UInt_t i = 0; i < myPlot.size(); ++i) {
    h1.push_back( (TH1F*) f1->Get(Form("%s_region1EB",myPlot[i].getHistName().c_str()) ) );
    h2.push_back( (TH1F*) f1->Get(Form("%s_region2EB",myPlot[i].getHistName().c_str()) ) );
    if (!h1.back() || !h2.back()) {
      cout << "Error: could not get histogram '" << myPlot[i].getHistName() << "'. End of programme" << endl;
      exit(EXIT_FAILURE);
    }
    h1.back()->SetDirectory(0);
    h2.back()->SetDirectory(0);
  }
  f1->Close();  

  adjustSettings_CMS_lumi(outDir);

  for (UInt_t i = 0; i < myPlot.size(); ++i) {
    drawTH1pair(h1[i], h2[i], myPlot[i].getXaxisName(), "a.u.", myPlot[i].getCanvasName(), outDir, 
		"EB: |#eta| < 1.0", "EB: |#eta| > 1.0", "in/out EB", -1.0, false); 
  }    

}



//==========================================

void compareDistributions(const string& outDir = "",			  
			  const Bool_t fromSameFile = true,
			  const string& inputFile = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC/iter_5/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_epsilonPlots.root",
			  const string& inputFile2 = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_nxtal9both_ext1_fromIter3/iter_1/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_nxtal9both_ext1_fromIter3_epsilonPlots.root"		
			  ) 
{

  if (fromSameFile) compareDistributionsSameFile(outDir, inputFile);
  else              compareDistributionsTwoFile(outDir, inputFile, inputFile2);

}
