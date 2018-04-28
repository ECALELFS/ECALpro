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
#include <TProfile2D.h>
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
#include <fstream>

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

#include "./ICfunctions.h"

using namespace std;

// macro to plot maps of IC folded in one SM (actually averaging the IC in the 36 SM)
// pass:
// 1) the name of the folder where to store the output plots
// 2) the input file name where the maps are stored (the file is probably on EOS, use root://eoscms//eos/cms/...)

//====================================================

void foldBarrelMapInSM(const string& outDir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC/iter_6/2DMaps/ICmaps/norm1etaring/",
		       const string& inputFile = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC/iter_6/2DMaps/ICmaps/norm1etaring/IC_fromECALpro_EcalBarrel_ic_2d_norm1etaring.root",
		       const string& mapNameInFile = "IC_fromECALpro_EcalBarrel_ic_2d_norm1etaring",
		       const Double_t mapMin = 0.95,
		       const Double_t mapMax = 1.05,
		       const Bool_t iphiOnXaxis = true,
		       const Bool_t iphiOnXaxisSM = false
		       ) 
{

  system(Form("mkdir -p %s",outDir.c_str()));
  system(Form("cp /afs/cern.ch/user/m/mciprian/public/index.php %s",outDir.c_str()));


  TH1::SetDefaultSumw2();

  gStyle->SetPalette(55, 0);  // 55:raibow palette ; 57: kBird (blue to yellow) ; 107 kVisibleSpectrum ; 77 kDarkRainBow                                               
  gStyle->SetNumberContours(100); // default is 20 

  // ofstream badICxtalFile("xtal_with_IC_equal_0.txt");
  // if (not badICxtalFile.is_open()) {
  //   std::cout << "Error: could not open file xtal_with_IC_equal_0.txt" << std::endl;
  //   exit(EXIT_FAILURE);
  // }

  string xAxisName_foldSM = iphiOnXaxisSM ? "i#phi" : "i#eta";
  string yAxisName_foldSM = iphiOnXaxisSM ? "i#eta" : "i#phi";

  string xAxisName_allEB = iphiOnXaxis ? "i#phi" : "i#eta";
  string yAxisName_allEB = iphiOnXaxis ? "i#eta" : "i#phi";

  Int_t xsizeCanvas = 700;
  Int_t ysizeCanvas = 600;

  TFile* f = TFile::Open(inputFile.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error opening file \"" << inputFile << "\".\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }

  TH2 *mapEB = NULL;  // I expect it to have iphi on X and ieta on Y (IC map from DrawIC.py, not from TH2F in *_calibMap.root)

  mapEB = (TH2*) f->Get(mapNameInFile.c_str());
  if (!mapEB) {
    cout << "Error: could not get '" << mapNameInFile << "' histogram. End of programme" << endl;
    exit(EXIT_FAILURE);
  }
  mapEB->SetDirectory(0);
  f->Close();

  TH2F *mapEB_original = new TH2F("mapEB_original",Form("original EB map"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  TH2F *mapEB_norm1etaRing = new TH2F("mapEB_norm1etaRing",Form("original EB map normalized to 1 in #eta-ring"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  TH2F *mapEB_SM = new TH2F("mapEB_SM",Form("EB map folded in SM"), 85, 0.5, 85.5 , 20, 0.5, 20.5);
  TH2F *mapEB_SM_allEB = new TH2F("mapEB_SM_allEB",Form("EB map folded in SM, repeated in full EB"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  TH2F *mapEB_divided_mapEB_SM_allEB = new TH2F("mapEB_divided_mapEB_SM_allEB",Form("EB map divided by folded map"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  TH2F *mapEB_norm1_eachModule = new TH2F("mapEB_norm1_eachModule",Form("EB map normalized to 1 in each module"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  TH2F *mapEB_original_Over_norm1eachModuleSMallEB = new TH2F("mapEB_original_Over_norm1eachModuleSMallEB",Form("EB map divided by folded map (after norm. to 1 in each module)"), 360, 0.5, 360.5, 171, -85.5, 85.5);


  TH2F *mapEB_original_Over_norm1eachModuleSMallEB_renorm1etaRing = new TH2F("mapEB_original_Over_norm1eachModuleSMallEB_renorm1etaRing",Form("EB map divided by folded map (after norm. to 1 in each module), renorm.to 1 in #eta-ring"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  
  TH2F *mapEB_norm1_eachModule_SM = new TH2F("mapEB_norm1_eachModule_SM",Form("EB map norm. to 1 in each module, folded in SM"), 85, 0.5, 85.5 , 20, 0.5, 20.5);
  TH2F *mapEB_norm1_eachModule_SM_allEB = new TH2F("mapEB_norm1_eachModule_SM_allEB",Form("EB map norm. to 1 in each module folded in SM, repeated in full EB"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  TH2F *mapEB_norm1_eachModule_SM_allEB_plusMinusSeparate = new TH2F("mapEB_norm1_eachModule_SM_allEB_plusMinusSeparate",Form("EB map norm. to 1 in each module folded in SM (EB+ and EB- separately)"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  TH2F *mapEB_norm1_eachModule_SM_allEB_SamePhi = new TH2F("mapEB_norm1_eachModule_SM_allEB_SamePhi",Form("EB map norm. to 1 in each module folded in SM (EB+,EB- same #phi)"), 360, 0.5, 360.5, 171, -85.5, 85.5);

  TH2F *ratio_SM_over_norm1_eachModule_SM = new TH2F("ratio_SM_over_norm1_eachModule_SM",Form("ratio of folded map w/ and w/o norm. to 1 in each module"), 360, 0.5, 360.5, 171, -85.5, 85.5);

  TH2F *ratio_norm1_eachModule_foldSM_allEB_plusMinusSeparate = new TH2F("ratio_norm1_eachModule_foldSM_allEB_plusMinusSeparate",Form("ratio of folded map norm. to 1 in #eta-ring: all SM / EB+- separate"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  TH2F *ratio_norm1_eachModule_foldSM_allEB_SamePhi = new TH2F("ratio_norm1_eachModule_foldSM_allEB_SamePhi",Form("ratio of folded map norm. to 1 in #eta-ring: all SM / EB+- same #phi"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  TH2F *ratio_norm1_eachModule_foldSM_plusMinusSeparate_SamePhi = new TH2F("ratio_norm1_eachModule_foldSM_plusMinusSeparate_SamePhi",Form("ratio of folded map norm. to 1 in #eta-ring: EB+- separate / EB+- same #phi"), 360, 0.5, 360.5, 171, -85.5, 85.5);

  // fold not with EBDetId::ic() but putting same iphi between EB+ and EB-
  TH2F *mapEB_norm1_eachModule_SM_samePhi = new TH2F("mapEB_norm1_eachModule_SM_samePhi",Form("EB map norm. to 1 in each module, folded in SM (EB+ & EB- same #phi)"), 85, 0.5, 85.5 , 20, 0.5, 20.5);


  TH1F* distrIC_original = new TH1F("distrIC_original","original IC",200,0.9,1.1);
  TH1F* distrIC_norm1_eachModule_SM_allEB = new TH1F("distrIC_norm1_eachModule_SM_allEB","IC normalized to 1 in each module and folded in SM",200,0.9,1.1);
  TH1F* distrIC_correctedTT = new TH1F("distrIC_correctedTT","IC after TT correction",200,0.9,1.1);
  TH1F* distrIC_norm1_eachModule = new TH1F("distrIC_norm1_eachModule","IC normalized to 1 in each module",200,0.9,1.1);

  TProfile* profilePhi_mapEB_norm1_eachModule_SM = new TProfile("profilePhi_mapEB_norm1_eachModule_SM","#phi profile of EB map folded in SM (norm. to 1 in each module)",20,0.5,20.5);
  TProfile* profilePhi_mapEB_norm1_eachModule_SM_samePhi = new TProfile("profilePhi_mapEB_norm1_eachModule_SM_samePhi","#phi profile of EB map folded in SM (norm. to 1 in each module, EB+ & EB- same #phi)",20,0.5,20.5);

  copyMapAllEB(mapEB_original, mapEB, iphiOnXaxis);
  copyMapAllEB(mapEB_norm1etaRing, mapEB, iphiOnXaxis);
  copyMapAllEB(mapEB_norm1_eachModule, mapEB, iphiOnXaxis);

  normalizeEBMapTo1_inEachModule(mapEB_norm1_eachModule, iphiOnXaxis);
  normalizeEBMapTo1_inEtaRing(mapEB_norm1etaRing, iphiOnXaxis);

  foldEBMapInSM(mapEB, mapEB_SM, mapEB_SM_allEB, iphiOnXaxis, iphiOnXaxisSM, true);
  foldEBMapInSM(mapEB_norm1_eachModule, mapEB_norm1_eachModule_SM, mapEB_norm1_eachModule_SM_allEB, iphiOnXaxis, iphiOnXaxisSM, true);
  makeICprofileIphiFromMap(profilePhi_mapEB_norm1_eachModule_SM, mapEB_norm1_eachModule_SM, false, true, true);
  foldEBMapInSM(mapEB_norm1_eachModule, mapEB_norm1_eachModule_SM_samePhi, mapEB_norm1_eachModule_SM_allEB_SamePhi, iphiOnXaxis, iphiOnXaxisSM, false);
  makeICprofileIphiFromMap(profilePhi_mapEB_norm1_eachModule_SM_samePhi, mapEB_norm1_eachModule_SM_samePhi, false, true, true);

  foldEBMapInSM_plusMinusSeparate(mapEB_norm1_eachModule, mapEB_norm1_eachModule_SM_allEB_plusMinusSeparate, iphiOnXaxis, true);

  mapEB_divided_mapEB_SM_allEB->Divide(mapEB_original, mapEB_SM_allEB);
  mapEB_original_Over_norm1eachModuleSMallEB->Divide(mapEB_original, mapEB_norm1_eachModule_SM_allEB);
  ratio_SM_over_norm1_eachModule_SM->Divide(mapEB_SM_allEB, mapEB_norm1_eachModule_SM_allEB);

  ratio_norm1_eachModule_foldSM_allEB_plusMinusSeparate->Divide(mapEB_norm1_eachModule_SM_allEB, mapEB_norm1_eachModule_SM_allEB_plusMinusSeparate);
  ratio_norm1_eachModule_foldSM_allEB_SamePhi->Divide(mapEB_norm1_eachModule_SM_allEB, mapEB_norm1_eachModule_SM_allEB_SamePhi);
  ratio_norm1_eachModule_foldSM_plusMinusSeparate_SamePhi->Divide(mapEB_norm1_eachModule_SM_allEB_plusMinusSeparate, mapEB_norm1_eachModule_SM_allEB_SamePhi);

  copyMapAllEB(mapEB_original_Over_norm1eachModuleSMallEB_renorm1etaRing, mapEB_original_Over_norm1eachModuleSMallEB, iphiOnXaxis);
  normalizeEBMapTo1_inEtaRing(mapEB_original_Over_norm1eachModuleSMallEB_renorm1etaRing, iphiOnXaxis);

  makeICdistributionFromMap(distrIC_original, mapEB_original, true);
  makeICdistributionFromMap(distrIC_norm1_eachModule_SM_allEB, mapEB_norm1_eachModule_SM_allEB, true);
  makeICdistributionFromMap(distrIC_correctedTT, mapEB_original_Over_norm1eachModuleSMallEB, true);
  makeICdistributionFromMap(distrIC_norm1_eachModule, mapEB_norm1_eachModule, true);

  TH2F *mapEB_original_Over_norm1eachModuleSMallEB_renorm1etaRing_foldSM = new TH2F("mapEB_original_Over_norm1eachModuleSMallEB_renorm1etaRing_foldSM",Form("EB map divided by folded map (after norm. to 1 in each module), renorm.to 1 in #eta-ring, fold SM"), 85, 0.5, 85.5, 20, 0.5, 20.5);
  TH2F *mapEB_original_Over_norm1eachModuleSMallEB_renorm1etaRing_foldSM_allEB = new TH2F("mapEB_original_Over_norm1eachModuleSMallEB_renorm1etaRing_foldSM_allEB",Form("EB map divided by folded map (after norm. to 1 in each module), renorm.to 1 in #eta-ring, fold SM all EB"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  foldEBMapInSM(mapEB_original_Over_norm1eachModuleSMallEB_renorm1etaRing, mapEB_original_Over_norm1eachModuleSMallEB_renorm1etaRing_foldSM, mapEB_original_Over_norm1eachModuleSMallEB_renorm1etaRing_foldSM_allEB, iphiOnXaxis, iphiOnXaxisSM, true);

  /////////////////////
  // plots here

  xsizeCanvas = 600;
  ysizeCanvas = 2. * xsizeCanvas * mapEB_SM->GetNbinsY() / mapEB_SM->GetNbinsX() + 0.1 *xsizeCanvas;
  drawMap(mapEB_SM, xAxisName_foldSM, yAxisName_foldSM, mapNameInFile+"_foldSM", outDir, mapMin, mapMax, xsizeCanvas, ysizeCanvas);
  drawMap(mapEB_norm1_eachModule_SM, xAxisName_foldSM, yAxisName_foldSM, mapNameInFile+"_norm1_eachModule_foldSM", outDir, max(mapMin,0.975), min(mapMax,1.025), xsizeCanvas, ysizeCanvas);  
  drawMap(mapEB_norm1_eachModule_SM_samePhi, xAxisName_foldSM, yAxisName_foldSM, mapNameInFile+"_norm1_eachModule_foldSM_samePhi", outDir, max(mapMin,0.975), min(mapMax,1.025), xsizeCanvas, ysizeCanvas);
  drawMap(mapEB_original_Over_norm1eachModuleSMallEB_renorm1etaRing_foldSM, xAxisName_foldSM, yAxisName_foldSM, mapNameInFile+"_divided_foldSMAfterNorm1eachModule_renorm1etaRing_foldSM", outDir, max(mapMin,0.975), min(mapMax,1.025), xsizeCanvas, ysizeCanvas);

  // change size
  xsizeCanvas = 1200;
  ysizeCanvas = 1.0 * xsizeCanvas * 171. / 360. + 0.1 *xsizeCanvas;

  // for maps normalized, custom the range

  drawMap(mapEB_original, xAxisName_allEB, yAxisName_allEB, mapNameInFile+"_original" , outDir, mapMin, mapMax, xsizeCanvas, ysizeCanvas);
  drawMap(mapEB_norm1etaRing, xAxisName_allEB, yAxisName_allEB, mapNameInFile+"_norm1etaRing" , outDir, max(mapMin,0.95), min(mapMax,1.05), xsizeCanvas, ysizeCanvas);
  drawMap(mapEB_SM_allEB, xAxisName_allEB, yAxisName_allEB, mapNameInFile+"_foldSM_allEB" , outDir, mapMin, mapMax, xsizeCanvas, ysizeCanvas);
  drawMap(mapEB_norm1_eachModule_SM_allEB, xAxisName_allEB, yAxisName_allEB, mapNameInFile+"_norm1_eachModule_foldSM_allEB" , outDir, max(mapMin,0.95), min(mapMax,1.05), xsizeCanvas, ysizeCanvas);
  drawMap(mapEB_norm1_eachModule_SM_allEB_plusMinusSeparate, xAxisName_allEB, yAxisName_allEB, mapNameInFile+"_norm1_eachModule_foldSM_allEB_plusMinusSeparate" , outDir, max(mapMin,0.95), min(mapMax,1.05), xsizeCanvas, ysizeCanvas);
  drawMap(mapEB_norm1_eachModule_SM_allEB_SamePhi, xAxisName_allEB, yAxisName_allEB, mapNameInFile+"_norm1_eachModule_foldSM_allEB_SamePhi" , outDir, max(mapMin,0.95), min(mapMax,1.05), xsizeCanvas, ysizeCanvas);
  drawMap(mapEB_divided_mapEB_SM_allEB, xAxisName_allEB, yAxisName_allEB, mapNameInFile+"_divided" , outDir, max(mapMin,0.95), min(mapMax,1.05), xsizeCanvas, ysizeCanvas);
  drawMap(mapEB_original_Over_norm1eachModuleSMallEB, xAxisName_allEB, yAxisName_allEB, mapNameInFile+"_divided_foldSMAfterNorm1eachModule" , outDir, mapMin, mapMax, xsizeCanvas, ysizeCanvas);
  drawMap(mapEB_norm1_eachModule, xAxisName_allEB, yAxisName_allEB, mapNameInFile+"_norm1_eachModule" , outDir, max(mapMin,0.95), min(mapMax,1.05), xsizeCanvas, ysizeCanvas);
  drawMap(ratio_SM_over_norm1_eachModule_SM, xAxisName_allEB, yAxisName_allEB, mapNameInFile+"_ratio_foldWithWithoutNorm1eachModule" , outDir, 0.995, 1.005, xsizeCanvas, ysizeCanvas);
  drawMap(ratio_norm1_eachModule_foldSM_allEB_plusMinusSeparate, xAxisName_allEB, yAxisName_allEB, mapNameInFile+"_ratio_foldSM_allEB_over_plusMinusSeparate" , outDir, 0.975, 1.025, xsizeCanvas, ysizeCanvas);  
  drawMap(ratio_norm1_eachModule_foldSM_allEB_SamePhi, xAxisName_allEB, yAxisName_allEB, mapNameInFile+"_ratio_foldSM_allEB_over_SamePhi" , outDir, 0.975, 1.025, xsizeCanvas, ysizeCanvas);
  drawMap(ratio_norm1_eachModule_foldSM_plusMinusSeparate_SamePhi, xAxisName_allEB, yAxisName_allEB, mapNameInFile+"_ratio_foldSM_plusMinusSeparate_over_SamePhi" , outDir, 0.975, 1.025, xsizeCanvas, ysizeCanvas);

  drawMap(mapEB_original_Over_norm1eachModuleSMallEB_renorm1etaRing, xAxisName_allEB, yAxisName_allEB, mapNameInFile+"_divided_foldSMAfterNorm1eachModule_renorm1etaRing" , outDir, max(mapMin,0.95), min(mapMax,1.05), xsizeCanvas, ysizeCanvas);
  drawMap(mapEB_original_Over_norm1eachModuleSMallEB_renorm1etaRing_foldSM_allEB, xAxisName_allEB, yAxisName_allEB, mapNameInFile+"_divided_foldSMAfterNorm1eachModule_renorm1etaRing_foldSM_allEB" , outDir, max(mapMin,0.975), min(mapMax,1.025), xsizeCanvas, ysizeCanvas);

  drawDistribution(distrIC_original, "intercalibration constants", "Events", mapNameInFile+"_original_1D", outDir, mapMin, mapMax); 
  drawDistribution(distrIC_norm1_eachModule_SM_allEB, "intercalibration constants", "Events", mapNameInFile+"_norm1_eachModule_foldSM_allEB_1D", outDir, mapMin, mapMax); 
  drawDistribution(distrIC_correctedTT, "intercalibration constant", "Events", mapNameInFile+"_correctedTT_1D", outDir, mapMin, mapMax); 
  drawDistribution(distrIC_norm1_eachModule, "intercalibration constant", "Events", mapNameInFile+"_norm1_eachModule_1D", outDir, mapMin, mapMax); 

  drawDistribution(profilePhi_mapEB_norm1_eachModule_SM, "i#phi", "mean IC", mapNameInFile+"_norm1_eachModule_foldSM_profileIphi", outDir, 0.5, 20.5); 
  drawDistribution(profilePhi_mapEB_norm1_eachModule_SM_samePhi, "i#phi", "mean IC", mapNameInFile+"_norm1_eachModule_foldSM_samePhi_profileIphi", outDir, 0.5, 20.5); 

  // delete mapEB_SM;
  // delete mapEB_SM_allEB;
  // delete mapEB_original;
  // delete mapEB_norm1etaRing;
  // delete mapEB_divided_mapEB_SM_allEB;
  // delete mapEB_norm1_eachModule;
  // delete mapEB_norm1_eachModule_SM;
  // delete mapEB_norm1_eachModule_SM_allEB;
  // delete mapEB_original_Over_norm1eachModuleSMallEB;
  // delete mapEB_original_Over_norm1eachModuleSMallEB_renorm1etaRing;
  // delete mapEB_norm1_eachModule_SM_samePhi;
  // delete mapEB_norm1_eachModule_SM_allEB_plusMinusSeparate;
  // delete mapEB_norm1_eachModule_SM_allEB_SamePhi;
  // delete ratio_norm1_eachModule_foldSM_allEB_plusMinusSeparate;
  // delete ratio_norm1_eachModule_foldSM_allEB_SamePhi;

  // delete distrIC_original;
  // delete distrIC_norm1_eachModule_SM_allEB;
  // delete distrIC_correctedTT;
  // delete distrIC_norm1_eachModule;

  // delete profilePhi_mapEB_norm1_eachModule_SM;
  // delete profilePhi_mapEB_norm1_eachModule_SM_samePhi;

  //badICxtalFile.close();

}
