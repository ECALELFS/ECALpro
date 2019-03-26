#include "./ICfunctions.h"

using namespace std;

void makeICmapAndDump(//const string& outDir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot_Legacy/AlCaP0_AllRun2017_condor/iter_1/2DMaps/ICmaps/IC_work/",
		      //const string& inputFile = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero_Run2/mciprian/AlCaP0_AllRun2017_condor/iter_1/AlCaP0_AllRun2017_condor_calibMap.root",
		      const string& outDir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/ICplot/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_ext1_fromIter6/iter_6/2DMaps/ICmaps/IC_work/",
		      const string& inputFile = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_ext1_fromIter6/iter_6/AlCaP0_Run2017_DE_run304366_ContCorrEoverEtrueScaledToV2MC_ext1_fromIter6_calibMap.root",
		      const string& outICdumpFileName = "dumpIC_norm1etaRing.dat",
		      const string& canvasNamePrefix = "calibMap_EB",
		      const Double_t mapMin = 0.95,
		      const Double_t mapMax = 1.05,
		      const Bool_t iphiOnXaxis = true,
		      const Bool_t iphiOnXaxisSM = false,
		      const Int_t all0_EB1_EE2 = 0,
		      const Bool_t excludeMod2EBm16 = false
		      ) 

{

  system(Form("mkdir -p %s",outDir.c_str()));
  system(Form("cp /afs/cern.ch/user/m/mciprian/public/index.php %s",outDir.c_str()));


  TH1::SetDefaultSumw2();

  gStyle->SetPalette(55, 0);  // 55:raibow palette ; 57: kBird (blue to yellow) ; 107 kVisibleSpectrum ; 77 kDarkRainBow                          
  gStyle->SetNumberContours(101); // default is 20 

  ofstream outICdumpFile((outDir+outICdumpFileName).c_str());
  if (not outICdumpFile.is_open()) {
    std::cout << "Error: could not open file " << outICdumpFileName << std::endl;
    exit(EXIT_FAILURE);
  }

  string xAxisName_foldSM = iphiOnXaxisSM ? "i#phi" : "i#eta";
  string yAxisName_foldSM = iphiOnXaxisSM ? "i#eta" : "i#phi";

  string xAxisName_allEB = iphiOnXaxis ? "i#phi" : "i#eta";
  string yAxisName_allEB = iphiOnXaxis ? "i#eta" : "i#phi";

  Int_t xsizeCanvas = 700;
  Int_t ysizeCanvas = 600;
  // change size
  xsizeCanvas = 1200;
  ysizeCanvas = 1.0 * xsizeCanvas * 171. / 360. + 0.1 *xsizeCanvas;

  Int_t xsizeCanvas_SM = 700;
  Int_t ysizeCanvas_SM = 600;
  // change size
  xsizeCanvas_SM = 1200;
  ysizeCanvas_SM = 2.0 * xsizeCanvas * 20. / 85. + 0.1 *xsizeCanvas;


  TH2F* mapEB_IC = new TH2F("mapEB_IC","EB map not normalized",360, 0.5, 360.5, 171, -85.5, 85.5);
  TH2F* mapEB_IC_err = new TH2F("mapEB_IC_err","EB uncertainty",360, 0.5, 360.5, 171, -85.5, 85.5);


  if (all0_EB1_EE2 != 2) {
    runOnTree(mapEB_IC, mapEB_IC_err, inputFile, 0, true);

    drawMap(mapEB_IC, xAxisName_allEB, yAxisName_allEB, "IC_CalibMapEB" , outDir, 0.9, 1.1, xsizeCanvas, ysizeCanvas);
    drawMap(mapEB_IC_err, xAxisName_allEB, yAxisName_allEB, "IC_err_CalibMapEB" , outDir, 0.0, 0.005, xsizeCanvas, ysizeCanvas, 50);
  }

  TH2F* mapEEp_IC = new TH2F("mapEEp_IC","EE+ map not normalized", 100, 0.5, 100.5, 100, 0.5, 100.5);
  TH2F* mapEEp_IC_err = new TH2F("mapEEp_IC_err","EE+ uncertainty", 100, 0.5, 100.5, 100, 0.5, 100.5);
  TH2F* mapEEm_IC = new TH2F("mapEEm_IC","EE- map not normalized", 100, 0.5, 100.5, 100, 0.5, 100.5);
  TH2F* mapEEm_IC_err = new TH2F("mapEEm_IC_err","EE- uncertainty", 100, 0.5, 100.5, 100, 0.5, 100.5);

  TH2F* mapEEp_IC_norm1etaRing = new TH2F("mapEEp_IC_norm1etaRing","EE+ map normalized to 1 in #eta-ring", 100, 0.5, 100.5, 100, 0.5, 100.5);
  TH2F* mapEEm_IC_norm1etaRing = new TH2F("mapEEm_IC_norm1etaRing","EE- map normalized to 1 in #eta-ring", 100, 0.5, 100.5, 100, 0.5, 100.5);

  if (all0_EB1_EE2 != 1) {
    runOnTree(mapEEp_IC, mapEEp_IC_err, inputFile, +1, true);
    runOnTree(mapEEm_IC, mapEEm_IC_err, inputFile, -1, true);

    copyMapEE(mapEEp_IC_norm1etaRing, mapEEp_IC);
    copyMapEE(mapEEm_IC_norm1etaRing, mapEEm_IC);
    normalizeEEMapTo1_inEtaRing(mapEEp_IC_norm1etaRing);
    normalizeEEMapTo1_inEtaRing(mapEEm_IC_norm1etaRing);

    drawMap(mapEEp_IC, "iX", "iY", "IC_CalibMapEEp" , outDir, 0.75, 1.25, 700, 700);
    drawMap(mapEEp_IC_err, "iX", "iY", "IC_err_CalibMapEEp" , outDir, 0.0, 0.02, 700, 700, 51);

    drawMap(mapEEm_IC, "iX", "iY", "IC_CalibMapEEm" , outDir, 0.75, 1.25, 700, 700);
    drawMap(mapEEm_IC_err, "iX", "iY", "IC_err_CalibMapEEm" , outDir, 0.0, 0.02, 700, 700, 51);

    drawMap(mapEEp_IC_norm1etaRing, "iX", "iY", "IC_CalibMapEEp_norm1etaRing" , outDir, 0.9, 1.1, 700, 700);
    drawMap(mapEEm_IC_norm1etaRing, "iX", "iY", "IC_CalibMapEEm_norm1etaRing" , outDir, 0.9, 1.1, 700, 700);

    //mapEEp_IC_norm1etaRing->SaveAs((outDir+"calibMap_EEp_norm1etaRing.root").c_str());
    //mapEEm_IC_norm1etaRing->SaveAs((outDir+"calibMap_EEm_norm1etaRing.root").c_str());
  }
  //////////////////////////


  // normalize original map to 1 in eta-ring (TT pattern still here)
  TH2F *mapEB_norm1etaRing = new TH2F("mapEB_norm1etaRing",Form("EB map normalized to 1 in #eta-ring"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  if (all0_EB1_EE2 != 2) {
    copyMapAllEB(mapEB_norm1etaRing, mapEB_IC, iphiOnXaxis);  
    normalizeEBMapTo1_inEtaRing(mapEB_norm1etaRing, iphiOnXaxis, excludeMod2EBm16); 
    drawMap(mapEB_norm1etaRing, xAxisName_allEB, yAxisName_allEB, canvasNamePrefix+"_norm1etaRing" , outDir, mapMin, mapMax, xsizeCanvas, ysizeCanvas);
  }
  // normalize original map to 1 in each module (TT pattern still here), prepare for folding
  TH2F *mapEB_norm1eachModule = new TH2F("mapEB_norm1eachModule",Form("EB map normalized to 1 in each module"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  if (all0_EB1_EE2 != 2) {
    copyMapAllEB(mapEB_norm1eachModule, mapEB_IC, iphiOnXaxis);  
    normalizeEBMapTo1_inEachModule(mapEB_norm1eachModule, iphiOnXaxis, excludeMod2EBm16);
    drawMap(mapEB_norm1eachModule, xAxisName_allEB, yAxisName_allEB, canvasNamePrefix+"_norm1eachModule" , outDir, mapMin, mapMax, xsizeCanvas, ysizeCanvas);
  }
  // fold map normalized to 1 in each module
  TH2F *mapEB_norm1eachModule_foldSM = new TH2F("mapEB_norm1eachModule_foldSM",Form("EB map norm. to 1 in each module, folded in SM"), 85, 0.5, 85.5 , 20, 0.5, 20.5);
  TH2F *mapEB_norm1eachModule_foldSM_allEB = new TH2F("mapEB_norm1eachModule_foldSM_allEB",Form("EB map norm. to 1 in each module folded in SM, repeated in full EB"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  TH2F *mapEB_norm1eachModule_foldSM_samePhi = new TH2F("mapEB_norm1eachModule_foldSM_samePhi",Form("EB map norm. to 1 in each module, folded in SM (EB+,EB- same #phi)"), 85, 0.5, 85.5 , 20, 0.5, 20.5);
  TH2F *mapEB_norm1eachModule_foldSM_allEB_samePhi = new TH2F("mapEB_norm1eachModule_foldSM_allEB_samePhi",Form("EB map norm. to 1 in each module folded in SM (EB+,EB- same #phi)"), 360, 0.5, 360.5, 171, -85.5, 85.5);

  TH2F *mapEB_norm1eachModule_foldSM_allEB_plusMinusSeparate = new TH2F("mapEB_norm1eachModule_foldSM_allEB_plusMinusSeparate",Form("EB map norm. to 1 in each module folded in SM, EB+,EB- separately"), 360, 0.5, 360.5, 171, -85.5, 85.5);

  if (all0_EB1_EE2 != 2) {
    foldEBMapInSM(mapEB_norm1eachModule, mapEB_norm1eachModule_foldSM, mapEB_norm1eachModule_foldSM_allEB, iphiOnXaxis, iphiOnXaxisSM);
    foldEBMapInSM(mapEB_norm1eachModule, mapEB_norm1eachModule_foldSM_samePhi, mapEB_norm1eachModule_foldSM_allEB_samePhi, iphiOnXaxis, iphiOnXaxisSM, false);
    foldEBMapInSM_plusMinusSeparate(mapEB_norm1eachModule, mapEB_norm1eachModule_foldSM_allEB_plusMinusSeparate, iphiOnXaxis);

    drawMap(mapEB_norm1eachModule_foldSM, xAxisName_foldSM, yAxisName_foldSM, canvasNamePrefix+"_norm1eachModule_foldSM", outDir, mapMin, mapMax, xsizeCanvas_SM, ysizeCanvas_SM);
    drawMap(mapEB_norm1eachModule_foldSM_samePhi, xAxisName_foldSM, yAxisName_foldSM, canvasNamePrefix+"_norm1eachModule_foldSM_samePhi", outDir, mapMin, mapMax, xsizeCanvas_SM, ysizeCanvas_SM);
    
    drawMap(mapEB_norm1eachModule_foldSM_allEB, xAxisName_allEB, yAxisName_allEB, canvasNamePrefix+"_norm1eachModule_foldSM_allEB" , outDir, mapMin, mapMax, xsizeCanvas, ysizeCanvas);
    drawMap(mapEB_norm1eachModule_foldSM_allEB_samePhi, xAxisName_allEB, yAxisName_allEB, canvasNamePrefix+"_norm1eachModule_foldSM_allEB_samePhi" , outDir, mapMin, mapMax, xsizeCanvas, ysizeCanvas);
    drawMap(mapEB_norm1eachModule_foldSM_allEB_plusMinusSeparate, xAxisName_allEB, yAxisName_allEB, canvasNamePrefix+"_norm1eachModule_foldSM_allEB_plusMinusSeparate" , outDir, mapMin, mapMax, xsizeCanvas, ysizeCanvas);
  }

  TH2F *mapEB_original_Over_norm1eachModuleFoldSMallEB_norm1etaRing = new TH2F("mapEB_original_Over_norm1eachModuleFoldSMallEB_norm1etaRing",Form("EB map divided by folded map (norm. to 1 in each module), norm. to 1 in #eta-ring"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  if (all0_EB1_EE2 != 2) {
    divideEBmap(mapEB_original_Over_norm1eachModuleFoldSMallEB_norm1etaRing, mapEB_IC, mapEB_norm1eachModule_foldSM_allEB, true, -1.0); 
    normalizeEBMapTo1_inEtaRing(mapEB_original_Over_norm1eachModuleFoldSMallEB_norm1etaRing, iphiOnXaxis, excludeMod2EBm16);
    drawMap(mapEB_original_Over_norm1eachModuleFoldSMallEB_norm1etaRing, xAxisName_allEB, yAxisName_allEB, canvasNamePrefix+"_divided_foldSMafterNorm1eachModule_norm1etaRing" , outDir, mapMin, mapMax, xsizeCanvas, ysizeCanvas);
  }

  TH2F *mapEB_original_Over_norm1eachModuleFoldSMallEB_samePhi_norm1etaRing = new TH2F("mapEB_original_Over_norm1eachModuleFoldSMallEB_samePhi_norm1etaRing",Form("EB map divided by folded map (norm. to 1 in each module, EB+,EB- same #phi), norm. to 1 in #eta-ring"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  if (all0_EB1_EE2 != 2) {
    divideEBmap(mapEB_original_Over_norm1eachModuleFoldSMallEB_samePhi_norm1etaRing, mapEB_IC, mapEB_norm1eachModule_foldSM_allEB_samePhi, true, -1.0);
    normalizeEBMapTo1_inEtaRing(mapEB_original_Over_norm1eachModuleFoldSMallEB_samePhi_norm1etaRing, iphiOnXaxis, excludeMod2EBm16);
    drawMap(mapEB_original_Over_norm1eachModuleFoldSMallEB_samePhi_norm1etaRing, xAxisName_allEB, yAxisName_allEB, canvasNamePrefix+"_divided_foldSMafterNorm1eachModuleSamePhi_norm1etaRing" , outDir, mapMin, mapMax, xsizeCanvas, ysizeCanvas);
  }

  TH2F *mapEB_original_Over_norm1eachModuleFoldSMallEB_plusMinusSeparate_norm1etaRing = new TH2F("mapEB_original_Over_norm1eachModuleFoldSMallEB_plusMinusSeparate_norm1etaRing",Form("EB map divided by folded map (norm. to 1 in each module, EB+,EB- separately), norm. to 1 in #eta-ring"), 360, 0.5, 360.5, 171, -85.5, 85.5);
  if (all0_EB1_EE2 != 2) {
    divideEBmap(mapEB_original_Over_norm1eachModuleFoldSMallEB_plusMinusSeparate_norm1etaRing, mapEB_IC, mapEB_norm1eachModule_foldSM_allEB_plusMinusSeparate, true, -1.0);
    normalizeEBMapTo1_inEtaRing(mapEB_original_Over_norm1eachModuleFoldSMallEB_plusMinusSeparate_norm1etaRing, iphiOnXaxis, excludeMod2EBm16);
    drawMap(mapEB_original_Over_norm1eachModuleFoldSMallEB_plusMinusSeparate_norm1etaRing, xAxisName_allEB, yAxisName_allEB, canvasNamePrefix+"_divided_foldSMafterNorm1eachModulePlusMinusSeparate_norm1etaRing" , outDir, mapMin, mapMax, xsizeCanvas, ysizeCanvas);
    //mapEB_original_Over_norm1eachModuleFoldSMallEB_plusMinusSeparate_norm1etaRing->SaveAs((outDir+canvasNamePrefix+"_divided_foldSMafterNorm1eachModulePlusMinusSeparate_norm1etaRing.root").c_str());
    checkICnormalizedTo1_inEtaRing(mapEB_original_Over_norm1eachModuleFoldSMallEB_plusMinusSeparate_norm1etaRing, true, true);
  }


  TProfile* profileEB_phi_final = new TProfile("profileEB_phi_final","i#phi profile of IC in EB (norm. to 1 in #eta-ring, TT corrected folding EB+/EB- separately)",360, 0.5, 360.5);
  makeICprofileIphiFromMap(profileEB_phi_final, mapEB_original_Over_norm1eachModuleFoldSMallEB_plusMinusSeparate_norm1etaRing, true, true, false);
  drawDistribution(profileEB_phi_final, "i#phi", "mean IC", canvasNamePrefix+"_divided_foldSMafterNorm1eachModulePlusMinusSeparate_norm1etaRing_iphiProfile", outDir, 0.5, 360.5, 700, 500, true);

  Double_t IC_value = 0.0;
  Double_t IC_error = 0.0; // stat. only for now
  Int_t bin = 0;
  TH2* map_IC_value = mapEB_original_Over_norm1eachModuleFoldSMallEB_plusMinusSeparate_norm1etaRing;
  TH2* map_IC_error = mapEB_IC_err;
  Int_t iz = 0; // 0 for EB

  outICdumpFile << "#--------------------------------------" << endl;   
  outICdumpFile << "# ieta(ix)  iphi(iy)  iz  IC  stat.unc." << endl;   
  outICdumpFile << "#--------------------------------------" << endl;   

  if (all0_EB1_EE2 != 2) {

    for (Int_t iphi = 1; iphi <= 360; ++iphi) {

      for (Int_t ieta = -85; ieta <= 85; ++ieta) {
      
	if (ieta == 0) continue;

	if (iphiOnXaxis) bin = map_IC_value->FindFixBin(iphi,ieta);
	else             bin = map_IC_value->FindFixBin(ieta,iphi);
	IC_value = map_IC_value->GetBinContent(bin); 
	if (IC_value < 0) IC_value = 1;
	IC_error = (map_IC_error->GetBinContent(bin) > EPSILON) ? map_IC_error->GetBinContent(bin) : 999; 

	outICdumpFile << right << setw(5) << ieta << " " 
		      << right << setw(5) << iphi << " "
		      << right << setw(3) << iz << "   "
		      << left << setw(9) << IC_value << " "
		      << left << setw(9) << IC_error << " "
		      << endl;
      
      }

    }

  }

  if (all0_EB1_EE2 != 1) {

    map_IC_value = mapEEp_IC_norm1etaRing;
    map_IC_error = mapEEp_IC_err;
    iz = 1; // 0 for EE+

    for (Int_t ix = 1; ix <= 100; ++ix) {

      for (Int_t iy = 1; iy <= 100; ++iy) {

	if (not EEDetId::validDetId(ix, iy, iz)) {
	  continue;
	  // IC_value = -1;
	  // IC_error = 999;
	} else {      
	  bin = map_IC_value->FindFixBin(ix,iy);
	  IC_value = map_IC_value->GetBinContent(bin); 
	  if (IC_value < 0) IC_value = 1;
	  IC_error = (map_IC_error->GetBinContent(bin) > EPSILON) ? map_IC_error->GetBinContent(bin) : 999; 
	}

	outICdumpFile << right << setw(5) << ix << " " 
		      << right << setw(5) << iy << " "
		      << right << setw(3) << iz << "   "
		      << left << setw(9) << IC_value << " "
		      << left << setw(9) << IC_error << " "
		      << endl;
      
      }

    }

    map_IC_value = mapEEm_IC_norm1etaRing;
    map_IC_error = mapEEm_IC_err;
    iz = -1; // 0 for EE-

    for (Int_t ix = 1; ix <= 100; ++ix) {

      for (Int_t iy = 1; iy <= 100; ++iy) {
      
	if (not EEDetId::validDetId(ix, iy, iz)) {
	  continue;
	  // IC_value = -1;
	  // IC_error = 999;
	} else {
	  bin = map_IC_value->FindFixBin(ix,iy);
	  IC_value = map_IC_value->GetBinContent(bin); 
	  if (IC_value < 0) IC_value = 1;
	  IC_error = (map_IC_error->GetBinContent(bin) > EPSILON) ? map_IC_error->GetBinContent(bin) : 999;
	}	

	outICdumpFile << right << setw(5) << ix << " " 
		      << right << setw(5) << iy << " "
		      << right << setw(3) << iz << "   "
		      << left << setw(9) << IC_value << " "
		      << left << setw(9) << IC_error << " "
		      << endl;
      
      }

    }

  }


  // no real need to delete, at the end of the script they will be automatically deleted
  // delete mapEB_IC;
  // delete mapEB_IC_err;
  // delete mapEB_norm1etaRing;
  // delete mapEB_norm1eachModule;
  // delete mapEB_norm1eachModule_foldSM;
  // delete mapEB_norm1eachModule_foldSM_allEB;
  // delete mapEB_original_Over_norm1eachModuleFoldSMallEB_norm1etaRing;

  // delete mapEEp_IC;
  // delete mapEEp_IC_err;
  // delete mapEEm_IC;
  // delete mapEEm_IC_err;


  outICdumpFile.close();

  // now saves maps for later usage:
  cout << "Saving maps in file: " << outDir << "calibrationMaps.root" << endl;
  string fullname = outDir + "calibrationMaps.root";
  TFile* f = TFile::Open(fullname.c_str(),"recreate");
  if (!f || !f->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error opening file \"" << inputFile << "\".\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }
  mapEEp_IC_norm1etaRing->Write("calibMap_EEp");  
  mapEEm_IC_norm1etaRing->Write("calibMap_EEm");  
  mapEB_original_Over_norm1eachModuleFoldSMallEB_plusMinusSeparate_norm1etaRing->Write("calibMap_EB");
  f->Close();
  delete f;

}
