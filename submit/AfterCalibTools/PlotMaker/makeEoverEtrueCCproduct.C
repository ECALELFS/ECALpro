#include "./ICfunctions.h"

using namespace std;

// macro to plot product of two containment corrections maps in EB 
// the usefulness is that you can use a CC map to start producing a new CC map and then 
// make the product of the CC to mimic an iterative CC derivation (which is difficult to implement because the CC are for each photon
// and when you create the clusters, which is the moment when the constants of previous iteration would be applied, you don't know if
// you a have a photon 1 or 2 yet
// pass:
// 1) the name of the folder where to store the output plots
// 2) the input file name where the maps are stored (the file is probably on EOS, use root://eoscms//eos/cms/...)
 
// the map in the calibMap file has ieta on x axis, but we produce the map with iphi on x axis (for full EB)
// we keep ieta on x axis for a single SM

void realDrawMapProduct(const string& outDir = "",
			const string& inputFile = "",
			const string& inputFile2 = "",
			const Int_t nPhoton = 1,
			const Double_t mapMin = 1.0,
			const Double_t mapMax = 1.12
		      )   
{

  TH1::SetDefaultSumw2();

  gStyle->SetPalette(55, 0);  // 55:raibow palette ; 57: kBird (blue to yellow) ; 107 kVisibleSpectrum ; 77 kDarkRainBow          
  gStyle->SetNumberContours(100); // default is 20 

  TFile* f = TFile::Open(inputFile.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error opening file \"" << inputFile << "\".\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }

  TH2F *mapEB = NULL;
  mapEB = (TH2F*) f->Get((nPhoton == 1) ? "calibMap_EB" : "calibMap_EB_g2");
  if (!mapEB) {
    cout << "Error: could not get EB histogram. End of programme" << endl;
    exit(EXIT_FAILURE);
  }
  mapEB->SetDirectory(0);
  f->Close();


  TFile* f2 = TFile::Open(inputFile2.c_str(),"READ");
  if (!f2 || !f2->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error opening file \"" << inputFile2 << "\".\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }

  TH2F *mapEB2 = NULL;
  mapEB2 = (TH2F*) f2->Get((nPhoton == 1) ? "calibMap_EB" : "calibMap_EB_g2");
  if (!mapEB2) {
    cout << "Error: could not get EB histogram 2. End of programme" << endl;
    exit(EXIT_FAILURE);
  }
  mapEB2->SetDirectory(0);
  f2->Close();

  Double_t mapEB_binContent = 0.0;
  Int_t bin_in = 0;
  Int_t bin_out = 0;

  TH2F* mapEB_new = new TH2F("mapEB_new","", 360, 0.5, 360.5, 171, -85.5, 85.5);
  TH2F* mapEB2_new = new TH2F("mapEB2_new","", 360, 0.5, 360.5, 171, -85.5, 85.5);

  for (Int_t iphi = 1; iphi <= 360; iphi++) {

    for (Int_t ieta = -85; ieta <= 85; ieta++) {

      bin_in = mapEB->FindFixBin(ieta,iphi);
      bin_out = mapEB_new->FindFixBin(iphi,ieta);
      mapEB_binContent = mapEB->GetBinContent( bin_in );
      mapEB_new->SetBinContent(bin_out, mapEB_binContent);
      mapEB_binContent = mapEB2->GetBinContent( bin_in );
      mapEB2_new->SetBinContent(bin_out, mapEB_binContent);

    }

  }

  TH2F* hProduct = new TH2F("hProduct","CC product", 360, 0.5, 360.5, 171, -85.5, 85.5);
  copyMapAllEB(hProduct,mapEB_new,true);
  hProduct->Multiply(mapEB2_new);

  // this map is useful if the one in all EB is obtained folding all SM, so we plot one single SM to make it better visible
  TH2F* hProduct_SM = new TH2F("hProduct_SM","CC product, single SM", 85, 0.5, 85.5, 20, 0.5, 20.5);

  for (Int_t iphi = 1; iphi <= 20; ++iphi) {

    for (Int_t ieta = 1; ieta <= 85; ++ieta) {
      
      bin_in = hProduct->FindFixBin(iphi,ieta);
      bin_out = hProduct_SM->FindFixBin(ieta, iphi);
      mapEB_binContent = hProduct->GetBinContent(bin_in);
      hProduct_SM->SetBinContent(bin_out,mapEB_binContent);

    }
    
  }

  //EB
  Int_t xsizeCanvas = 1200;
  Int_t ysizeCanvas = 1.0 * xsizeCanvas * mapEB_new->GetNbinsY() / mapEB_new->GetNbinsX() + 0.1 *xsizeCanvas;

  // string xAxisName = iphiOnXaxis ? "i#phi" : "i#eta";
  // string yAxisName = iphiOnXaxis ? "i#eta" : "i#phi";
  drawMap(hProduct, "i#phi", "i#eta", Form("calibMap_EB_g%d_product",nPhoton), outDir, mapMin, mapMax, xsizeCanvas, ysizeCanvas);

  xsizeCanvas = 1200;
  ysizeCanvas = 1.0 * xsizeCanvas * 171. / 360. + 0.1 *xsizeCanvas;
  drawMap(hProduct_SM, "i#eta", "i#phi", Form("calibMap_EB_g%d_product_singleSM",nPhoton), outDir, mapMin, mapMax, xsizeCanvas, ysizeCanvas);

  TProfile* hProduct_SM_iphiProfile = new TProfile("hProduct_SM_iphiProfile","Containment corrections in SM: i#phi profile",20, 0.5, 20.5);
  makeICprofileIphiFromMap(hProduct_SM_iphiProfile, hProduct_SM, false, false, true);
  drawDistribution(hProduct_SM_iphiProfile, "i#phi", "mean CC", Form("calibMap_EB_g%d_product_singleSM_iphiProfile",nPhoton), outDir, 0.5, 20.5, 700, 600, true);

  string outFileOpenMode = (nPhoton == 1) ? "RECREATE" : "UPDATE";
  TFile* fout = TFile::Open(Form("%s/ContainmentCorrections_EoverEtrue.root",outDir.c_str()), outFileOpenMode.c_str());
  if (!fout || !fout->IsOpen()) {
    cout << "*******************************" << endl;
    cout << "Error when opening output file.\nApplication will be terminated." << endl;
    cout << "*******************************" << endl;
    exit(EXIT_FAILURE);
  }

  TH2F* hProduct_ietaOnX = (TH2F*) mapEB->Clone();  // save object with same convention of maps in CC file, i.e., ieta on X axis
  hProduct_ietaOnX->Multiply(mapEB2);
  hProduct_ietaOnX->Write((nPhoton == 1) ? "calibMap_EB" : "calibMap_EB_g2");
  fout->Close();


  delete mapEB_new;
  delete mapEB2_new;
  delete hProduct;
  delete hProduct_SM;
  delete hProduct_SM_iphiProfile;

}


void makeEoverEtrueCCproduct(const string& outDir = "/afs/cern.ch/user/m/mciprian/www/pi0calib/CC_EoverEtrue/product_CC/pi0Gun_MC_EoverEtrue_foldSM_v4_iter2/",
			     const string& inputFile1 = "root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero2017/mciprian/pi0Gun_MCV1_EoverEtrue_foldSM_EoverEtrueCC_iter2_notScaled/iter_0/pi0Gun_MCV1_EoverEtrue_foldSM_EoverEtrueCC_iter2_notScaled_calibMap.root",
			     const string& inputFile2 = "/afs/cern.ch/user/m/mciprian/www/pi0calib/CC_EoverEtrue/product_CC/pi0Gun_MC_EoverEtrue_foldSM_v4_iter1/ContainmentCorrections_EoverEtrue.root",
			     const Double_t mapMin = 1.0,
			     const Double_t mapMax = 1.12 
			     ) 
{

  // maps in calibMap.root file have ieta on X axis, but it could change
  // in this macro we draw the map with iphi on X axis if full EB, and opposite if 1 SM

  system(Form("mkdir -p %s",outDir.c_str()));
  system(Form("cp /afs/cern.ch/user/m/mciprian/public/index.php %s",outDir.c_str()));

  realDrawMapProduct(outDir, inputFile1, inputFile2, 1, mapMin, mapMax);
  realDrawMapProduct(outDir, inputFile1, inputFile2, 2, mapMin, mapMax);
  

}
