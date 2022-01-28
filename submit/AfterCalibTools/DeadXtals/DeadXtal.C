#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <cmath>
#include <cstdlib>
#include <sstream>
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

//Usage: .x DeadXtal.C+("root://eoscms//eos/cms/store/group/alca_ecalcalib/lpernie/ALL_2012C_pi0_NewTag_01/iter_0/2012C_epsilonPlots.root","plots/")
void DeadXtal( TString Path, TString OutDir ){

  //Starting
  cout<<"Starting with: DeadXtal.C"<<endl;
  TFile* fin = TFile::Open( Path.Data() );
  if(!fin) { cout << "Invalid file: " << Path.Data() << " .. try again" << endl; return; }
  std::string calibMapFileName = string(Path.Data());
  std::string strToReplace = "epsilonPlots";
  calibMapFileName.replace(calibMapFileName.find(strToReplace.c_str()),strToReplace.size(),"calibMap");

  TFile* MapFile = TFile::Open( calibMapFileName.c_str() );
  TTree *calibEB = (TTree*) MapFile->Get("calibEB");
  TTree *calibEE = (TTree*) MapFile->Get("calibEE");
  if(!MapFile) { cout << "Invalid file: MapFile .. try again" << endl; return; }
  system( (string("mkdir -p ") + OutDir.Data() ).c_str() );
  TString nameFile = OutDir + "/h_DeadXtal.root";
  TFile* output = new TFile( nameFile.Data(), "RECREATE" );
  output->cd();
  TH2F* rms_EB  = new TH2F("rms_EB","DeadXtal #phi on x #eta on y",MAX_IPHI, MIN_IPHI, MAX_IPHI, 2*MAX_IETA+1, -MAX_IETA-0.5, MAX_IETA+0.5 );
  TH2F* rms_EB_r= new TH2F("rms_EB_r","DeadXtal #eta on x #phi on y", 2*MAX_IETA+1, -MAX_IETA-0.5, MAX_IETA+0.5, MAX_IPHI, MIN_IPHI, MAX_IPHI);
  TH2F* rms_EEp = new TH2F("rms_EEp","DeadXtal iX on x iY on y (EEp)",100,0.5,100.5,100,0.5,100.5);
  TH2F* rms_EEm = new TH2F("rms_EEm","DeadXtal iY on x iY on y (EEm)",100,0.5,100.5,100,0.5,100.5);

  //EB
  cout<<"Now EB..."<<endl;
  std::vector<int> Xtal_EB; Xtal_EB.clear();
  for(int nEB=0; nEB<61200; nEB++){
    //ostringstream convert; convert << nEB;
    string convert = std::to_string(nEB);
    TString name = "Barrel/epsilon_EB_iR_" + convert;
    TH1F * h1 = (TH1F*) fin->Get( name.Data() );
    if(h1->GetEntries()==0){
	Xtal_EB.push_back( nEB );
    }
  }
  //EE
  cout<<"Now EE..."<<endl;
  std::vector<int> Xtal_EE; Xtal_EE.clear();
  for(int nEE=0; nEE<14648; nEE++){
    //ostringstream convert; convert << nEE;
    string convert = std::to_string(nEE);
    TString name = "Endcap/epsilon_EE_iR_" + convert;
    TH1F * h1 = (TH1F*) fin->Get( name.Data() );
    if(h1->GetEntries()==0){
	Xtal_EE.push_back( nEE );
    }
  }
  //DrawPlots
  cout<<"Drawing EB."<<endl;
  Int_t ieta, iphi, hashedIndex; 
  calibEB->SetBranchAddress( "ieta_", &ieta);
  calibEB->SetBranchAddress( "iphi_", &iphi);
  calibEB->SetBranchAddress( "hashedIndex_", &hashedIndex);
  std::map<int,std::vector<int>> EBMap; EBMap.clear();
  Long64_t nentriesEB = calibEB->GetEntriesFast();
  for(Long64_t iEntry=0; iEntry<nentriesEB; iEntry++){
    calibEB->GetEntry(iEntry);
    std::vector<int> EtaPhi; EtaPhi.clear();
    EtaPhi.push_back( ieta ); EtaPhi.push_back( iphi );
    EBMap[hashedIndex] = EtaPhi;
  }
  for(int i=0; i<int(Xtal_EB.size()); i++){
    int BinEta = EBMap[Xtal_EB[i]][0] + 86;
    int BinPhi = EBMap[Xtal_EB[i]][1] + 1;
    rms_EB->SetBinContent( BinPhi, BinEta, 1 );
    rms_EB_r->SetBinContent( BinEta, BinPhi, 1 );
  }
  cout<<"And finally drawing EE."<<endl;
  Int_t ix, iy, zside, hashedIndex2; 
  calibEE->SetBranchAddress( "ix_", &ix);
  calibEE->SetBranchAddress( "iy_", &iy);
  calibEE->SetBranchAddress( "zside_", &zside);
  calibEE->SetBranchAddress( "hashedIndex_", &hashedIndex2);
  std::map<int,std::vector<int>> EEMap; EEMap.clear();
  Long64_t nentriesEE = calibEE->GetEntriesFast();
  for(Long64_t iEntry=0; iEntry<nentriesEE; iEntry++){
    calibEE->GetEntry(iEntry);
    std::vector<int> XYeZ; XYeZ.clear();
    XYeZ.push_back( ix ); XYeZ.push_back( iy ); XYeZ.push_back( zside );
    EEMap[hashedIndex2] = XYeZ;
  }
  for(int i=0; i<int(Xtal_EE.size()); i++){
    int BinX = EEMap[Xtal_EE[i]][0]+1;
    int BinY = EEMap[Xtal_EE[i]][1]+1;
    if( EEMap[Xtal_EE[i]][2]<0 ) rms_EEm->SetBinContent( BinX, BinY, 1 );
    if( EEMap[Xtal_EE[i]][2]>0 ) rms_EEp->SetBinContent( BinX, BinY, 1 );
  }
  //END
  cout<<"Done!!"<<endl;
  output->cd();
  rms_EB->Write();
  rms_EB_r->Write();
  rms_EEp->Write();
  rms_EEm->Write();
  output->Close();
}
