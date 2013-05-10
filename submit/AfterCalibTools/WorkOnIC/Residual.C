#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <set>
#include <cmath>
#include <cstdlib>
//Root Stuff
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TROOT.h"

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

//#define DEBUG
//#define DEBUG2
//#define DEBUG3

//Usage:
//.x Residual.C+("/store/group/alca_ecalcalib/lpernie/","ALL_2012B_NormSelect_NormFit_Merged_02", "11","11", "./2012B/", "2012merg_",true)
void Residual( string inputFile0, string inputFile , string NiterEB,string NiterEE,  string OutPath, string Tag, bool isPi0 )
{

    bool debug = false, debug2 = false, debug3 = false;
#ifdef DEBUG
    debug = true;
#endif
#ifdef DEBUG2
    debug2 = true;
#endif
#ifdef DEBUG3
    debug3 = true;
#endif

    //Out File
    system( (string("mkdir -p ") + OutPath ).c_str());
    string fileName = "";
    fileName = OutPath + "/IC_ResidsualFree.root";
    TFile* fout = new TFile(fileName.c_str(),"RECREATE");
    //for cross check
    TTree* Final_tree = new TTree("calibEB","Tree of EB Inter-calibration constants");
    Int_t   ieta1_, iphi1_;
    Float_t coeff1_;
    Final_tree->Branch( "coeff1_", &coeff1_, "coeff1_/F");
    Final_tree->Branch( "ieta1_", &ieta1_, "ieta1_/I");
    Final_tree->Branch( "iphi1_", &iphi1_, "iphi1_/I");  

    //In File
    string inputFile_mineEB = "root://eoscms//eos/cms" + inputFile0 + inputFile + "/iter_"+ NiterEB + "/" + Tag + "calibMap.root";
    string inputFile_mineEE = "root://eoscms//eos/cms" + inputFile0 + inputFile + "/iter_"+ NiterEE + "/" + Tag + "calibMap.root";
    if(isPi0) inputFile             = "root://eoscms//eos/cms/store/group/alca_ecalcalib/lpernie/ALL_2010_WithNEWSelection_02/iter_13/calibMap.root";
    else      inputFile             = "root://eoscms//eos/cms/store/group/alca_ecalcalib/lpernie/ALL_2010_forResid_01/iter_15/calibMap.root";
    TFile* f_base = TFile::Open( inputFile.c_str() );
    if(!f_base) {
	  cout << "Invalid file: " << inputFile << " .. try again" << endl;
	  return;
    }
    TFile* f_mineEB = TFile::Open( inputFile_mineEB.c_str() );
    if(!f_mineEB) {
	  cout << "Invalid file: " << inputFile_mineEB << " .. try again" << endl;
	  return;
    }
    TFile* f_mineEE = TFile::Open( inputFile_mineEE.c_str() );
    if(!f_mineEE) {
	  cout << "Invalid file: " << inputFile_mineEE << " .. try again" << endl;
	  return;
    }
    //Histos
    TCanvas* myc1 = new TCanvas("myc1", "CMS", 600, 600);
    gStyle->SetOptStat(0);
    TH2F* h = (TH2F*) f_mineEB->Get( "calibMap_EB" );
    TH2F* h_base = (TH2F*) f_base->Get( "calibMap_EB" );
    TH2F* h_out = new TH2F("calibMap_EB","EB calib coefficients: #eta on x, #phi on y", 2*MAX_IETA+1, -MAX_IETA-0.5, MAX_IETA+0.5, MAX_IPHI, MIN_IPHI-0.5, MAX_IPHI+0.5 );

    cout<<"Computing Residuals...  "<<endl;
    //F2_a: correction for iPhi; iEta<55
    float f1_iEta[85]={0}; 
    float mediaA[85]={0};
    float f2a_iPhi[20]={0}, f2b_iPhi[20]={0};
    float mediaB[20]={0}, mediaC[20]={0};
    //float f2E_iPhi[20]={0}, mediaE[20]={0};
    float meanEta(0),meanEtaT(0);
    float meanPhi1(0),meanPhi1T(0);
    float meanPhi2(0),meanPhi2T(0);
    for(int i=MIN_IETA; i< 2*MAX_IETA+2; i++){
	  for(int j=MIN_IPHI; j<MAX_IPHI+1; j++){
		if( h_base->GetBinContent(i,j)!=1. && i!=MAX_IETA+1 ){
		    f1_iEta[ (TMath::Abs( i-(MAX_IETA+1) ))-1 ] += h_base->GetBinContent( i, j );
		    mediaA[ (TMath::Abs( i-(MAX_IETA+1) ))-1 ]++;
		    meanEta+= h_base->GetBinContent( i, j ); meanEtaT++;
		    if((TMath::Abs( i-(MAX_IETA+1) )) < 55 ){
			  f2a_iPhi[ j%20 ] += h_base->GetBinContent( i, j ); 
			  mediaB[j%20]++; 
			  meanPhi1+= h_base->GetBinContent( i, j );
			  meanPhi1T++;
		    }
		    else                                    { 
			  f2b_iPhi[ j%20 ] += h_base->GetBinContent( i, j ); 
			  mediaC[j%20]++;
			  meanPhi2+= h_base->GetBinContent( i, j );
			  meanPhi2T++;}
		}
	  }
    }
    for(int i=0; i<85; i++){
	  f1_iEta[ i ] = (meanEta/meanEtaT)/(f1_iEta[ i ]/mediaA[ i ]);
	  if( debug ) cout<<"My weightA: "<<mediaA[ i ]<<" "<<f1_iEta[ i ]<<endl;
    }

    for(int i=0; i<20; i++){ 
	  f2a_iPhi[ i ] = (meanPhi1/meanPhi1T)/(f2a_iPhi[ i ]/mediaB[ i ]);
	  f2b_iPhi[ i ] = (meanPhi2/meanPhi2T)/(f2b_iPhi[ i ]/mediaC[ i ]);
	  if( debug ) cout<<"My weightB: "<<mediaB[ i ]<<" "<<f2a_iPhi[ i ]<<endl;
	  if( debug ) cout<<"My weightC: "<<mediaC[ i ]<<" "<<f2b_iPhi[ i ]<<endl;
    }

    // EB Residual correction
    cout<<"Let's start the corrections...  "<<endl;

    float med=0., med2=0., medT=0.;
    float media2[85]={0}, media2X[85]={0};
    for(int i=MIN_IETA; i< 2*MAX_IETA+2; i++){
	  for(int j=MIN_IPHI; j<MAX_IPHI+1; j++){
		if( h->GetBinContent(i,j)!=1. && i!=MAX_IETA+1 ){
		    float    corr  = f1_iEta[ (TMath::Abs( i-(MAX_IETA+1) ))-1 ];
		    if( (TMath::Abs(i-(MAX_IETA+1))) < 55 ) corr *= f2a_iPhi[(j%20)];
		    else                                    corr *= f2b_iPhi[(j%20)];
		    if( debug2 && i==MIN_IETA ) cout<<i<<" "<<j<<" ->Correction: "<<f1_iEta[ (TMath::Abs( i-(MAX_IETA+1) ))-1 ]<<" * "<<f2a_iPhi[j%20]<<" * "<<f2b_iPhi[j%20]<<" so...: "<<corr<<" (to be applied at: "<<h->GetBinContent(i,j)<<" )"<<endl;
		    h_out->SetBinContent( i, j, h->GetBinContent(i,j)*corr );
		    coeff1_ = h->GetBinContent(i,j)*corr;
		    medT++; med+=h->GetBinContent(i,j)*corr; med2+=h->GetBinContent(i,j);
		    ieta1_ = (TMath::Abs( i-(MAX_IETA+1) )) ; iphi1_ = j;
		    Final_tree->Fill();
		}
		else h_out->SetBinContent( i, j, h->GetBinContent(i,j));
	  }
	  if(i!=MAX_IETA+1 && debug2 ) cout<<"BIN: "<<med2/medT<<"  "<<med/medT<<endl;
	  media2[ (TMath::Abs( i-(MAX_IETA+1) ))-1 ]=med/medT; media2X[ (TMath::Abs( i-(MAX_IETA+1) ))-1 ]=(TMath::Abs( i-(MAX_IETA+1) ));
	  med=0.; medT=0.; med2=0;
    }
    TGraph *gr = new TGraph(85,media2X,media2);
    gr->Draw("lta");
    string saving = OutPath + "OldResidual_Eta.png";
    myc1->SaveAs(saving.c_str());

    if(debug3) cout<<""<<endl;
    //Just Rewrite EE now
    cout<<"Now rewrite EE...  "<<endl;
    TH2F* hep = (TH2F*) f_mineEE->Get( "calibMap_EEp" );
    TH2F* hem = (TH2F*) f_mineEE->Get( "calibMap_EEm" );

    fout->cd();
    h_out->Write();
    hep->Write();
    hem->Write();
    Final_tree->Write();
    fout->Close();

    delete fout; delete h_out;
    delete hep;  delete hem;
    cout<<"Done... You finished the Residual Correction!!"<<endl;
}
