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

//.x MakeRatioIC.C+( "root://eoscms//eos/cms//store/caf/user/lpernie/ALL_2010ArelevantFiles_NOXTALWRONG_01/iter_4/calibMap.root",
//"root://eoscms//eos/cms/store/group/alca_ecalcalib/lpernie/ALL_2010_WithNEWSelection_02/iter_4/calibMap.root", "./Ratio_IC_2010/","Old_Selec_2012_selct.png" )
void MakeRatioIC( string inputFile0, string inputFile1, string OutPath, string OutPath1 )
{

    bool debug = false;
#ifdef DEBUG
    debug = true;
#endif

    //Out File
    system( (string("mkdir -p ") + OutPath ).c_str());
    string fileName = "";

    //In File
    TFile* f0 = TFile::Open( inputFile0.c_str() );
    if(!f0) {
	  cout << "Invalid file: " << inputFile0 << " .. try again" << endl;
	  return;
    }
    TFile* f1 = TFile::Open( inputFile1.c_str() );
    if(!f1) {
	  cout << "Invalid file: " << inputFile1 << " .. try again" << endl;
	  return;
    }

    TCanvas* myc1 = new TCanvas("myc1", "CMS", 600, 600);
    gStyle->SetOptStat(1111);
    TH2F* h0 = (TH2F*) f0->Get( "calibMap_EB" );
    TH2F* h1 = (TH2F*) f1->Get( "calibMap_EB" );

    TH1F *h_Div = new TH1F("h_Div","Ratio Between IC",100,0.99,1.02);

    cout<<"Computing Ratio...  "<<endl;

    for(int i=MIN_IETA; i< 2*MAX_IETA+2; i++){
	  for(int j=MIN_IPHI; j<MAX_IPHI+1; j++){
		if( h0->GetBinContent(i,j)!=1. && h1->GetBinContent(i,j)!=1. && i!=MAX_IETA+1 ){
		    if(debug) cout<<h0->GetBinContent(i,j)<<" / "<<h1->GetBinContent(i,j)<<endl;
		    h_Div->Fill(h0->GetBinContent(i,j)/h1->GetBinContent(i,j) );
		}
		else if(debug) cout<<"Bin is 1"<<endl;
	  }
    }
    myc1->cd();
    h_Div->Draw();
    string out = OutPath + OutPath1;
    myc1->SaveAs( out.c_str() );

    delete h_Div;
    delete myc1;
    cout<<"Done... You finished the Ratio of IC!!"<<endl;
}
