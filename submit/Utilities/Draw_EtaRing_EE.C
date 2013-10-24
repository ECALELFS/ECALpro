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
#include "TProfile.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::map;
using std::pair;
using std::stringstream;

struct iXiYtoRing {
  int iX;
  int iY;
  int sign;
  int Ring;
  iXiYtoRing() : iX(0), iY(0), sign(0), Ring(-1) { }
} GiveRing;

int GetRing(int x, int y, vector<iXiYtoRing> VectRing);

//Usage:
//.x Draw_EtaRing_EE.C+
void Draw_EtaRing_EE( bool EB=false, bool EE=true, bool noplot = true ){

  cout<<"Let's start with Parsing!"<<endl;
  //PARSING
  ifstream file;
  file.open("../common/Endc_x_y_ring.txt", ifstream::in);
  vector<iXiYtoRing> VectRing;

  while ( !file.eof() ) {
    string Line;
    getline( file, Line);
    string value;
    stringstream MyLine(Line);

    char * cstr, *p;
    cstr = new char [Line.size()+1];
    strcpy (cstr, Line.c_str());
    p=strtok (cstr," ");
    int i(0);
    while (p!=NULL){
	if(i==0)  GiveRing.iX = atoi(p);
	if(i==1)  GiveRing.iY = atoi(p);
	if(i==2)  GiveRing.sign = atoi(p);
	if(i==3){ 
	  GiveRing.Ring = atoi(p);
	  VectRing.push_back(GiveRing);
	}
	p=strtok(NULL," ");
	i++;
    }
    delete[] cstr;  
  }
  if(!noplot){
    //EE
    cout<<"Loop on the EE histo."<<endl;
    TH2F* h_EE = new TH2F("h_EE","EE etaring", 100, 0., 100., 100, 0., 100. );

    for(int x=0; x<100;x++){
	for(int y=0; y<100;y++){
	  int ring = GetRing( x, y, VectRing);
	  if( ring==-1 )     h_EE->SetBinContent(x,y,0);
	  else if(ring%2==0) h_EE->SetBinContent(x,y,100);
	  else               h_EE->SetBinContent(x,y,50);
	}
    }
    TCanvas* myc1 = new TCanvas("myc1", "CMS", 600, 600);
    h_EE->Draw("colz");
    myc1->SaveAs("EE_EtaRing.png");
    cout<<"Drawing..."<<endl;
    cout<<"Finally: THE END!"<<endl;
  }
  //Check on the Histo
  TFile* f_Ec  = TFile::Open( "../AfterCalibTools/WorkOnIC/prova_7Iter/ECAL/IC_SetAverageTo1_EtaRing.root" );
  if(!f_Ec) cout<<"Problem with File: EcalPro"<<endl;
  TFile* f_Ca  = TFile::Open( "../AfterCalibTools/WorkOnIC/prova_7Iter/CAL/IC_SetAverageTo1_EtaRing.root" );
  if(!f_Ca) cout<<"Problem with File: EtaRing"<<endl;
  TFile* f_Div = TFile::Open( "../AfterCalibTools/WorkOnIC/prova_7Iter/2012C_ECALCC_Vs_CalCC_7Iter_EtaRing1.root" );
  if(!f_Div) cout<<"Problem with File: division"<<endl;

  if(EB){
  cout<<"Now EB: "<<endl;
    TCanvas* myc1 = new TCanvas("myc1", "CMS", 600, 600);
    gStyle->SetOptStat(0);
    for(int nh=0; nh<=3; nh++){
	TProfile *pr_h = new TProfile("pr_h","profile", 171, 0., 171., 0.9, 1.1); pr_h->GetXaxis()->SetTitle("iEtaBin");
	TH2F* h;
	if(nh==0) h = (TH2F*) f_Ec->Get( "calibMap_EB" );
	if(nh==1) h = (TH2F*) f_Ca->Get( "calibMap_EB" );
	if(nh==2) h = (TH2F*) f_Div->Get( "Ratio_EBMap" );
	if(nh==3){h = (TH2F*) f_Ec->Get( "calibMap_EB" ); TH2F* h_2 = (TH2F*) f_Ca->Get( "calibMap_EB" ); h->Divide(h_2); }
	float Corr[171]={0.};
	float Norm[360]={0.};
	for(int iEta=0+1; iEta<171+1; iEta++){
	  for(int iPhi=0+1; iPhi<360+1; iPhi++){
	    if( iEta!=86 && h->GetBinContent(nh!=2?iEta:iPhi,nh!=2?iPhi:iEta)!=0 && h->GetBinContent(nh!=2?iEta:iPhi,nh!=2?iPhi:iEta)!=1 ){
		Corr[iEta-1] += h->GetBinContent(nh!=2?iEta:iPhi,nh!=2?iPhi:iEta);
		Norm[iEta-1] += 1.;
		if(nh==2) pr_h->Fill( float(iEta) , h->GetBinContent(nh!=2?iEta:iPhi,nh!=2?iPhi:iEta) );
		//cout<<iPhi<<" "<<iEta<<" : "<<h->GetBinContent(nh!=2?iEta:iPhi,nh!=2?iPhi:iEta)<<"  "<<Corr[iEta-1]<<"  "<<Norm[iEta-1]<<endl;
	    }
	  }
	}
	if(nh==2){ pr_h->Draw("colz"); myc1->SaveAs("Profile.png");}
	if(nh==0) cout<<"Now ECALpro"<<endl;
	if(nh==1) cout<<"Now Caltech"<<endl;
	if(nh==2) cout<<"Now Comparison"<<endl;
	if(nh==3) cout<<"Now the division"<<endl;
	for(int i=0; i<171;i++){
	  if(Norm[i]>0) cout<<Corr[i]/Norm[i]<<endl;
	}

    }
  }
  if(EE){
    cout<<"Now EE: "<<endl;
    for(int nh=0; nh<=3; nh++){
	TH2F* h;
	if(nh==0) h = (TH2F*) f_Ec->Get( "calibMap_EEp" ); 
	if(nh==1) h = (TH2F*) f_Ca->Get( "calibMap_EEp" ); 
	if(nh==2) h = (TH2F*) f_Div->Get( "Ratio_EEpMap" ); 
	if(nh==3){h = (TH2F*) f_Ec->Get( "calibMap_EEp" ); TH2F* h_2 = (TH2F*) f_Ca->Get( "calibMap_EEp" ); h->Divide(h_2); }
	float Corr[40]={0.};
	float Norm[40]={0.};
	for(int x=0; x<100;x++){
	  for(int y=0; y<100;y++){
	    int ring = GetRing( x, y, VectRing);
	    if( h->GetBinContent(x+1,y+1)!=0 && h->GetBinContent(x+1,y+1)!=1 ){
		//cout<<x<<" "<<y<<" "<<ring<<" "<<h->GetBinContent(x+1,y+1)<<endl;
		Corr[ring] += h->GetBinContent(x+1,y+1);
		Norm[ring] += 1.;
	    }
	  }
	}
	if(nh==0) cout<<"Now ECALpro"<<endl;
	if(nh==1) cout<<"Now Caltech"<<endl;
	if(nh==2) cout<<"Now Comparison"<<endl;
	if(nh==3) cout<<"Now the division"<<endl;
	for(int i=0; i<40;i++){
	  if(Norm[i]>0) cout<<Corr[i]/Norm[i]<<endl;
	}
    }
  }
}

//Find the ring
int GetRing(int x, int y, vector<iXiYtoRing> VectRing ){

  int index(0);
  bool NotFounf = false;
  for( size_t i=0; i<VectRing.size(); i++){
    if( VectRing[i].iX != x || VectRing[i].iY != y ) index++;
    if( VectRing[i].iX == x && VectRing[i].iY == y ) break;
    if( i==VectRing.size()-1 )                       NotFounf = true;
  }
  if(!NotFounf) return VectRing[index].Ring;
  else return -1;

}
