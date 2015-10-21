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
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"

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
//#define DEBUG4

struct iXiYtoRing {
  int iX;
  int iY;
  int sign;
  int Ring;
  iXiYtoRing() : iX(0), iY(0), sign(0), Ring(-1) { }
} GiveRing;

int GetRing(int x, int y, vector<iXiYtoRing> VectRing, bool debug3);

//Usage:
//.x SetIC1toiEta.C+(false,"ALL_Neutrino_Pt2to20_AVE40BX25_FoldEtaRing_03/Absolute_IC_AlphaStudies_2012.root","ALL_Neutrino_Pt2to20_AVE40BX25_FoldEtaRing_03","",0,"")
//.x SetIC1toiEta.C+(true(global)?false(etaRIng),"2011/ABSIC_ResidsualFree.root", "./2011/" "2012D_ETA_ERR/Error_Stat_2012D.root")
//.x SetIC1toiEta.C+(false,"Compare2012B/2012B_GlobalMY_Vs_Caltech.root", "Compare2012B/", 2 ) [0=not for comparison mine-caltech, 1 for comp., is mine, 2 for comp. is clatech]
void SetIC1toiEta( bool AllEta, string inputFile,  TString OutPath, string PutStatError="", int isCal=0 , string IOV=""){
  //AllEta: true -> mediated to 1 globally; false -> mediated to 1 per EtaRing
  //inputFile -> root file with histos
  //OutPath -> output folder

  bool debug = false, debug2 = false, debug3 = false, debug4 = false;
#ifdef DEBUG
  debug = true;
#endif
#ifdef DEBUG2
  debug2 = true;
#endif
#ifdef DEBUG3
  debug3 = true;
#endif
#ifdef DEBUG4
  debug4 = true;
#endif

  //Output File
  system( (string("mkdir -p ") + OutPath ).Data());
  string fileName = "";    
  if(AllEta){
    if(isCal==0)  fileName = OutPath + "/IC_SetAverageTo1_global.root";
    if(isCal==1)  fileName = OutPath + "/ICForComp_SetAverageTo1_global_Mine.root";
    if(isCal==2)  fileName = OutPath + "/ICForComp_SetAverageTo1_global_Caltech.root";
  }
  else{  
    if(isCal==0)  fileName = OutPath + "/IC_SetAverageTo1_EtaRing.root";
    if(isCal==1)  fileName = OutPath + "/ICForComp_SetAverageTo1_EtaRing_Mine.root";
    if(isCal==2)  fileName = OutPath + "/ICForComp_SetAverageTo1_EtaRing_Caltech.root";
  }
  TFile* fout = new TFile(fileName.c_str(),"RECREATE");
  fstream  f_IC;
  TString NameTxt="";
  if(AllEta){
    if(isCal==0)  NameTxt = "/IC_Ecal_global";
    if(isCal==1)  NameTxt = "/ICForComp_Ecal_global_Mine";
    if(isCal==2)  NameTxt = "/ICForComp_Ecal_global_Caltech";
  }
  else  {
    if(isCal==0)  NameTxt = "/IC_Ecal_EtaRing";
    if(isCal==1)  NameTxt = "/ICForComp_Ecal_EtaRing_Mine";
    if(isCal==2)  NameTxt = "/ICForComp_Ecal_EtaRing_Caltech";
  }
  if(PutStatError == "") NameTxt += "_noErr";
  if(IOV != "") NameTxt += "_" + IOV;
  f_IC.open( (OutPath + NameTxt) + ".txt", ios::out);
  if( !f_IC ){ cout << "Impossible to open file.txt."; exit(1); }

  //Open File
  TFile* f = TFile::Open( inputFile.c_str() );
  if(!f) {
    cout << "Invalid file: " << inputFile << " .. try again" << endl;
    return;
  }
  //Opern Error File
  TFile* f_err;
  if(PutStatError != "" ) f_err = TFile::Open( PutStatError.c_str() );
  TH2F* h_errEB; TH2F* h_errEEm; TH2F* h_errEEp;
  if(PutStatError != "" ){
    h_errEB  = (TH2F*) f_err->Get( "Sig_EB" );
    h_errEEm = (TH2F*) f_err->Get( "Sig_EEm" );
    h_errEEp = (TH2F*) f_err->Get( "Sig_EEp" );
  }

  if(AllEta){
    TString      name="calibMap_EB";
    if(isCal==1) name="calibMap_EB_1";
    if(isCal==2) name="calibMap_EB_2";
    TH2F* h = (TH2F*) f->Get( name.Data() );
    TH2F* h_out = new TH2F("calibMap_EB","EB calib coefficients: #eta on x, #phi on y", 2*MAX_IETA+1, -MAX_IETA-0.5, MAX_IETA+0.5, MAX_IPHI, MIN_IPHI-0.5, MAX_IPHI+0.5 );
    cout<<"Let's start with EB!"<<endl;

    //Loop on iETA
    double Ic_tmp(0.);
    int    Ic_tot(0.);
    for(int i=MIN_IETA; i< 2*MAX_IETA+2; i++){ //from =1 to <172
	for(int j=MIN_IPHI; j<MAX_IPHI+1; j++){  //from =1 to <361
	  if( h->GetBinContent(i,j)==0 && i !=  MAX_IETA+1 )  cout<<"WARNING: IC in EB is ZERO!!"<<endl;
	  if( h->GetBinContent(i,j) != 1. && h->GetBinContent(i,j) != 0. ){
	    Ic_tmp += h->GetBinContent(i,j);
	    Ic_tot ++;
	  }
	  if(debug) cout<<"IC: iETA "<<i<<" iPHI "<<j<<" = "<<h->GetBinContent(i,j)<<endl;
	}
    }
    cout<<"EB: iC!=1 are: "<<Ic_tot<<endl;
    cout<<"Ic_tmp: "<< Ic_tmp <<" Ic_tot: "<<Ic_tot<<" ----> Weight: "<<Ic_tmp/Ic_tot<<endl;
    Ic_tmp/=Ic_tot;
    //Now Correct it:
    for(int i=MIN_IETA; i< 2*MAX_IETA+2; i++){
	for(int j=MIN_IPHI; j<MAX_IPHI+1; j++){
	  if( h->GetBinContent(i,j) != 1 && h->GetBinContent(i,j) != 0 ) h_out->SetBinContent( i, j, h->GetBinContent(i,j)*(1/Ic_tmp));
	  else                                                           h_out->SetBinContent( i, j, h->GetBinContent(i,j) );
	}
    }
    //Cross check:
    double Ic_tmp1(0.);
    int    Ic_tot1(0.);
    for(int i=MIN_IETA; i< 2*MAX_IETA+2; i++){
	for(int j=MIN_IPHI; j<MAX_IPHI+1; j++){
//	  if( i != MAX_IETA+1 ){
//	    Ic_tmp1 += h_out->GetBinContent(i,j);
//	    Ic_tot1++;
//	  }
	  if( h_out->GetBinContent(i,j) != 1 &&  i != MAX_IETA+1){ Ic_tmp1 += h_out->GetBinContent(i,j); Ic_tot1++;  }
	  if( i != MAX_IETA+1 ){
	    float Error = 0;
	    if(PutStatError != "") Error = h_errEB->GetBinContent(i,j);
	    if(h_out->GetBinContent(i,j)==1. ){  f_IC << i-(MAX_IETA+1) << " "<< j << " 0 " << "-1." << " 999."<<endl;}
	    else                              {  f_IC << i-(MAX_IETA+1) << " "<< j << " 0 " << h_out->GetBinContent(i,j) << " " << Error <<endl;}
	  }
	}
    }
    cout<<"FIRST Cross check: Mean is "<<float(Ic_tmp1/float(Ic_tot1))<<endl;

    //Save Final File
    fout->cd();
    h_out->Write();

    //EEp
    name="calibMap_EEp";
    if(isCal==1) name="calibMap_EEp_1";
    if(isCal==2) name="calibMap_EEp_2";
    TH2F* hep = (TH2F*) f->Get( name.Data() );
    TH2F* hmap_EEp = new TH2F("calibMap_EEp","EE+ calib coefficients",100,0.5,100.5,100,0.5,100.5);
    hmap_EEp->GetXaxis()->SetTitle("ix");  hmap_EEp->GetYaxis()->SetTitle("iy");
    double Ic_tmpE(0.);
    int    Ic_totE(0.);
    for(int x=0; x<100;x++){
	for(int y=0; y<100;y++){
	  if( hep->GetBinContent(x+1,y+1) != 0 && hep->GetBinContent(x+1,y+1) != 1){
	    Ic_tmpE+=hep->GetBinContent(x+1,y+1);
	    Ic_totE++;
	  }
	}
    }
    cout<<"Ic_tmpE: "<< Ic_tmpE <<" Ic_totE: "<<Ic_totE<<" ----> Weight: "<<Ic_tmpE/Ic_totE<<endl;
    Ic_tmpE/=Ic_totE;
    //Now Correct it:
    for(int x=0; x<100;x++){
	for(int y=0; y<100;y++){
	  if( hep->GetBinContent(x+1,y+1) != 0 && hep->GetBinContent(x+1,y+1) != 1){
	    hmap_EEp->SetBinContent( x+1, y+1, hep->GetBinContent(x+1,y+1)/Ic_tmpE );
	  }
	  else hmap_EEp->SetBinContent( x+1, y+1, hep->GetBinContent(x+1,y+1) );
	}
    }
    //Cross check:
    double Ic_tmpE1(0.);
    int    Ic_totE1(0.);
    for(int x=0; x<100;x++){
	for(int y=0; y<100;y++){
	  if( hmap_EEp->GetBinContent(x+1,y+1) != 0 && hmap_EEp->GetBinContent(x+1,y+1) != 1){
	    Ic_tmpE1+=hmap_EEp->GetBinContent(x+1,y+1); Ic_totE1++;
	  }
	  if( hmap_EEp->GetBinContent(x+1,y+1) != 0){
	    float Error = 0;
	    if(PutStatError != "") Error = h_errEEp->GetBinContent(x+1,y+1);
	    if(hmap_EEp->GetBinContent(x+1,y+1)!=1 ){ f_IC << x+1 << " "<< y+1 << " 1 " << hmap_EEp->GetBinContent(x+1,y+1) << " " << Error <<endl;}
	    else                                    { f_IC << x+1 << " "<< y+1 << " 1 " << "-1." << " 999."<<endl;}
	  }
	}
    }
    cout<<"SECOND Cross check: Mean is "<<float(Ic_tmpE1/float(Ic_totE1))<<endl;

    //EEm
    name="calibMap_EEm";
    if(isCal==1) name="calibMap_EEm_1";
    if(isCal==2) name="calibMap_EEm_2";
    TH2F* hem = (TH2F*) f->Get( name.Data() );
    TH2F* hmap_EEm = new TH2F("calibMap_EEm","EE- calib coefficients",100,0.5,100.5,100,0.5,100.5);
    hmap_EEm->GetXaxis()->SetTitle("ix");  hmap_EEm->GetYaxis()->SetTitle("iy");
    double Ic_tmpEm(0.);
    int    Ic_totEm(0.);
    for(int x=0; x<100;x++){
	for(int y=0; y<100;y++){
	  if( hem->GetBinContent(x+1,y+1) != 0 && hem->GetBinContent(x+1,y+1) != 1){
	    Ic_tmpEm+=hem->GetBinContent(x+1,y+1);
	    Ic_totEm++;
	  }
	}
    }
    cout<<"Ic_tmpEm: "<< Ic_tmpEm <<" Ic_totEm: "<<Ic_totEm<<" ----> Weight: "<<Ic_tmpEm/Ic_totEm<<endl;
    Ic_tmpEm/=Ic_totEm;
    //Now Correct it:
    for(int x=0; x<100;x++){
	for(int y=0; y<100;y++){
	  if( hem->GetBinContent(x+1,y+1) != 0 && hem->GetBinContent(x+1,y+1) != 1){
	    hmap_EEm->SetBinContent( x+1, y+1, hem->GetBinContent(x+1,y+1)/Ic_tmpEm );
	  }
	  else hmap_EEm->SetBinContent( x+1, y+1, hem->GetBinContent(x+1,y+1) );
	}
    }
    //Cross check:
    double Ic_tmpE1m(0.);
    int    Ic_totE1m(0.);
    for(int x=0; x<100;x++){
	for(int y=0; y<100;y++){
	  if( hmap_EEm->GetBinContent(x+1,y+1) != 0 && hmap_EEm->GetBinContent(x+1,y+1) != 1){
	    Ic_tmpE1m+=hmap_EEm->GetBinContent(x+1,y+1); Ic_totE1m++;
	  }
	  if( hmap_EEm->GetBinContent(x+1,y+1) != 0){
	    float Error = 0;
	    if(PutStatError != "") Error = h_errEEm->GetBinContent(x+1,y+1);
	    if(hmap_EEm->GetBinContent(x+1,y+1)!=1){ f_IC << x+1 << " "<< y+1 << " -1 " << hmap_EEm->GetBinContent(x+1,y+1) << " " << Error <<endl;}
	    else                                   { f_IC << x+1 << " "<< y+1 << " -1 " << "-1." << " 999."<<endl;}
	  }
	}
    }
    cout<<"THIRD Cross check: Mean is "<<float(Ic_tmpE1m/float(Ic_totE1m))<<endl;
    //Write Output
    fout->cd();
    hmap_EEm->Write();
    hmap_EEp->Write();
  }
  else{

    TString      name="calibMap_EB";
    if(isCal==1) name="calibMap_EB_1";
    if(isCal==2) name="calibMap_EB_2";
    TH2F* h = (TH2F*) f->Get( name.Data() );
    TH2F* h_out = new TH2F("calibMap_EB","EB calib coefficients: #eta on x, #phi on y", 2*MAX_IETA+1, -MAX_IETA-0.5, MAX_IETA+0.5, MAX_IPHI, MIN_IPHI-0.5, MAX_IPHI+0.5 );

    cout<<"Let's start with EB!"<<endl;

    //Loop on iETA
    int nEB_1(0);
    for(int i=MIN_IETA; i< 2*MAX_IETA+2; i++){ //+1 prima 
	double Ic_tmp(0.);
	int    Ic_tot(0.);

	//Average on iPhi
	for(int j=MIN_IPHI; j<MAX_IPHI+1; j++){
	  if( h->GetBinContent(i,j) != 1 ){
	    if( h->GetBinContent(i,j)==0 && i !=  MAX_IETA+1 ) cout<<"WARNING: IC in EB is ZERO!!"<<endl;
	    Ic_tmp += h->GetBinContent(i,j);
	    Ic_tot ++;
	  }
	  if(debug) cout<<"IC: iETA "<<i<<" iPHI "<<j<<" = "<<h->GetBinContent(i,j)<<endl;
	}

	if(debug) cout<<"--> Sum of IC: "<<Ic_tmp<<" The average is: "<<Ic_tmp/Ic_tot<<endl;
	Ic_tmp/=Ic_tot;

	for(int j=MIN_IPHI; j<MAX_IPHI+1; j++){
	  if( h->GetBinContent(i,j) != 1 ) h_out->SetBinContent( i, j, h->GetBinContent(i,j)*(1/Ic_tmp));
	  else                             h_out->SetBinContent( i, j, h->GetBinContent(i,j) );
	  if(debug) cout<<"NEW IC: iETA "<<i<<" iPHI "<<j<<" =>  "<<h->GetBinContent(i,j)<<" * "<<Ic_tmp<<" = "<<h->GetBinContent(i,j)*Ic_tmp<<endl;
	  if( i != MAX_IETA+1 ){
	    float Error = 0;
	    if(PutStatError != "") Error = h_errEB->GetBinContent(i,j);
	    if(h->GetBinContent(i,j)==1 ){ nEB_1++;  f_IC << i-(MAX_IETA+1) << " "<< j << " 0 " << "-1." << " 999."<<endl;}
	    else                                  {  f_IC << i-(MAX_IETA+1) << " "<< j << " 0 " << h->GetBinContent(i,j)*(1/Ic_tmp) << " " << Error <<endl;}
	  }
	}
    }
    cout<<"EB: iC==1 are: "<<nEB_1<<endl;
    //Save Final File
    fout->cd();
    h_out->Write();

    cout<<"Now EE! Befor Parsing..."<<endl;
    //Now the EE-----------------------------------------------------

    //PARSING
    ifstream file;
    file.open("/afs/cern.ch/work/l/lpernie/ECALpro/gitHubCalib/CMSSW_5_3_6/src/CalibCode/submit/AfterCalibTools/WorkOnIC/InputFile/Endc_x_y_ring.txt", ifstream::in);
    vector<iXiYtoRing> VectRing;

    while ( !file.eof() ) {
	string Line;
	getline( file, Line);
	if(debug2) cout<<Line<<endl;
	string value;
	stringstream MyLine(Line);

	char * cstr, *p;
	cstr = new char [Line.size()+1];
	strcpy (cstr, Line.c_str());
	p=strtok (cstr," ");
	int i(0);
	while (p!=NULL){
	  if(debug2) cout<<"--->"<< i<<" " << p << endl;
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

    if(debug2){
	cout<<"vec size "<<VectRing.size()<<endl;
	for(size_t j=0; j<VectRing.size(); j++ ){
	  cout<< VectRing[j].iX<<" "<<VectRing[j].iY<<" "<<VectRing[j].sign<<" "<<VectRing[j].Ring<<endl;
	}
    }

    cout<<"and now EE plus"<<endl;

    //EE PLUS
    name="calibMap_EEp";
    if(isCal==1) name="calibMap_EEp_1";
    if(isCal==2) name="calibMap_EEp_2";
    TH2F* hep = (TH2F*) f->Get( name.Data() );
    TH2F* hmap_EEp = new TH2F("calibMap_EEp","EE+ calib coefficients",100,0.5,100.5,100,0.5,100.5);
    hmap_EEp->GetXaxis()->SetTitle("ix");
    hmap_EEp->GetYaxis()->SetTitle("iy");
    float Ringweight_p[40]={0.};
    float RingNum_p[40]={0.};

    for(int x=0; x<100;x++){
	for(int y=0; y<100;y++){
	  if( hep->GetBinContent(x+1,y+1) != 0 && hep->GetBinContent(x+1,y+1) != 1){
	    if(debug3) cout<<"Cerco: x "<<x<<" y "<<y<<" "<<hep->GetBinContent(x+1,y+1)<<endl;
	    int ring = GetRing( x, y, VectRing, debug3);
	    Ringweight_p[ring]+=hep->GetBinContent(x+1,y+1);
	    RingNum_p[ring]++;
	  }
	}
    }
    //Make the average
    for(int j=0; j<40;j++){ if(debug3) cout<<"Index: "<<j<<" Ringweight" <<Ringweight_p[j]<<" RingNum "<<RingNum_p[j]<<endl; if(RingNum_p[j]!=0) Ringweight_p[j]/=RingNum_p[j];  }

    float VAR[50]={0.}, TOT[50]={0.};
    //Fill the new Histo
    int nEEp_1(0);
    for(int x=0; x<100;x++){
	for(int y=0; y<100;y++){
	  if( hep->GetBinContent(x+1,y+1) != 0 && hep->GetBinContent(x+1,y+1) != 1){
	    int ring = GetRing( x, y, VectRing, debug3);

	    VAR[ring] += hep->GetBinContent(x+1,y+1)*(1/Ringweight_p[ring]);
	    TOT[ring]++;

	    hmap_EEp->SetBinContent( x+1, y+1, hep->GetBinContent(x+1,y+1)*(1/Ringweight_p[ring]) );
	    if(debug3) cout<<"EEp riempi con: "<<hep->GetBinContent(x+1,y+1)<<" * "<<1./Ringweight_p[ring]<<endl;
	  }
	  else hmap_EEp->SetBinContent( x+1, y+1, hep->GetBinContent(x+1,y+1) );
	  if( hep->GetBinContent(x+1,y+1) != 0){
	    float Error = 0;
	    if(PutStatError != "") Error = h_errEEp->GetBinContent(x+1,y+1);
	    if( hep->GetBinContent(x+1,y+1)!=1 ){ f_IC << x+1 << " "<< y+1 << " 1 " << hmap_EEp->GetBinContent(x+1,y+1) << " " << Error <<endl;}
	    else                      { nEEp_1++; f_IC << x+1 << " "<< y+1 << " 1 " << "-1." << " 999."<<endl;}
	  }
	}
    }
    cout<<"EEp: iC==1 are: "<<nEEp_1<<endl;
    //VEDI SE SONO MEDIATE
    for(int i=0; i<50;i++){
	if( TOT[i]!=0 && fabs(VAR[i]/TOT[i]-1.)>0.001 ) cout<<"WARNING: MEDIA EE+ "<<i<<"): "<<VAR[i]/TOT[i]<<endl;
    }

    //Save Final File
    fout->cd();
    hmap_EEp->Write();

    cout<<"then finally EE minus"<<endl;
    //EE MINUS
    name="calibMap_EEm";
    if(isCal==1) name="calibMap_EEm_1";
    if(isCal==2) name="calibMap_EEm_2";
    TH2F* hem = (TH2F*) f->Get( name.Data() );
    TH2F* hmap_EEm = new TH2F("calibMap_EEm","EE- calib coefficients",100,0.5,100.5,100,0.5,100.5);
    hmap_EEm->GetXaxis()->SetTitle("ix");
    hmap_EEm->GetYaxis()->SetTitle("iy");
    float Ringweight_m[40]={0.};
    float RingNum_m[40]={0.};

    for(int x=0; x<100;x++){
	for(int y=0; y<100;y++){
	  if( hem->GetBinContent(x+1,y+1) != 0 && hem->GetBinContent(x+1,y+1) != 1){
	    if(debug3) cout<<"Cerco: x "<<x<<" y "<<y<<" "<<hem->GetBinContent(x+1,y+1)<<endl;
	    int ring = GetRing( x, y, VectRing, debug3);
	    Ringweight_m[ring]+=hem->GetBinContent(x+1,y+1);
	    RingNum_m[ring]++;
	  }
	}
    }
    //Make the average
    for(int j=0; j<40;j++){ if(debug3) cout<<"Index: "<<j<<" Ringweight" <<Ringweight_m[j]<<" RingNum "<<RingNum_m[j]<<endl; if(RingNum_m[j]!=0) Ringweight_m[j]/=RingNum_m[j];  }

    float VAR1[50]={0.}, TOT1[50]={0.};
    //Fill the new Histo
    int nEEm_1(0);
    for(int x=0; x<100;x++){
	for(int y=0; y<100;y++){
	  if( hem->GetBinContent(x+1,y+1) != 0 && hem->GetBinContent(x+1,y+1) != 1){
	    int ring = GetRing( x, y, VectRing, debug3);

	    VAR1[ring] += hem->GetBinContent(x+1,y+1)*(1/Ringweight_m[ring]);
	    TOT1[ring]++;

	    hmap_EEm->SetBinContent( x+1, y+1, hem->GetBinContent(x+1,y+1)*(1/Ringweight_m[ring]) );
	    if(debug3) cout<<"EEm riempi con: "<<hem->GetBinContent(x+1,y+1)<<" * "<<1./Ringweight_m[ring]<<endl;
	  }
	  else hmap_EEm->SetBinContent( x+1, y+1, hem->GetBinContent(x+1,y+1) );
	  if( hem->GetBinContent(x+1,y+1) != 0){
	    float Error = 0;
	    if(PutStatError != "") Error = h_errEEm->GetBinContent(x+1,y+1);
	    if(hem->GetBinContent(x+1,y+1)!=1 ){ f_IC << x+1 << " "<< y+1 << " -1 " << hmap_EEm->GetBinContent(x+1,y+1) << " " << Error <<endl;}
	    else                     { nEEm_1++; f_IC << x+1 << " "<< y+1 << " -1 " << "-1." << " 999."<<endl;}
	  }
	}
    }
    cout<<"EEm: iC==1 are: "<<nEEm_1<<endl;

    //VEDI SE SONO MEDIATE
    for(int i=0; i<50;i++){
	if( TOT1[i]!=0 && fabs(VAR1[i]/TOT1[i]-1.)>0.001 ) cout<<"WARNING: MEDIA EE- "<<i<<"): "<<VAR1[i]/TOT1[i]<<endl;
    }

    //Save Final File
    fout->cd();
    hmap_EEm->Write();

    cout<<"Finish with last cross check!"<<endl;
    //Last Control check
    if(debug4){
	//EB
	int EBN_one(0), EBN_null(0), EBN_tot(0);
	for(int i=MIN_IETA; i< 2*MAX_IETA+1; i++){
	  for(int j=MIN_IPHI; j<MAX_IPHI+1; j++){
	    if(h_out->GetBinContent(i,j)!=0) cout<<"EB Diff IC: "<<h->GetBinContent(i,j) - h_out->GetBinContent(i,j)<<" Ratio: "<<h->GetBinContent(i,j) / h_out->GetBinContent(i,j)<<endl;
	    if(h->GetBinContent(i,j)==1) EBN_one++;
	    if(h->GetBinContent(i,j)==0) EBN_null++;
	    EBN_tot++;
	  }
	}
	cout<<"EBN_tot: "<<EBN_tot<<" EBN_one: "<<EBN_one<<" EBN_null: "<<EBN_null<<endl;
	//EE
	int EEmN_one(0), EEpN_one(0), EEmN_null(0), EEpN_null(0), EE_tot(0);
	for(int x=0; x<100;x++){
	  for(int y=0; y<100;y++){
	    if(hmap_EEm->GetBinContent( x+1, y+1 )!=0)
		cout<<"EE m Diff IC: "<< hem->GetBinContent(x+1,y+1) - hmap_EEm->GetBinContent( x+1, y+1 )<<" Ratio: "<< hem->GetBinContent(x+1,y+1) / hmap_EEm->GetBinContent( x+1, y+1 ) <<endl;
	    if(hmap_EEp->GetBinContent( x+1, y+1 )!=0)
		cout<<"EE p Diff IC: "<< hep->GetBinContent(x+1,y+1) - hmap_EEp->GetBinContent( x+1, y+1 )<<" Ratio: "<< hep->GetBinContent(x+1,y+1) / hmap_EEp->GetBinContent( x+1, y+1 ) <<endl;
	    if(hem->GetBinContent(x+1,y+1)==1) EEmN_one++;
	    if(hep->GetBinContent(x+1,y+1)==1) EEpN_one++;
	    if(hem->GetBinContent(x+1,y+1)==0) EEmN_null++;
	    if(hep->GetBinContent(x+1,y+1)==0) EEpN_null++;
	    EE_tot++;
	  }
	}
	cout<<"EE_tot M: "<<EE_tot<<" EEmN_one: "<<EEmN_one<<" EEmN_null: "<<EEmN_null<<endl;
	cout<<"EE_tot P: "<<EE_tot<<" EEpN_one: "<<EEpN_one<<" EEpN_null: "<<EEpN_null<<endl;
    }
  }//if not global mean

  cout<<"Finally: THE END!"<<endl;
  fout->Close();
  f_IC.close();

}


//Find the ring
int GetRing(int x, int y, vector<iXiYtoRing> VectRing, bool debug3){

  int index(0);
  for( size_t i=0; i<VectRing.size(); i++){
    if(  VectRing[i].iX != x || VectRing[i].iY != y ) index++;
    if(  VectRing[i].iX == x && VectRing[i].iY == y ) break;
  }
  if(debug3) cout<<"Ho trovato: "<<VectRing[index].iX<<" e "<<VectRing[index].iY<<endl;
  return VectRing[index].Ring;

}
