#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "TTree.h"
#include "TLatex.h"
#include "TMath.h"
#include "TBranch.h"
#include "TFile.h"
#include "TStyle.h"
#include "TPad.h"
#include "TProfile.h"
#include "TPave.h"
#include "TPaveText.h"
#include "TString.h"
#include "TPaveStats.h"

using namespace std;

//.x Step3_BestBin.C+("ALL_MINBIAS_UNCAL_L1_NOL1FILTER/Fstep2_notSorted_EB_pi0.txt", 8382, "ALL_MINBIAS_UNCAL_L1_NOL1FILTER/Fstep1_BINcut_EB_pi0.txt",8400)
void Step3_BestBin(TString Input, int nLine, TString Input2, int nLine2){

  //Create an array with all the values
  vector< float > values_effL; values_effL.clear();
  vector< float > values_sbL; values_sbL.clear();
  vector< float > values_muL; values_muL.clear();
  vector< float > values_chiL; values_chiL.clear();
  vector< int >   values_binL; values_binL.clear();
  vector< float > values_effH; values_effH.clear();
  vector< float > values_sbH; values_sbH.clear();
  vector< float > values_muH; values_muH.clear();
  vector< float > values_chiH; values_chiH.clear();
  vector< int >   values_binH; values_binH.clear();
  string Line, check;
  ifstream infile;
  infile.open (Input.Data());
  int nL=1;
  while(!infile.eof()){
    if(nL==nLine) break;
    nL++;
    getline(infile,Line); 
    string buf; 
    stringstream ss(Line); 
    vector<string> tokens; tokens.clear();
    while (ss >> buf)  tokens.push_back(buf);
    if( tokens[0]=="L_BIN" ){
	values_effL.push_back( atof(tokens[9].c_str()) );
	values_sbL.push_back(  atof(tokens[3].c_str()) );
	values_muL.push_back(  atof(tokens[5].c_str()) );
	values_chiL.push_back( atof(tokens[7].c_str()) );
	values_binL.push_back( atof(tokens[1].c_str()) );
    }
    if( tokens[0]=="H_BIN" ){
	values_effH.push_back( atof(tokens[9].c_str()) );
	values_sbH.push_back(  atof(tokens[3].c_str()) );
	values_muH.push_back(  atof(tokens[5].c_str()) );
	values_chiH.push_back( atof(tokens[7].c_str()) );
	values_binH.push_back( atof(tokens[1].c_str()) );
    }

  }
  infile.close();
  int numBadFit_L(0), numBadFit_H(0);
  //Find the best values Low Eta
  float tmp1 = -1.;   int bin1 = -1.;
  float tmp2 = 100.;  int bin2 = -1.;
  for(unsigned int i(0); i<values_effL.size(); i++){
    //if(values_sbL[i] > tmp1 && values_chiL[i]<0.08 ){
    if(values_sbL[i] > tmp1 && values_chiL[i]<1. ){
	tmp1 = values_sbL[i];
	bin1 = i; cout<<"-->LOW_ETA LOW S/B bin:  "<<values_sbL[i]<<" bin "<<values_binL[bin1]<<endl;
    }
    //if(values_muL[i] < tmp2 && values_chiL[i]<0.08 ){
    if(values_muL[i] < tmp2 && values_chiL[i]<1. ){
	tmp2 = values_muL[i];
	bin2 = i; cout<<"-->LOW_ETA HIGH Mu bin:  "<<values_muL[i]<<" bin "<<values_binL[bin2]<<endl;
    }
    if( values_chiL[i]>0.08 ) numBadFit_L++;
  }
  //Find the best values High Eta
  tmp1 = -1.;   int bin3 = -1.;
  tmp2 = 100.;  int bin4 = -1.;
  for(unsigned int i(0); i<values_effH.size(); i++){
    if(values_sbH[i] > tmp1 && values_chiH[i]<0.019 && values_muH[i]>0.0002 ){
	tmp1 = values_sbH[i];
	bin3 = i; cout<<"-->HIGH_ETA LOW S/B bin:  "<<values_sbH[i]<<" bin "<<values_binH[bin3]<<endl;
    }
    if(values_muH[i] < tmp2 && values_chiH[i]<0.019 && values_muH[i]>0.00206 ){
	tmp2 = values_muH[i];
	bin4 = i; cout<<"-->HIGH_ETA HIGH Mu bin:  "<<values_muH[i]<<" bin "<<values_binH[bin4]<<endl;
    }
    if( values_chiH[i]>0.03 ) numBadFit_H++;
  }
  //Print them
  cout<<"LOW Eta"<<endl;
  cout<<"S/B) Bin: "<<values_binL[bin1]<<" sb: "<<values_sbL[bin1]<<" smu: "<<values_muL[bin1]<<" eff: "<<values_effL[bin1]<<" chi: "<<values_chiL[bin1]<<endl;
  cout<<"M/U) Bin: "<<values_binL[bin2]<<" sb: "<<values_sbL[bin2]<<" smu: "<<values_muL[bin2]<<" eff: "<<values_effL[bin2]<<" chi: "<<values_chiL[bin2]<<endl;
  cout<<"HIGH Eta"<<endl;
  cout<<"S/B) Bin: "<<values_binH[bin3]<<" sb: "<<values_sbH[bin3]<<" smu: "<<values_muH[bin3]<<" eff: "<<values_effH[bin3]<<" chi: "<<values_chiH[bin3]<<endl;
  cout<<"M/U) Bin: "<<values_binH[bin4]<<" sb: "<<values_sbH[bin4]<<" smu: "<<values_muH[bin4]<<" eff: "<<values_effH[bin4]<<" chi: "<<values_chiH[bin4]<<endl;
  cout<<"We have "<<numBadFit_L<<" bad fit in L and "<<numBadFit_L<<" in H."<<endl;
  cout<<endl;
  //Look for the cuts used
  ifstream infile2;
  infile2.open (Input2.Data());
  string Line2, check2;
  int nL2=1;
  while(!infile2.eof()){
    if(nL2==nLine2) break;
    nL2++;
    getline(infile2,Line);   
    string buf;
    stringstream ss(Line);
    vector<string> tokens; tokens.clear();
    while (ss >> buf)  tokens.push_back(buf);
    if( atof(tokens[1].c_str()) == values_binL[bin1]  ) cout<<"LOW: S/B Selection: "<<Line<<endl;
    if( atof(tokens[1].c_str()) == values_binL[bin2]  ) cout<<"LOW: SigmaMu/Mu Selection:"<<Line<<endl;
    if( atof(tokens[1].c_str()) == values_binH[bin3]  ) cout<<"HIGH: S/B Selection: "<<Line<<endl;
    if( atof(tokens[1].c_str()) == values_binH[bin4]  ) cout<<"HIGH: SigmaMu/Mu Selection:"<<Line<<endl;
  }
}
