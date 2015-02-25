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
  vector< float > values_eff; values_eff.clear();
  vector< float > values_sb; values_sb.clear();
  vector< float > values_mu; values_mu.clear();
  vector< float > values_chi; values_chi.clear();
  vector< int > values_bin; values_bin.clear();
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
    values_eff.push_back( atof(tokens[9].c_str()) );
    values_sb.push_back(  atof(tokens[3].c_str()) );
    values_mu.push_back(  atof(tokens[5].c_str()) );
    values_chi.push_back( atof(tokens[7].c_str()) );
    values_bin.push_back( atof(tokens[1].c_str()) );
  }
  infile.close();
  //Find the best values
  float tmp1 = -1.;   int bin1 = -1.;
  float tmp2 = 100.;  int bin2 = -1.;
  for(unsigned int i(0); i<values_eff.size(); i++){
    if(values_sb[i] > tmp1){
	tmp1 = values_sb[i];
	bin1 = i; cout<<"--> S/B bin:  "<<values_sb[i]<<" bin "<<values_bin[bin1]<<endl;
    }
    if(values_mu[i] < tmp2){
	tmp2 = values_mu[i];
	bin2 = i; cout<<"--> Mu bin:  "<<values_mu[i]<<" bin "<<values_bin[bin2]<<endl;
    }
  }
  //Print them
  cout<<"S/B) Bin: "<<values_bin[bin1]<<" sb: "<<values_sb[bin1]<<" smu: "<<values_mu[bin1]<<" eff: "<<values_eff[bin1]<<" chi: "<<values_chi[bin1]<<endl;
  cout<<"M/U) Bin: "<<values_bin[bin2]<<" sb: "<<values_sb[bin2]<<" smu: "<<values_mu[bin2]<<" eff: "<<values_eff[bin2]<<" chi: "<<values_chi[bin2]<<endl;
  //Look for the cuts used
  ifstream infile2;
  infile2.open (Input2.Data());
  string Line2, check2;
  int nL2=1;
  while(!infile2.eof()){
    if(nL2==nLine) break;
    nL2++;
    getline(infile2,Line);   
    string buf;
    stringstream ss(Line);
    vector<string> tokens; tokens.clear();
    while (ss >> buf)  tokens.push_back(buf);
    if( atof(tokens[1].c_str()) == values_bin[bin1]  ) cout<<"S/B Selection: "<<Line<<endl;
    if( atof(tokens[1].c_str()) == values_bin[bin2]  ) cout<<"SigmaMu/Mu Selection:"<<Line<<endl;
  }
}
