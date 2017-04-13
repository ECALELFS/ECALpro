#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>     //to use ostringstream to convert numbers to string in c++
#include <vector>
#include <cstring>
#include <cmath>
#include <math.h>

#include <TTree.h>
#include <TChain.h>
#include <TROOT.h>
#include <TFile.h>

#include "calibAnaEcalEB.h"
#include "calibAnaEcalEE.h"

using namespace std;

Int_t main(int argc, char* argv[]) {

  string fileToChain = "";
  string path(argv[1]);
  string eosPath = "root://eoscms//eos/cms" + path;
  string dirName(argv[2]);
  string iter_num(argv[3]);
  string tagName(argv[4]);
  string EBorEE(argv[5]);
  string wwwPath(argv[6]);
  string Pi0orEta(argv[7]);

  TChain *chain;
  string Result;          // string which will contain the number used as the loop iterator below (to create the chian of files)

  string iterNumber = "iter_" + iter_num;
  calibAnaEcal *ana = NULL;

  if (EBorEE == "EB") {
    
    chain = new TChain("calibEB");

    for (Int_t i = 0; i < 31; i++) {
      ostringstream convert;   // stream used for the conversion                                                                                                       
      convert << i;      // insert the textual representation of 'i' in the characters in the stream                     
      Result = convert.str();
      fileToChain = eosPath + dirName + "/iter_" + iter_num + "/" + tagName + "Barrel_" + Result + "_calibMap.root";
      chain->Add(TString(fileToChain.c_str()));
    }

    ana = new calibAnaEcalEB(chain);

  } else {

    chain = new TChain("calibEE");

    for (Int_t i = 0; i < 8; i++) {
      ostringstream convert;   // stream used for the conversion                            
      convert << i;      // insert the textual representation of 'i' in the characters in the stream                  
      Result = convert.str();
      fileToChain = eosPath + dirName + "/iter_" + iter_num + "/" + tagName + "Endcap_" + Result + "_calibMap.root";
      chain->Add(TString(fileToChain.c_str()));
    }

    ana = new calibAnaEcalEE(chain);

  }

  if (ana != NULL) {
    ana->setEBorEE(EBorEE);
    ana->setPi0orEta(Pi0orEta);
    //  ana->setWhichEE("");
    ana->setDirName(dirName);
    ana->setIterNumber(iterNumber);
    ana->setWwwPath(wwwPath);
    ana->Loop();
  }    

  return 0;


}
