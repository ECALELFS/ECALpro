#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TClass.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVirtualFitter.h>

#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h                                                                                  
#include <cstdio>
#include <cmath>
#include <sstream> //to use ostringstream to convert numbers to string in c++                                                          
 
using namespace std;

void countHistoInFile(const string& fname = "", const string& folder = "", const string& histoNameMatch = "") { 

  TKey *key; 
  TFile *f = TFile::Open(fname.c_str(), "READ"); 
  if (!f || f->IsZombie()) { 
    cout << "Unable to open " << fname << " for reading..." <<endl; 
    return; 
  } 
  
  TDirectory *dir = NULL;
  if (folder != "") {
    if (f->GetKey(folder.c_str())) {
      dir = f->GetDirectory(folder.c_str());
      dir->cd();
    } else {
      cout << "Error: could not find folder " << folder << ". Exit" << endl;
      exit(EXIT_FAILURE);
    }
  } else dir = f;

  Int_t total = 0; 
  TIter next((TList *) dir->GetListOfKeys()); 

  while ((key = (TKey *)next())) { 
    TClass *cl = gROOT->GetClass(key->GetClassName()); 
    if (cl->InheritsFrom("TH1")) { 
      // the following line is not needed if you only want 
      // to count the histograms 
      TH1 *h = (TH1 *)key->ReadObj(); 
      string hname(h->GetName());
      if (histoNameMatch != "" and hname.find(histoNameMatch) != string::npos) {
	//cout << "Histo found: " << h->GetName() << " - " << h->GetTitle() << endl; 
	total++; 
      }
    } 
  } 
  cout << "Found " << total << " Histograms" << endl; 

}
