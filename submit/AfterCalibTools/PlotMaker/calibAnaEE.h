//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun  1 01:49:10 2016 by ROOT version 6.06/01
// from TTree calibEE/Tree of EE Inter-calibration constants
// found on file: pi0data_2016B_Run273730Barrel_0_calibMap.root
//////////////////////////////////////////////////////////

#ifndef calibAnaEE_h
#define calibAnaEE_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <string>

// Header file for the classes stored in the TTree if any.

class calibAnaEE {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           ix;
   Int_t           iy;
   Int_t           zside;
   Int_t           sc;
   Int_t           isc;
   Int_t           ic;
   Int_t           iquadrant;
   Int_t           hashedIndex;
   Int_t           iter;
   Float_t         coeff;
   Float_t         Signal;
   Float_t         Backgr;
   Float_t         Chisqu;
   Float_t         Ndof;
   Float_t         fit_mean;
   Float_t         fit_mean_err;
   Float_t         fit_sigma;
   Float_t         fit_Snorm;
   Float_t         fit_b0;
   Float_t         fit_b1;
   Float_t         fit_b2;
   Float_t         fit_b3;
   Float_t         fit_Bnorm;

   // List of branches
   TBranch        *b_ix;   //!
   TBranch        *b_iy;   //!
   TBranch        *b_zside;   //!
   TBranch        *b_sc;   //!
   TBranch        *b_isc;   //!
   TBranch        *b_ic;   //!
   TBranch        *b_iquadrant;   //!
   TBranch        *b_hashedIndex;   //!
   TBranch        *b_iter;   //!
   TBranch        *b_coeff;   //!
   TBranch        *b_Signal;   //!
   TBranch        *b_Backgr;   //!
   TBranch        *b_Chisqu;   //!
   TBranch        *b_Ndof;   //!
   TBranch        *b_fit_mean;   //!
   TBranch        *b_fit_mean_err;   //!
   TBranch        *b_fit_sigma;   //!
   TBranch        *b_fit_Snorm;   //!
   TBranch        *b_fit_b0;   //!
   TBranch        *b_fit_b1;   //!
   TBranch        *b_fit_b2;   //!
   TBranch        *b_fit_b3;   //!
   TBranch        *b_fit_Bnorm;   //!

   calibAnaEE(TTree *tree=0);
   virtual ~calibAnaEE();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string whichEE, std::string dirName, std::string iterNumber);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef calibAnaEE_cxx
calibAnaEE::calibAnaEE(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pi0data_2016B_Run273730Barrel_0_calibMap.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("pi0data_2016B_Run273730Barrel_0_calibMap.root");
      }
      f->GetObject("calibEE",tree);

   }
   Init(tree);
}

calibAnaEE::~calibAnaEE()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t calibAnaEE::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t calibAnaEE::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void calibAnaEE::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ix", &ix, &b_ix);
   fChain->SetBranchAddress("iy", &iy, &b_iy);
   fChain->SetBranchAddress("zside", &zside, &b_zside);
   fChain->SetBranchAddress("sc", &sc, &b_sc);
   fChain->SetBranchAddress("isc", &isc, &b_isc);
   fChain->SetBranchAddress("ic", &ic, &b_ic);
   fChain->SetBranchAddress("iquadrant", &iquadrant, &b_iquadrant);
   fChain->SetBranchAddress("hashedIndex", &hashedIndex, &b_hashedIndex);
   fChain->SetBranchAddress("iter", &iter, &b_iter);
   fChain->SetBranchAddress("coeff", &coeff, &b_coeff);
   fChain->SetBranchAddress("Signal", &Signal, &b_Signal);
   fChain->SetBranchAddress("Backgr", &Backgr, &b_Backgr);
   fChain->SetBranchAddress("Chisqu", &Chisqu, &b_Chisqu);
   fChain->SetBranchAddress("Ndof", &Ndof, &b_Ndof);
   fChain->SetBranchAddress("fit_mean", &fit_mean, &b_fit_mean);
   fChain->SetBranchAddress("fit_mean_err", &fit_mean_err, &b_fit_mean_err);
   fChain->SetBranchAddress("fit_sigma", &fit_sigma, &b_fit_sigma);
   fChain->SetBranchAddress("fit_Snorm", &fit_Snorm, &b_fit_Snorm);
   fChain->SetBranchAddress("fit_b0", &fit_b0, &b_fit_b0);
   fChain->SetBranchAddress("fit_b1", &fit_b1, &b_fit_b1);
   fChain->SetBranchAddress("fit_b2", &fit_b2, &b_fit_b2);
   fChain->SetBranchAddress("fit_b3", &fit_b3, &b_fit_b3);
   fChain->SetBranchAddress("fit_Bnorm", &fit_Bnorm, &b_fit_Bnorm);
   Notify();
}

Bool_t calibAnaEE::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void calibAnaEE::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t calibAnaEE::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef calibAnaEE_cxx
