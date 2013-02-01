// #define private public
// #define protected public

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TH1.h"
#include "RooPoisson.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooNDKeysPdf.h"
#include "TRandom.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"
#include "TBranchElement.h"
#include "TTreeFormula.h"
#include "TROOT.h"
#include "Cintex/Cintex.h"
#include "GBRForest.h"
#include <string.h>
#endif

//#include "TMVARegGui.C"

using namespace TMVA;
using namespace RooFit;



void applyenergy() {

  ROOT::Cintex::Cintex::Enable();   

  
  //printf("include: %s\n",gSystem->GetIncludePath());
  //return;
  
  Long64_t maxentries = -1;
  
  


  //TFile *fmc = new TFile("/home/bendavid/cms/hist/hgg-v0-Sept1/local/filefi/merged/hgg-v0_s11-h120gg-gf-v11-pu_noskim.root","READ");
  //TFile *fmc = new TFile("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_f11--h121gg-gf-v14b-pu_noskim.root","READ");
  TFile *fmc = new TFile("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_f11-zjets-v14b-pu_noskim.root","READ");  
  //TFile *fmc = new TFile("/scratch/bendavid/cms/hist/hgg-v0/t2mit/filefi/025/f11--h121gg-gf-v14b-pu/hgg-v0_f11--h121gg-gf-v14b-pu_noskim_0000.root","READ");
  TDirectory *dir = (TDirectory*)fmc->FindObjectAny("PhotonTreeWriterPresel");
  TTree *hmcph = (TTree*)dir->Get("hPhotonTree");

  TDirectory *dirsingle = (TDirectory*)fmc->FindObjectAny("PhotonTreeWriterSingle");
  TTree *hmcsingleph = (TTree*)dirsingle->Get("hPhotonTree");

  TFile *fmcele = new TFile("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_f11-zjets-v14b-pu_noskim.root","READ");  
  TDirectory *direle = (TDirectory*)fmcele->FindObjectAny("PhotonTreeWriterE");
  TTree *hmcele = (TTree*)direle->Get("hPhotonTree");
  TDirectory *direlesingle = (TDirectory*)fmcele->FindObjectAny("PhotonTreeWriterE");
  TTree *hmcelesingle = (TTree*)direlesingle->Get("hPhotonTreeSingle");

  
  TFile *fdele = new TFile("/scratch/bendavid/cms/hist/hgg-v0/MergedDel2011J16.root","READ");  
  TDirectory *dirdele = (TDirectory*)fdele->FindObjectAny("PhotonTreeWriterE");
  TTree *hdele = (TTree*)dirdele->Get("hPhotonTree");  
  TDirectory *dirdelesingle = (TDirectory*)fdele->FindObjectAny("PhotonTreeWriterE");
  TTree *hdelesingle = (TTree*)dirdelesingle->Get("hPhotonTreeSingle");  
  
  TFile *fgbropt = new TFile("fgbrtraintest.root","READ");
  const GBRForest *gbropt = (GBRForest*)fgbropt->Get("gbrtrain");      
  
  std::vector<std::string> *varlist = (std::vector<std::string>*)fgbropt->Get("varlist");

  
  std::vector<std::string> *varlisteb = varlist;
  std::vector<std::string> *varlistee = varlist;
  
  const GBRForest *gbr = 0;
  const GBRForest *gbreb = gbropt;
  const GBRForest *gbree = gbropt;
  
  UInt_t nvarseb = varlisteb->size();
  UInt_t nvarsee = varlistee->size();
  
  Float_t *vals = 0;
  
  Float_t *valseb = new Float_t[nvarseb];
  Float_t *valsee = new Float_t[nvarsee];
 
  
  TFile *fmvacor = new TFile("fmvacor.root","RECREATE");
  
  TTree *hmvacorph = new TTree("hmvacorph","");
  TTree *hmvacorele = new TTree("hmvacorele","");
  TTree *hmvacorelesingle = new TTree("hmvacorelesingle","");

  TTree *hmvacordele = new TTree("hmvacordele","");
  TTree *hmvacordelesingle = new TTree("hmvacordelesingle","");
  
  TTree *hmvacorqcdsingle = new TTree("hmvacorqcdsingle","");
  
  hmvacorele->SetAutoFlush(-1000000000);
  hmvacorelesingle->SetAutoFlush(-1000000000);

  hmvacordele->SetAutoFlush(-1000000000);
  hmvacordelesingle->SetAutoFlush(-1000000000);  
  
  
  Float_t massmvacor=0.;
  Float_t massmvacorerr=0.;
  Float_t massmvacorerrlo=0.;
  Float_t massmvacorerrhi=0.;
  Float_t ph1emvacor=0.;
  Float_t ph1emvacorerr=0.;
  Float_t ph1emvacorerrlo=0.;
  Float_t ph1emvacorerrhi=0.;
  Float_t ph1bdt = 0.;
  Float_t ph1bdtvar = 0.;
  Int_t   ph1dcoridx=0;
  Float_t ph1emvadcor=0.;
  Float_t ph2emvacor=0.;
  Float_t ph2emvacorerr=0.;
  Float_t ph2emvacorerrlo=0.;
  Float_t ph2emvacorerrhi=0.;
  Float_t ph2bdt = 0.;
  Float_t ph2bdtvar = 0.;
  Int_t   ph2dcoridx=0;
  Float_t ph2emvadcor=0.;  
  
  Float_t phemvacor=0.;
  Float_t phemvacorerr=0.;
  Float_t phbdt=0.;
  Float_t phbdtvar=0.;
  Int_t phdcoridx = 0;
  Float_t phemvadcor=0;
  
  Float_t phregtarget = 0.;
 
  for (UInt_t isample=1; isample<3; ++isample) {
    TTree *hmc = 0;
    TTree *hmvacor = 0;
    TTree *hmcsingle = 0;
    TTree *hmvacorsingle = 0;
    if (isample==0) {
      hmc = hmcph;
      hmvacor = hmvacorph;
      //hmcsingle = hmcqcdsingle;
      //hmvacorsingle = hmvacorqcdsingle;
    }
    else if (isample==1) {
      hmc = hmcele;
      hmvacor = hmvacorele;
      hmcsingle = hmcelesingle;
      hmvacorsingle = hmvacorelesingle;
    }    
    else if (isample==2) {
      hmc = hdele;
      hmvacor = hmvacordele;
      
      hmcsingle = hdelesingle;
      hmvacorsingle = hmvacordelesingle;
      
    } 

//     std::vector<TTreeFormula*> forms1;
//     std::vector<TTreeFormula*> forms2;
//     std::vector<TTreeFormula*> formssingle;
    TTreeFormula **formseb1 = new TTreeFormula*[nvarseb];
    TTreeFormula **formseb2 = new TTreeFormula*[nvarseb];
    TTreeFormula **formsebsingle = new TTreeFormula*[nvarseb];
    for (UInt_t ivar=0; ivar<varlisteb->size(); ++ivar) {
      TString expression = varlisteb->at(ivar);      
      expression.ReplaceAll("ph.genz","vtxZ");
      TString expr1(expression);
      expr1.ReplaceAll("ph.","ph1.");
      TString expr2(expression);
      expr2.ReplaceAll("ph.","ph2.");
      printf("expr = %s, expr1 = %s, expr2 = %s\n",expression.Data(),expr1.Data(),expr2.Data());
      formseb1[ivar] = new TTreeFormula(expr1,expr1,hmc);
      formseb2[ivar] = new TTreeFormula(expr2,expr2,hmc);
      if (hmcsingle) formsebsingle[ivar] = new TTreeFormula(expression,expression,hmcsingle);
    }

    TTreeFormula **formsee1 = new TTreeFormula*[nvarsee];
    TTreeFormula **formsee2 = new TTreeFormula*[nvarsee];
    TTreeFormula **formseesingle = new TTreeFormula*[nvarsee];
    for (UInt_t ivar=0; ivar<varlistee->size(); ++ivar) {
      TString expression = varlistee->at(ivar);    
      expression.ReplaceAll("ph.genz","vtxZ");      
      if (expression=="(1.0-(!ismc)*0.072)*ph.scpse/ph.scrawe") expression = "ph.scpse/ph.scrawe";
      TString expr1(expression);
      expr1.ReplaceAll("ph.","ph1.");
      TString expr2(expression);
      expr2.ReplaceAll("ph.","ph2.");
      printf("expr = %s, expr1 = %s, expr2 = %s\n",expression.Data(),expr1.Data(),expr2.Data());
      formsee1[ivar] = new TTreeFormula(expr1,expr1,hmc);
      formsee2[ivar] = new TTreeFormula(expr2,expr2,hmc);
      if (hmcsingle) formseesingle[ivar] = new TTreeFormula(expression,expression,hmcsingle);
    }
    
    TString denebexpr1 = "ph1.scrawe";
    TString denebexpr2 = "ph2.scrawe";
    TString denebexprsingle = "ph.scrawe";
    TTreeFormula *denebform1 = new TTreeFormula(denebexpr1,denebexpr1,hmc);
    TTreeFormula *denebform2 = new TTreeFormula(denebexpr2,denebexpr2,hmc);
    TTreeFormula *denebformsingle = 0;
    if (hmcsingle) denebformsingle = new TTreeFormula(denebexprsingle,denebexprsingle,hmcsingle);

//     TString deneeexpr1 = "ph1.scrawe + (1.0-(!ismc)*0.072)*ph1.scpse";
//     TString deneeexpr2 = "ph2.scrawe + (1.0-(!ismc)*0.072)*ph2.scpse";
   
    TString deneeexpr1 = "ph1.scrawe + ph1.scpse";
    TString deneeexpr2 = "ph2.scrawe + ph2.scpse";   
    
    TString deneeexprsingle = "ph.scrawe + ph.scpse";
    TTreeFormula *deneeform1 = new TTreeFormula(deneeexpr1,deneeexpr1,hmc);
    TTreeFormula *deneeform2 = new TTreeFormula(deneeexpr2,deneeexpr2,hmc);
    TTreeFormula *deneeformsingle = 0;
    if (hmcsingle) deneeformsingle = new TTreeFormula(deneeexprsingle,deneeexprsingle,hmcsingle);    
    
    TTreeFormula *costhetaform = new TTreeFormula("costheta","costheta",hmc);

    TString isbexpr1 = "ph1.isbarrel";
    TString isbexpr2 = "ph2.isbarrel";
    TString isbexprsingle = "ph.isbarrel";
    TTreeFormula *isbform1 = new TTreeFormula(isbexpr1,isbexpr1,hmc);
    TTreeFormula *isbform2 = new TTreeFormula(isbexpr2,isbexpr2,hmc);
    TTreeFormula *isbformsingle = 0;
    if (hmcsingle) isbformsingle = new TTreeFormula(isbexprsingle,isbexprsingle,hmcsingle);    
    
    hmvacor->Branch("massmvacor",&massmvacor,"massmvacor/F");
    hmvacor->Branch("massmvacorerr",&massmvacorerr,"massmvacorerr/F");
    //hmvacor->Branch("massmvacorerrlo",&massmvacorerrlo,"massmvacorerrlo/F");
    //hmvacor->Branch("massmvacorerrhi",&massmvacorerrhi,"massmvacorerrhi/F");
    
    hmvacor->Branch("ph1.emvacor",&ph1emvacor,"ph1.emvacor/F");
    hmvacor->Branch("ph1.emvacorerr",&ph1emvacorerr,"ph1.emvacorerr/F");
    hmvacor->Branch("ph1.bdt",&ph1bdt,"ph1.bdt/F");
    hmvacor->Branch("ph1.bdtvar",&ph1bdtvar,"ph1.bdtvar/F");    
    hmvacor->Branch("ph1.dcoridx",&ph1dcoridx,"ph1.dcoridx/I"); 
    hmvacor->Branch("ph1.emvadcor",&ph1emvadcor,"ph1.emvadcor/F");    
    
    //hmvacor->Branch("ph1.emvacorerrlo",&ph1emvacorerrlo,"ph1.emvacorerrlo/F");
    //hmvacor->Branch("ph1.emvacorerrhi",&ph1emvacorerrhi,"ph1.emvacorerrhi/F");
    hmvacor->Branch("ph2.emvacor",&ph2emvacor,"ph2.emvacor/F");
    hmvacor->Branch("ph2.emvacorerr",&ph2emvacorerr,"ph2.emvacorerr/F");  
    hmvacor->Branch("ph2.bdt",&ph2bdt,"ph2.bdt/F");
    hmvacor->Branch("ph2.bdtvar",&ph2bdtvar,"ph2.bdtvar/F");    
    hmvacor->Branch("ph2.dcoridx",&ph2dcoridx,"ph2.dcoridx/I");    
    hmvacor->Branch("ph2.emvadcor",&ph2emvadcor,"ph2.emvadcor/F");    
    
    
    //hmvacor->Branch("ph2.emvacorerrlo",&ph2emvacorerrlo,"ph2.emvacorerrlo/F");
    //hmvacor->Branch("ph2.emvacorerrhi",&ph2emvacorerrhi,"ph2.emvacorerrhi/F");
    
    if (hmvacorsingle) {
      hmvacorsingle->Branch("ph.emvacor",&phemvacor,"ph.emvacor/F");
      hmvacorsingle->Branch("ph.emvacorerr",&phemvacorerr,"ph.emvacorerr/F");
      hmvacorsingle->Branch("ph.dcoridx",&phdcoridx,"ph.dcoridx/I");    
      hmvacorsingle->Branch("ph.emvadcor",&phemvadcor,"ph.emvadcor/F");            
      //hmvacor->Branch("ph.bdt",&ph1bdt,"ph.bdt/F");
      //hmvacor->Branch("ph.bdtvar",&ph1bdtvar,"ph.bdtvar/F");      
    }
    
    //TString method = "MLP method";
    //TString method = "BDT method";
    TString method = "BDTG method";
    //TString method = "PDEFoam method";
    
    for (Long64_t i=0; i<hmc->GetEntries(); ++i) {
      hmc->LoadTree(i);

      float den1, den2;
      
      bool isb1 = isbform1->EvalInstance();
      bool isb2 = isbform2->EvalInstance();

      if (isb1) {
        gbr = gbreb;
        //gbrvar = gbrvareb;
        //gbrdcor = gbrdcoreb;
        vals = valseb;
        den1 = denebform1->EvalInstance();
        for (UInt_t ivar=0; ivar<nvarseb; ++ivar) {
          valseb[ivar] = formseb1[ivar]->EvalInstance();
        }        
      }
      else {
        gbr = gbree;
        //gbrvar = gbrvaree;
        //gbrdcor = gbrdcoree;
        vals = valsee;
        den1 = deneeform1->EvalInstance();
        for (UInt_t ivar=0; ivar<nvarsee; ++ivar) {
          valsee[ivar] = formsee1[ivar]->EvalInstance();
        }              
      }
      

      phregtarget = gbr->GetResponse(vals);     
      ph1emvacor = phregtarget*den1;

      //printf("phregtarget = %5f, ph1emvacor = %5f\n",phregtarget,ph1emvacor);
      
      
      if (isb2) {
        gbr = gbreb;
        vals = valseb;        
        den2 = denebform2->EvalInstance();
        for (UInt_t ivar=0; ivar<nvarseb; ++ivar) {
          valseb[ivar] = formseb2[ivar]->EvalInstance();
        }        
      }
      else {
        gbr = gbree;
        vals = valsee;        
        den2 = deneeform2->EvalInstance();
        for (UInt_t ivar=0; ivar<nvarsee; ++ivar) {
          valsee[ivar] = formsee2[ivar]->EvalInstance();
        }              
      }

      phregtarget = gbr->GetResponse(vals);      
      ph2emvacor = phregtarget*den2;

      

      
      massmvacor = TMath::Sqrt(2.0*ph1emvacor*ph2emvacor*(1.0-costhetaform->EvalInstance()));
      //massmvacorerr = 0.5*massmvacor*TMath::Sqrt(ph1emvacorerr*ph1emvacorerr/ph1emvacor/ph1emvacor + ph2emvacorerr*ph2emvacorerr/ph2emvacor/ph2emvacor);
      //massmvacorerrlo = 0.5*massmvacor*TMath::Sqrt(ph1emvacorerrlo*ph1emvacorerrlo/ph1emvacor/ph1emvacor + ph2emvacorerrlo*ph2emvacorerrlo/ph2emvacor/ph2emvacor);
      //massmvacorerrhi = 0.5*massmvacor*TMath::Sqrt(ph1emvacorerrhi*ph1emvacorerrhi/ph1emvacor/ph1emvacor + ph2emvacorerrhi*ph2emvacorerrhi/ph2emvacor/ph2emvacor);
        
      hmvacor->Fill();
      
    }
    hmc->AddFriend(hmvacor);
    hmvacor->Write();
    
    if (hmcsingle) {
      for (Long64_t i=0; i<hmcsingle->GetEntries(); ++i) {
        hmcsingle->LoadTree(i);

        float den;
        bool isbsingle = isbformsingle->EvalInstance();

        if (isbsingle) {
          gbr = gbreb;
          vals = valseb;          
          den = denebformsingle->EvalInstance();
          for (UInt_t ivar=0; ivar<nvarseb; ++ivar) {
            valseb[ivar] = formsebsingle[ivar]->EvalInstance();
          }        
        }
        else {
          gbr = gbree;
          vals = valsee;          
          den = deneeformsingle->EvalInstance();
          for (UInt_t ivar=0; ivar<nvarsee; ++ivar) {
            valsee[ivar] = formseesingle[ivar]->EvalInstance();
          }              
        }

        phregtarget = gbr->GetResponse(vals);     
        phemvacor = phregtarget*den;

        
        hmvacorsingle->Fill();
        
      }
      hmcsingle->AddFriend(hmvacorsingle);
      hmvacorsingle->Write();      
      
    }
    
  }
  
  
  
  
  
//   
}
