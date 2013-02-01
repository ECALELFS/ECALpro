// $Id: EcalEnerCorr.cc,v 1.0 2011/06/30 18:45:47 montanin

//#ifndef EcalEnerCorr_H
//#define EcalEnerCorr_H

#include "CalibCode/CalibTools/interface/EcalEnerCorr.h"

// Root includes
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"

#include <sstream>
#include <iostream>
#include <vector>
#include <map>

using namespace std;

//========================================================================================

EcalEnerCorr::EcalEnerCorr( ) {
//========================================================================================
  //if(!evt) throw std::invalid_argument("null pointer for fwlite::EventContainer. No data. no party!");
  cout << "building EcalEnerCorr ..." ;
 
   fillEnergyEdgeEB();
   fillEnergyEdgeEE();
   fillEtaEdgeEE();

   noCorrectionsEB_ = true;
   noCorrectionsEE_ = true;

   vdebug=false;
   EBedge=1.479;
   EEedge=3.0;
   
   extCorrPointEE = new TH1F*[ENBINSEE];

   cout << "... done. " << endl;
}
//============================================

void EcalEnerCorr::setVerbosity(bool n) {
   vdebug = n;
}
int  EcalEnerCorr::getEnBinsEB(){
   return ENBINSEB;
}
int  EcalEnerCorr::getEnBinsEE(){
   return ENBINSEE;
}
int  EcalEnerCorr::getEtaBinsEB(){
   return ETABINSEB;
}
int  EcalEnerCorr::getEtaBinsEE(){
   return ETABINSEE;
}
//============================================   

void EcalEnerCorr::fillEnergyEdgeEB() 
{
   enBinBoundEB[0] = 0. ;
   enBinBoundEB[1] = 0.7 ;
   enBinBoundEB[2] = 1.0 ;
   enBinBoundEB[3] = 1.1 ;
   enBinBoundEB[4] = 1.2 ; 
   enBinBoundEB[5] = 1.3 ;
   enBinBoundEB[6] = 1.5;
   enBinBoundEB[7] = 2.0;
   enBinBoundEB[8] = 2.5;
   enBinBoundEB[9] = 3.0;
   enBinBoundEB[10] = 3.5;
   enBinBoundEB[11] = 5.0;
   enBinBoundEB[12] = 10.0;
}

void EcalEnerCorr::fillEnergyEdgeEE() 
{
   enBinBoundEE[0] = 0. ;
   enBinBoundEE[1] = 1.9 ;
   enBinBoundEE[2] = 2.2 ;
   enBinBoundEE[3] = 2.5 ;
   enBinBoundEE[4] = 3.0 ;
   enBinBoundEE[5] = 3.5 ;
   enBinBoundEE[6] = 4.0 ;
   enBinBoundEE[7] = 4.5;
   enBinBoundEE[8] = 5.0;
   enBinBoundEE[9] = 5.5;
   enBinBoundEE[10] = 6.5;
   enBinBoundEE[11] = 7.5;
   enBinBoundEE[12] = 8.5;
   enBinBoundEE[13] = 10.0; 
   enBinBoundEE[14] = 15.0; 
   enBinBoundEE[15] = 30.0;
}

void EcalEnerCorr::fillEtaEdgeEE(){
   etaBinBoundEE[0] = 1.479 ; //****
   etaBinBoundEE[1] = 1.51 ;
   etaBinBoundEE[2] = 1.53 ;
   etaBinBoundEE[3] = 1.56 ;
   etaBinBoundEE[4] = 1.59 ;  //****
   etaBinBoundEE[5] = 1.61 ;
   etaBinBoundEE[6] = 1.64 ;
   etaBinBoundEE[7] = 1.67 ;  //****
   etaBinBoundEE[8] = 1.70 ;
   etaBinBoundEE[9] = 1.73 ;
   etaBinBoundEE[10]= 1.77 ;
   etaBinBoundEE[11]= 1.81 ;  //****
   etaBinBoundEE[12]= 1.85 ;
   etaBinBoundEE[13]= 1.89 ;
   etaBinBoundEE[14]= 1.93 ;
   etaBinBoundEE[15]= 1.97 ;  //****
   etaBinBoundEE[16]= 2.02 ;
   etaBinBoundEE[17]= 2.07 ;
   etaBinBoundEE[18]= 2.11 ;
   etaBinBoundEE[19]= 2.16 ;  //****
   etaBinBoundEE[20]= 2.21 ;
   etaBinBoundEE[21]= 2.26 ;
   etaBinBoundEE[22]= 2.31 ;
   etaBinBoundEE[23]= 2.37 ;  //****
   etaBinBoundEE[24]= 2.43 ;
   etaBinBoundEE[25]= 2.50 ;
   etaBinBoundEE[26]= 2.57 ;
   etaBinBoundEE[27]= 2.64 ;  //****
   etaBinBoundEE[28]= 2.74 ;
   etaBinBoundEE[29]= 2.86 ;
   etaBinBoundEE[30]= 3.0 ;   //****
}

//================================================================================================

// energy bins in EB which have an unique fit function
bool EcalEnerCorr::uniqueFunc(int i){
   bool n=false;
   if ( i==0 || i==1 || i==10 || i==11) n = true;
   return n;
}

// definition of the Super Module borders
bool EcalEnerCorr::etaBorderS(int i){
   
   bool cond=false;
   if (i==1 || i==20 || i==21 || i==40 || i==41 || i==60 || i==61 || i==85 || 
       i==87|| i==111|| i==112|| i==131|| i==132|| i==151|| i==152|| i==171){
      cond= true;
   }
   return cond;
}

// definition of the Module borders
bool EcalEnerCorr::etaBorderM(int i){

   bool cond=false;
   int mod = i%5;
   if ((i<86)&&(mod == 1 || mod ==0)&&(!etaBorderS(i))){cond= true;}
   if ((i>86)&&(mod == 1 || mod ==2)&&(!etaBorderS(i))){cond= true;}
   return cond;
}

//================================================================================================

double EcalEnerCorr::getContainmentCorrectionsEB( double energy, int ieta ) {
//================================================================================================
   // this work with eta index

   double maxEcluster = 10.0;
   
   if(energy>maxEcluster) return 1.;

   int myieta = TMath::Abs(ieta);
   double corrFactor=1.0;
   int ien;  
 
    //cout << "en,ieta="<< energy << "," << ieta << " :: ";

   /// find energy energy bin for the photon
   for(ien=0; ien< ENBINSEB; ++ien) { 
     if(energy <= enBinBoundEB[ien+1]) break;
   }

   /// SuperModule borders
   if (etaBorderS(myieta+86))
   {
      //cout << "S - I valori che uso: energia = "<< ien<<" - eta = "<<myieta<<endl;
      corrFactor = extCorrFunc1SupModBorders[ien]->Eval(myieta);
   }
   /// central region
   else if (myieta<58)
   {
      if (uniqueFunc(ien))
      {
	 
	 //cout << "1.U - I valori che uso: energia = "<< ien<<" - eta = "<<myieta<<endl;
         corrFactor = extCorrFunc1Bulk[ien]->Eval(myieta);
      }
      else if (etaBorderM(myieta+86))
      {
	 
	 //cout << "1.M - I valori che uso: energia = "<< ien<<" - eta = "<<myieta<<endl;
         corrFactor = extCorrFunc1ModBorders[ien]->Eval(myieta);
      }
      else 
      {
	 
	 //cout << "1.B - I valori che uso: energia = "<< ien<<" - eta = "<<myieta<<endl;
         corrFactor = extCorrFunc1Bulk[ien]->Eval(myieta);
      }
   
   }
   /// side region
   else if (myieta>57)
   {
      if (uniqueFunc(ien) )
      {
	 
	 //cout << "2.U - I valori che uso: energia = "<< ien<<" - eta = "<<myieta<<endl;
         corrFactor = extCorrFunc2Bulk[ien]->Eval(myieta);
      }
      else  if (etaBorderM(myieta+86))
      {
	 
	 //cout << "2.M - I valori che uso: energia = "<< ien<<" - eta = "<<myieta<<endl;
         corrFactor = extCorrFunc2ModBorders[ien]->Eval(myieta);
      }
      else
      { 
	 
	 //cout << "2.B - I valori che uso: energia = "<< ien<<" - eta = "<<myieta<<endl;
         corrFactor = extCorrFunc2Bulk[ien]->Eval(myieta);
      }
   } 

   //cout << endl;

   return 1./corrFactor;
}


//================================================================================================

double EcalEnerCorr::getContainmentCorrectionsEE( double energy, double ieta ) {
//================================================================================================
// this work with the eta value

   double maxEcluster = 30.0;
   
   if(energy>maxEcluster)
      return 1.;

   double myieta = TMath::Abs(ieta);
   double corrFactor=1.0;

   int ien;  
 
   /// find energy energy bin for the first photon
   for(ien=0; ien< ENBINSEE; ++ien) { 
     if(energy <= enBinBoundEE[ien+1]) break;
   }
   
   if (ien==ENBINSEE-1) corrFactor=0.95;
   else corrFactor = extCorrFunctionEE[ien]->Eval(myieta);
   
   return 1./corrFactor;
} 
//================================================================================================

double EcalEnerCorr::getContainmentPointCorrectionsEE( double energy, double ieta ) {
//================================================================================================
// this work with the eta value

   double maxEcluster = 30.0;
   double NetaBinsD=getEtaBinsEE();
   double deltaEta = (EEedge - EBedge)/NetaBinsD;

   if(energy>maxEcluster)
      return 1.;

   double myieta = TMath::Abs(ieta);
   int ietaBin=(myieta - EBedge)/deltaEta;
   ietaBin+=1;
   //cout<<"qua vediamo cosa abbiamo   "<<myieta<<"  "<<ietaBin<<"  "<<energy<<endl;
   double corrFactor=1.0; 
   int ien; 
 
   /// find energy energy bin for the first photon
   for(ien=0; ien< ENBINSEE; ++ien) { 
     if(energy <= enBinBoundEE[ien+1]) break;
   }

   //if(!extVCorrPointEE[ien]){cout<<" Qui c'e' qualcosa che non va nel load"<< endl;}
   //if(extCorrPointEE[ien]){cout<<" Il load in teoria funziona"<< endl;}

   //double nonso=extVCorrPointEE[ien]->GetEntries();
   //cout<< "fammi vedere il bin in eta  "<< ietaBin<< "   e il bin in energia  "<< ien<< endl;
   // TFile *fifi = new TFile("fifi.root","RECREATE");
//    for (int i=0; i<ENBINSEE;i++){
//       fifi->cd();
//       //extCorrPointEE[i]->Write();
//       extCorrPointEE[i]->Write();
//       //cout<< "vediamo il ciclo   "<< i<<"   "<<extVCorrPointEE[i]->GetBinContent(1)<<endl;
//    }
//    fifi->Close();

   if (ien==ENBINSEE-1) {corrFactor=0.95;}
   else if (extCorrPointEE[ien]->GetBinContent(ietaBin)!= 0.){
      corrFactor = extCorrPointEE[ien]->GetBinContent(ietaBin);
      //cout<<"Qualcosa ha letto"<<endl;
   }
   else {corrFactor = 1.;
      //cout<<"ha letto un bello zero"<<endl;
   }
   
   return 1./corrFactor;
} 

//================================================================================================

double EcalEnerCorr::getContainmentMinvCorrections( double ieta ) {
//================================================================================================
// this work with the eta value of the Pi0

   double myieta = TMath::Abs(ieta);
   double corrFactor=135.0;
   double NetaBins=10;
   double deltaEta;
   int ietaBin;

   if (myieta < EBedge)
   {
      NetaBins=extCorrMinvFunctionEB->GetNbinsX();
      deltaEta = EBedge/NetaBins;
      ietaBin=myieta/deltaEta +1.;
      corrFactor = extCorrMinvFunctionEB->GetBinContent(ietaBin);
   } 
   else if (myieta < EEedge)
   {
      NetaBins=extCorrMinvFunctionEE->GetNbinsX();
      deltaEta = (EEedge - EBedge)/NetaBins;
      ietaBin=(myieta - EBedge)/deltaEta + 1.;
      corrFactor = extCorrMinvFunctionEE->GetBinContent(ietaBin);
   }
   return 135./corrFactor;
}

//================================================================================================

bool EcalEnerCorr::loadContainmentCorrectionsEB(const char* cfile){
//================================================================================================
   TFile* corrFunctionInputFile = TFile::Open(cfile);
   cout << "loading corrections EB from <" << cfile << "> ..." << endl;

   for(int i=0;i<ENBINSEB;++i)
   {
      stringstream buffer_enBin;
	             buffer_enBin << i;
      string name;
             name ="f1s_corr_vs_eta_en_";
             name += buffer_enBin.str();
      extCorrFunc1SupModBorders[i] = (TF1*)corrFunctionInputFile->Get(name.c_str()); 
      if(!extCorrFunc1SupModBorders[i])
      {
         cout << "ERROR: Correction File loaded, but correction TF1 f1s not found: en = " <<i << endl;
         return false;
      }

      name="f1m_corr_vs_eta_en_";
      name += buffer_enBin.str();
      extCorrFunc1ModBorders[i] = (TF1*)corrFunctionInputFile->Get(name.c_str()); 
      if(!extCorrFunc1ModBorders[i])
      {
         cout << "ERROR: Correction File loaded, but correction TF1 f1m not found: en = " <<i << endl;
         return false;
      }

      name="f2m_corr_vs_eta_en_";
      name += buffer_enBin.str();
      extCorrFunc2ModBorders[i] = (TF1*)corrFunctionInputFile->Get(name.c_str()); 
      if(!extCorrFunc2ModBorders[i])
      {
         cout << "ERROR: Correction File loaded, but correction TF1 f2m not found: en = " <<i << endl;
         return false;
      }

      if (uniqueFunc(i))
       {
          name="f1c_corr_vs_eta_en_";
          name += buffer_enBin.str();
          extCorrFunc1Bulk[i] = (TF1*)corrFunctionInputFile->Get(name.c_str()); 
          if(!extCorrFunc1Bulk[i])
          {
             cout << "ERROR: Correction File loaded, but correction TF1 f1c not found: en = " <<i << endl;
             return false;
          }
          
          name="f2c_corr_vs_eta_en_";
          name += buffer_enBin.str();
          extCorrFunc2Bulk[i] = (TF1*)corrFunctionInputFile->Get(name.c_str()); 
          if(!extCorrFunc2Bulk[i])
          {
             cout << "ERROR: Correction File loaded, but correction TF1 f2c not found: en = " <<i << endl;
             return false;
          }

       }
       else 
       {
          name="f1o_corr_vs_eta_en_";
          name += buffer_enBin.str();
          extCorrFunc1Bulk[i] = (TF1*)corrFunctionInputFile->Get(name.c_str()); 
          if(!extCorrFunc1Bulk[i])
          {
             cout << "ERROR: Correction File loaded, but correction TF1 f1o not found: en = " <<i << endl;
             return false;
          }
        
          name="f2o_corr_vs_eta_en_";
          name += buffer_enBin.str();
          extCorrFunc2Bulk[i] = (TF1*)corrFunctionInputFile->Get(name.c_str()); 
          if(!extCorrFunc2Bulk[i])
          {
             cout << "ERROR: Correction File loaded, but correction TF1 f2o not found: en = " <<i << endl;
             return false;
          }
       }
     
//        name="corrSM_vs_eta_en_";
//        name += buffer_enBin.str();
//        extCorrectionPointsS[i] = (TH1F*)corrFunctionInputFile->Get(name.c_str()); 
//        if(!extCorrectionPointsS[i])
//        {
//           cout << "ERROR: Correction File loaded, but correction TH1F for SM not found: en = " <<i << endl;
//           return false;
//        }
  
       // name="f1m_corr_vs_eta_en_";
//        name += buffer_enBin.str();
//        extCorrectionFunctionMinv[i] = (TF1*)corrFunctionInputFile->Get(name.c_str()); 
//        if(!extCorrectionFunctionMinv[i])
//        {
//           cout << "ERROR: Correction File loaded, but correction TF1 for Minv not found: en = " <<i << endl;
//           return false;
//        }
   }

   corrFunctionInputFile->Close();
   noCorrectionsEB_ = false;
   cout << "done." << endl;

   return true;
}

//================================================================================================

bool EcalEnerCorr::loadContainmentCorrectionsEE(const char* cfile){
//================================================================================================

   TFile* corrFunctionInputFile = TFile::Open(cfile);
   cout << "loading corrections EE from <" << cfile << "> ..." << endl;

   for(int i=0;i<ENBINSEE;++i)
   {
      stringstream buffer_enBin;
	             buffer_enBin << i;
      string name;
             name ="f_corr_vs_eta_en_";
             name += buffer_enBin.str();
      extCorrFunctionEE[i] = (TF1*)corrFunctionInputFile->Get(name.c_str()); 
      if(!extCorrFunctionEE[i])
      {
         cout << "ERROR: Correction File loaded, but correction TF1 f in EE not found: en = " <<i << endl;
         return false;
      }
   }

   corrFunctionInputFile->Close();
   noCorrectionsEE_ = false;
   cout << "done." << endl;

   return true;
}

//================================================================================================

bool EcalEnerCorr::loadContainmentPointCorrectionsEE(const char* cfile){
//================================================================================================

   TFile* corrFunctionInputFile = TFile::Open(cfile);
   cout << "loading point corrections EE from <" << cfile << "> ..." << endl;

   for(int i=0;i<ENBINSEE;++i)
   {
      stringstream bufferBin;
      bufferBin << i;
      string name;
      name ="corr_vs_eta_en_";
      name += bufferBin.str();
      extCorrPointEE[i] = (TH1F*)corrFunctionInputFile->Get(name.c_str());
      extCorrPointEE[i]->SetDirectory(0);
      if(!extCorrPointEE[i])
      {
         cout << "ERROR: Correction File loaded, but correction TH1F h in EE not found: en = " <<i << endl;
         return false;
      }
   }

   corrFunctionInputFile->Close();
   noCorrectionsEE_ = false;
   cout << "done." << endl;

   return true;
}


//================================================================================================
bool EcalEnerCorr::loadContainmentMinvCorrections(const char* cfile){
//================================================================================================

   TFile* corrMinvInputFile = TFile::Open(cfile);
   cout << "loading Minv corrections from <" << cfile << "> ..." << endl;
    
   extCorrMinvFunctionEB = (TH1F*)corrMinvInputFile->Get("funCorrPi0Mass_reco_vs_etaRec_EB_mu"); 
   extCorrMinvFunctionEB->SetDirectory(0);
   if(!extCorrMinvFunctionEB)
   {
      cout << "ERROR: Correction File loaded, but correction TH1F Minv in EB not found"<<endl;
      return false;
   }
   extCorrMinvFunctionEE = (TH1F*)corrMinvInputFile->Get("funCorrPi0Mass_reco_vs_etaRec_EE_mu"); 
   extCorrMinvFunctionEE->SetDirectory(0);
   if(!extCorrMinvFunctionEE)
   {
      cout << "ERROR: Correction File loaded, but correction TH1F Minv in EE not found"<<endl;
      return false;
   }

   corrMinvInputFile->Close();
   cout << "done." << endl;

   return true;
}

//#endif
