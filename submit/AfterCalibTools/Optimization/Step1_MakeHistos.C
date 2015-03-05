#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TFormula.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMarker.h>
#include <TChain.h>
#include <iostream>
#include <memory>
#include <string>
#include <map>
#include <vector>
#include "TTree.h"
#include "TLatex.h"
#include "TMath.h"
#include "TBranch.h"
#include "TFile.h"
#include "TStyle.h"
#include "TString.h"

#define NPI0MAX 15000

using namespace std;

//.x Step1_MakeHistos.C+("ALL_MINBIAS_UNCAL_L1_NOL1FILTER", true, true)
void Step1_MakeHistos( TString folder = "ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns", bool isEB=true, bool isPi0=false ){

  //Check options
  if(isEB)  cout<<"Running on Barrel"<<endl;
  else      cout<<"Running on Endcap"<<endl;
  if(isPi0) cout<<"Running on Pi0"<<endl;
  else      cout<<"Running on Eta"<<endl;
  //OutPut Files
  TString Comm = "mkdir -p " + folder;
  system( Comm.Data() );
  FILE *file_txt;
  TString name_txt;
  if(isEB){
    name_txt = "/Fstep1_BINcut_EB_pi0.txt";
    if(!isPi0) name_txt = "/Fstep1_BINcut_EB_eta.txt";
  }
  else{
    name_txt = "/Fstep1_BINcut_EE_pi0.txt";
    if(!isPi0) name_txt = "/Fstep1_BINcut_EE_eta.txt";
  }
  file_txt=fopen( (folder + name_txt).Data(),"w");
  TString nameOutput;
  if(isEB){
    nameOutput="/Fstep1_EB_pi0.root";
    if(!isPi0) nameOutput="/Fstep1_EB_eta.root";
  }
  else{
    nameOutput="/Fstep1_EE_pi0.root";
    if(!isPi0) nameOutput="/Fstep1_EE_eta.root";
  }
  TFile* OutFile = new TFile( (folder + nameOutput).Data(),"RECREATE");
  //Input File
  TChain *tree = new TChain("Tree_Optim","Tree_Optim");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_0.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_1.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_10.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_100.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_101.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_102.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_103.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_105.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_106.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_107.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_108.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_109.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_11.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_110.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_111.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_112.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_113.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_114.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_115.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_116.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_117.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_118.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_119.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_12.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_120.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_121.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_122.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_123.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_124.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_125.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_126.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_127.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_128.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_129.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_13.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_130.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_131.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_132.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_133.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_134.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_135.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_136.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_137.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_138.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_139.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_14.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_140.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_141.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_142.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_144.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_145.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_146.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_147.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_148.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_149.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_15.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_150.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_151.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_153.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_154.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_155.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_156.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_157.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_158.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_159.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_16.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_160.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_161.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_162.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_163.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_164.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_165.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_166.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_167.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_168.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_169.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_17.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_170.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_171.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_172.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_173.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_174.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_175.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_176.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_177.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_178.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_179.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_18.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_180.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_181.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_182.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_183.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_184.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_186.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_187.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_189.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_190.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_193.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_194.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_195.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_196.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_197.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_198.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_199.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_2.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_20.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_200.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_201.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_203.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_204.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_205.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_206.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_207.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_209.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_21.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_211.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_212.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_214.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_215.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_216.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_217.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_218.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_219.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_22.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_220.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_221.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_222.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_223.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_224.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_225.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_227.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_23.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_230.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_234.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_236.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_237.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_239.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_24.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_241.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_242.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_244.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_245.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_246.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_247.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_248.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_249.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_25.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_250.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_251.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_252.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_253.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_254.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_255.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_256.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_257.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_258.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_259.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_26.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_260.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_264.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_268.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_269.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_27.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_273.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_274.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_275.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_279.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_28.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_280.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_282.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_283.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_284.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_286.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_287.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_288.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_289.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_29.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_290.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_291.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_292.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_294.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_295.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_296.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_297.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_298.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_299.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_3.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_30.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_301.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_302.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_304.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_308.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_31.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_313.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_314.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_318.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_32.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_320.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_321.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_322.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_324.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_325.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_326.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_327.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_328.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_330.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_331.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_332.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_333.root");
  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_334.root");
  //tree->Add("root://eoscms//eos/cms/store/group/dpg_ecal/alca_ecalcalib/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EE_pi0/iter_0/epsilonPlots.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_0.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_1.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_10.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_100.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_101.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_102.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_103.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_105.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_106.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_107.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_108.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_109.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_11.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_110.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_111.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_112.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_113.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_114.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_115.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_116.root");
  //tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU25ns_EB_pi0/iter_0/EcalNtp_117.root");

//  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU50ns_EB_pi0_NewSeed/iter_0/EcalNtp_0.root");
//  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU50ns_EB_pi0_NewSeed/iter_0/EcalNtp_1.root");
//  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU50ns_EB_pi0_NewSeed/iter_0/EcalNtp_10.root");
//  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU50ns_EB_pi0_NewSeed/iter_0/EcalNtp_11.root");
//  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU50ns_EB_pi0_NewSeed/iter_0/EcalNtp_12.root");
//  tree->Add("root://eoscms//eos/cms/store/caf/user/lpernie/ALL_MINBIAS_UNCAL_L1_NOL1FILTER_40PU50ns_EB_pi0_NewSeed/iter_0/EcalNtp_13.root");
  Int_t event = tree->GetEntries();
  cout << "Number of events in tree: " << event << endl;
  //Check TTree
  Int_t currentTreeId = -1;
  Int_t oldTree = -1;
  for(Int_t ievent = 0; ievent < event; ievent++) {
    tree->LoadTree(ievent);
    currentTreeId = tree->GetTreeNumber();
    if(currentTreeId != oldTree) {
	cout << "Checking TTree: " << currentTreeId << endl;
	oldTree = currentTreeId;
    }
  }
  //Variables declaration 
  Int_t npi;
  tree->SetBranchAddress("STr2_NPi0_rec",&npi);
  Float_t  massPi0[NPI0MAX];
  tree->SetBranchAddress("STr2_mPi0_rec",&massPi0);
  Int_t    iseb[NPI0MAX];
  tree->SetBranchAddress("STr2_Pi0recIsEB",&iseb); // 1=barrrel, 2=endcap
  Float_t  ptpi0[NPI0MAX];
  tree->SetBranchAddress("STr2_ptPi0_rec",&ptpi0);
  Float_t  etapi0[NPI0MAX];
  tree->SetBranchAddress("STr2_etaPi0_rec",&etapi0);
  Float_t ptclu1[NPI0MAX];
  tree->SetBranchAddress("STr2_ptG1_rec",&ptclu1);
  Float_t ptclu2[NPI0MAX];
  tree->SetBranchAddress("STr2_ptG2_rec",&ptclu2);
  Int_t ncris1[NPI0MAX];
  tree->SetBranchAddress("STr2_n1CrisPi0_rec",&ncris1);
  Int_t ncris2[NPI0MAX];
  tree->SetBranchAddress("STr2_n2CrisPi0_rec",&ncris2);
  Float_t E_Es_e1_1[NPI0MAX];
  tree->SetBranchAddress("STr2_Es_e1_1",&E_Es_e1_1);
  Float_t E_Es_e2_1[NPI0MAX];
  tree->SetBranchAddress("STr2_Es_e2_1",&E_Es_e2_1);
  Float_t E_Es_e1_2[NPI0MAX];
  tree->SetBranchAddress("STr2_Es_e1_2",&E_Es_e1_2);
  Float_t E_Es_e2_2[NPI0MAX];
  tree->SetBranchAddress("STr2_Es_e2_2",&E_Es_e2_2);
  Float_t iso[NPI0MAX];
  tree->SetBranchAddress("STr2_IsoPi0_rec",&iso);
  Float_t s4s9_1[NPI0MAX];
  tree->SetBranchAddress("STr2_S4S9_1",&s4s9_1);
  Float_t s4s9_2[NPI0MAX];
  tree->SetBranchAddress("STr2_S4S9_2",&s4s9_2);
  //Range of the selection
  //Nxtal
  vector <Int_t> ncri1cut;
  for(Int_t j=4;j<9;j++) ncri1cut.push_back(j);
  vector <Int_t> ncri2cut;
  for(Int_t j=4;j<7;j++) ncri2cut.push_back(j);
  //Pt_clus
  Double_t pcs=isPi0?0.4:1.4;
  Double_t pcf=isPi0?1.0:2.0;
  Double_t pcp=0.2;
  vector <Float_t> ptclucut;
  for(Double_t j=pcs;j<=pcf;j+=pcp) ptclucut.push_back(j);
  //Pt Pi0
  Double_t pps=isPi0?1.2:3.0;
  Double_t ppf=isPi0?2.4:4.0;
  Double_t ppp=0.3;
  vector <Float_t> ptPi0cut;
  for(Double_t j=pps;j<=ppf;j+=ppp) ptPi0cut.push_back(j);
  //ES --> not used for now
  Double_t ess=0.0;
  Double_t esf=1.5;
  Double_t esp=0.3;
  vector <Float_t> elayercut; elayercut.clear();
  for(Double_t j=ess;j<=esf;j+=esp) elayercut.push_back(j);
  //S4S9
  Double_t s4i=isPi0?0.7:0.75;
  Double_t s4f=isPi0?0.9:0.95;
  Double_t s4p=0.05;
  vector <Float_t> s4s9cut;
  for(Double_t j=s4i;j<=s4f;j+=s4p) s4s9cut.push_back(j);
  //ISO
  Double_t isoi=0.0;
  Double_t isof=0.3;
  Double_t isop=0.05;
  vector <Float_t> isocut;
  for(Double_t j=isoi;j<=isof;j+=isop) isocut.push_back(j);
  //Print txt with Cut id and cits used
  int nBin(0);
  for(unsigned i=0; i<ncri1cut.size();i++ ){
    for(unsigned j=0; j<ncri2cut.size();j++ ){
	for(unsigned k=0; k<ptclucut.size();k++ ){
	  //for(unsigned h=0; h<elayercut.size();h++ ){
	    for(unsigned s=0; s<s4s9cut.size();s++ ){
		for(unsigned q=0; q<isocut.size();q++ ){
		  for(unsigned z=0; z<ptPi0cut.size();z++ ){
		    //fprintf(file_txt,"BIN=%i  ncri1cut %i  ncri2cut %i  ptclucut %.4f  elayercut %.4f s4s9 %.4f Iso %.4f PtPi0 %.4f \n", nBin, ncri1cut[i],  ncri2cut[j], ptclucut[k], elayercut[h], s4s9cut[s], isocut[q], ptPi0cut[z] );
		    fprintf(file_txt,"BIN %i ncri1cut %i ncri2cut %i ptclucut %.2f s4s9 %.2f Iso %.2f PtPi0 %.2f \n", nBin, ncri1cut[i],  ncri2cut[j], ptclucut[k], s4s9cut[s], isocut[q], ptPi0cut[z] );
		    nBin++;
		  }
		}
	    }
	  //}
	}
    }
  }
  cout<<"You are optimizing for: "<<nBin<<" bins."<<endl;
  float xmin=0.05, xmax=0.3;
  if(!isPi0){ xmin=0.2; xmax=1.;}
  TH2F *hmass        = new TH2F("hmass","Histo For optimization",nBin,0.5,nBin+0.5,100,xmin,xmax);
  TH2F *hmassEBL     = new TH2F("hmassEBL","Histo For optimization",nBin,0.5,nBin+0.5,100,xmin,xmax);
  TH2F *hmassEBH     = new TH2F("hmassEBH","Histo For optimization",nBin,0.5,nBin+0.5,100,xmin,xmax);
  TH2F *hmassEEL     = new TH2F("hmassEEL","Histo For optimization",nBin,0.5,nBin+0.5,100,xmin,xmax);
  TH2F *hmassEEH     = new TH2F("hmassEEH","Histo For optimization",nBin,0.5,nBin+0.5,100,xmin,xmax);
  TH2F *hmass_tot    = new TH2F("hmass_tot","Histo For optimization",nBin,0.5,nBin+0.5,100,xmin,xmax);
  TH2F *hmass_totEBL = new TH2F("hmass_totEBL","Histo For optimization",nBin,0.5,nBin+0.5,100,xmin,xmax);
  TH2F *hmass_totEBH = new TH2F("hmass_totEBH","Histo For optimization",nBin,0.5,nBin+0.5,100,xmin,xmax);
  TH2F *hmass_totEEL = new TH2F("hmass_totEEL","Histo For optimization",nBin,0.5,nBin+0.5,100,xmin,xmax);
  TH2F *hmass_totEEH = new TH2F("hmass_totEEH","Histo For optimization",nBin,0.5,nBin+0.5,100,xmin,xmax);

  int EBorEE(-1);
  if(isEB)  EBorEE=1;
  else      EBorEE=2;
  cout << "Start reading data ..." << endl;
  for(int ii=0; ii<event; ii++){
    tree->GetEntry(ii);
    if( ii%100000==0 ) cout<<"Ev: "<<ii<<" / "<<event<<endl;
    for(Int_t jj=0; jj<npi; jj++){
	if( fabs(etapi0[jj]) < 2.5 ){ //I only optimize at low #eta (this cut has effect only in EE)
	  //Start Grid of cuts
	  Int_t binTrue=0;
	  for(unsigned i=0; i<ncri1cut.size();i++ ){
	    for(unsigned j=0; j<ncri2cut.size();j++ ){
		for(unsigned k=0; k<ptclucut.size();k++ ){
		  //for(unsigned h=0; h<elayercut.size();h++ ){
		    for(unsigned s=0; s<s4s9cut.size();s++ ){
			for(unsigned q=0; q<isocut.size();q++ ){
			  for(unsigned z=0; z<ptPi0cut.size();z++ ){ 
			    //Filling Histos
//			    if( iseb[jj]==EBorEE && ncris1[jj]>ncri1cut[i] && ncris2[jj]>ncri2cut[j] && ptclu1[jj]>ptclucut[k] && ptclu2[jj]>ptclucut[k] && ptpi0[jj]>ptPi0cut[z]
//				  && s4s9_1[jj]>s4s9cut[s] && s4s9_2[jj]>s4s9cut[s] && iso[jj]>isocut[q] && ( (E_Es_e1_1[jj]+E_Es_e2_1[jj])>elayercut[h] || iseb[jj]==1) && ((E_Es_e1_2[jj]+E_Es_e2_2[jj])>elayercut[h] || iseb[jj]==1) ){
			    if( iseb[jj]==EBorEE && ncris1[jj]>ncri1cut[i] && ncris2[jj]>ncri2cut[j] && ptclu1[jj]>ptclucut[k] && ptclu2[jj]>ptclucut[k] && ptpi0[jj]>ptPi0cut[z]
				  && s4s9_1[jj]>s4s9cut[s] && s4s9_2[jj]>s4s9cut[s] && iso[jj]>isocut[q] ){
				hmass->Fill( binTrue+1, massPi0[jj] );
				if( fabs(etapi0[jj]) < 1.0 && iseb[jj]==1 ) hmassEBL->Fill( binTrue+1, massPi0[jj] );
				if( fabs(etapi0[jj]) > 1.0 && iseb[jj]==1 ) hmassEBH->Fill( binTrue+1, massPi0[jj] );
				if( fabs(etapi0[jj]) < 1.8 && iseb[jj]==2 ) hmassEEL->Fill( binTrue+1, massPi0[jj] );
				if( fabs(etapi0[jj]) > 1.8 && iseb[jj]==2 ) hmassEEH->Fill( binTrue+1, massPi0[jj] );
			    }
			    if( iseb[jj]==EBorEE ){
				hmass_tot->Fill( binTrue+1, massPi0[jj] );
				if( fabs(etapi0[jj]) < 1.0 && iseb[jj]==1 ) hmass_totEBL->Fill( binTrue+1, massPi0[jj] );
				if( fabs(etapi0[jj]) > 1.0 && iseb[jj]==1 ) hmass_totEBH->Fill( binTrue+1, massPi0[jj] );
				if( fabs(etapi0[jj]) < 1.8 && iseb[jj]==2 ) hmass_totEEL->Fill( binTrue+1, massPi0[jj] );
				if( fabs(etapi0[jj]) > 1.8 && iseb[jj]==2 ) hmass_totEEH->Fill( binTrue+1, massPi0[jj] );
			    }
			    binTrue++;
			  }
			}
		    }
		  //}
		}
	    }
	  }
	}
    }
  }
  //Saving Info
  cout << "Start writing histos ..." << endl;
  OutFile->cd(); 
  hmass->Write();
  hmassEBL->Write();
  hmassEBH->Write();
  hmassEEL->Write();
  hmassEEH->Write();
  hmass_tot->Write();
  hmass_totEBL->Write();
  hmass_totEBH->Write();
  hmass_totEEL->Write();
  hmass_totEEH->Write();
  OutFile->Write();
  OutFile->Close();
}
