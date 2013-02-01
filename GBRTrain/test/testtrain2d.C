#include "TFile.h"
#include "TTree.h"
#include "GBRTrainer2D.h"
#include "GBRForest2D.h"
#include "Cintex/Cintex.h"


void testtrain2d() {
 
  TFile *fmcele = new TFile("/scratch/bendavid/cms/hist/hgg-v0/merged/hgg-v0_f11-zjets-v14b-pu_noskim.root","READ");  
  TDirectory *direle = (TDirectory*)fmcele->FindObjectAny("PhotonTreeWriterE");
  TTree *hmcele = (TTree*)direle->Get("hPhotonTree");
  TDirectory *direlesingle = (TDirectory*)fmcele->FindObjectAny("PhotonTreeWriterE");
  TTree *hmcelesingle = (TTree*)direlesingle->Get("hPhotonTreeSingle");  
  
  
  GBRTrainer2D *train = new GBRTrainer2D;
  train->SetTree(hmcelesingle);
  train->SetTrainingCut("ph.ispromptgen && ph.isbarrel && ph.hasphoton && ph.pt>25 && evt%20==0");
  train->SetMinEvents(200);
  train->SetShrinkage(0.1);
    
 
  //regression target is the vector pT of the supercluster
  
  //simpler version in absolute detector coordinates
  //train->SetTargetXVar("(ph.genpt*cos(ph.genphi) - ph.scrawe*cos(ph.scphi)/cosh(ph.sceta))/(ph.scrawe/cosh(ph.sceta))");
  //train->SetTargetYVar("(ph.genpt*sin(ph.genphi) - ph.scrawe*sin(ph.scphi)/cosh(ph.sceta))/(ph.scrawe/cosh(ph.sceta))");

  //better version rotating in phi event by event to line up with the (uncorrected) supercluster phi by applying the
  //transformation phi->(phi-scphi)
  //expanded form is written out here, but you can see that the second term for the pll component is just 1
  //if the phi measurement would be perfect, then the perpendicular target component is strictly 0 and
  //the whole thing reduces nearly identically to the 1d response regression (modulo some places where the mean is computed instead of the median)
//   train->SetTargetXVar("(ph.genpt*cos(ph.genphi-ph.scphi) - ph.scrawe/cosh(ph.sceta))/(ph.scrawe/cosh(ph.sceta))");
//   train->SetTargetYVar("ph.genpt*sin(ph.genphi-ph.scphi)/(ph.scrawe/cosh(ph.sceta))");

  //even better version, adapted from the above, but removing the unnecessary subtraction of the uncorrected supercluster vector
  train->SetTargetXVar("ph.genpt*cos(ph.genphi-ph.scphi)/(ph.scrawe/cosh(ph.sceta))");
  train->SetTargetYVar("ph.genpt*sin(ph.genphi-ph.scphi)/(ph.scrawe/cosh(ph.sceta))");
  
  
  std::vector<std::string> *varsf = new std::vector<std::string>;
  
  varsf->push_back("ph.scrawe");
  varsf->push_back("ph.r9");
  varsf->push_back("ph.sceta");
  varsf->push_back("ph.scphi");
  
  varsf->push_back("ph.e5x5/ph.scrawe");
  varsf->push_back("ph.hovere");
  varsf->push_back("ph.scetawidth");
  varsf->push_back("ph.scphiwidth");
  varsf->push_back("ph.scnclusters");
  
  varsf->push_back("nVtx");
  

/*  varsf->push_back("ph.etaseed-ph.sceta");
  varsf->push_back("atan2(sin(ph.phiseed-ph.scphi),cos(ph.phiseed-ph.scphi))");
  varsf->push_back("ph.eseed/ph.scrawe");
  varsf->push_back("ph.e3x3seed/ph.eseed");
  varsf->push_back("ph.e5x5seed/ph.eseed");
  varsf->push_back("ph.sigietaietaseed");   
  varsf->push_back("ph.sigiphiphiseed");   
  varsf->push_back("ph.covietaiphiseed");
  varsf->push_back("ph.emaxseed/ph.eseed");
  varsf->push_back("log(ph.e2ndseed/ph.emaxseed)");
  varsf->push_back("log(ph.etopseed/ph.emaxseed)");
  varsf->push_back("log(ph.ebottomseed/ph.emaxseed)");
  varsf->push_back("log(ph.eleftseed/ph.emaxseed)");
  varsf->push_back("log(ph.erightseed/ph.emaxseed)");
  varsf->push_back("(ph.etopseed-ph.ebottomseed)/(ph.etopseed+ph.ebottomseed)");
  varsf->push_back("(ph.eleftseed-ph.erightseed)/(ph.eleftseed+ph.erightseed)");

  
  varsf->push_back("(ph.ebc2>0.0)*(ph.etabc2-ph.sceta)");
  varsf->push_back("(ph.ebc2>0.0)*atan2(sin(ph.phibc2-ph.scphi),cos(ph.phibc2-ph.scphi))");
  varsf->push_back("(ph.ebc2>0.0)*ph.ebc2/ph.scrawe");
  varsf->push_back("(ph.ebc2>0.0)*ph.e3x3bc2/ph.ebc2");
  varsf->push_back("(ph.ebc2>0.0)*ph.e5x5bc2/ph.ebc2");
  varsf->push_back("(ph.ebc2>0.0)*ph.sigietaietabc2");   
  varsf->push_back("(ph.ebc2>0.0)*ph.sigiphiphibc2");   
  varsf->push_back("(ph.ebc2>0.0)*ph.covietaiphibc2");
  varsf->push_back("(ph.ebc2>0.0)*ph.emaxbc2/ph.ebc2");
  varsf->push_back("(ph.ebc2>0.0)*log(ph.e2ndbc2/ph.emaxbc2)");
  varsf->push_back("(ph.ebc2>0.0)*log(ph.etopbc2/ph.emaxbc2)");
  varsf->push_back("(ph.ebc2>0.0)*log(ph.ebottombc2/ph.emaxbc2)");
  varsf->push_back("(ph.ebc2>0.0)*log(ph.eleftbc2/ph.emaxbc2)");
  varsf->push_back("(ph.ebc2>0.0)*log(ph.erightbc2/ph.emaxbc2)");
  varsf->push_back("(ph.ebc2>0.0)*(ph.etopbc2-ph.ebottombc2)/(ph.etopbc2+ph.ebottombc2)");
  varsf->push_back("(ph.ebc2>0.0)*(ph.eleftbc2-ph.erightbc2)/(ph.eleftbc2+ph.erightbc2)");
  

  varsf->push_back("(ph.ebclast>0.0)*(ph.etabclast-ph.sceta)");
  varsf->push_back("(ph.ebclast>0.0)*atan2(sin(ph.phibclast-ph.scphi),cos(ph.phibclast-ph.scphi))");
  varsf->push_back("(ph.ebclast>0.0)*ph.ebclast/ph.scrawe");
  varsf->push_back("(ph.ebclast>0.0)*ph.e3x3bclast/ph.ebclast");
  varsf->push_back("(ph.ebclast>0.0)*ph.e5x5bclast/ph.ebclast");
  varsf->push_back("(ph.ebclast>0.0)*ph.sigietaietabclast");   
  varsf->push_back("(ph.ebclast>0.0)*ph.sigiphiphibclast");   
  varsf->push_back("(ph.ebclast>0.0)*ph.covietaiphibclast");

  //2nd last basiccluster variables (generic barrel and endcap)
  varsf->push_back("(ph.ebclast2>0.0)*(ph.etabclast2-ph.sceta)");
  varsf->push_back("(ph.ebclast2>0.0)*atan2(sin(ph.phibclast2-ph.scphi),cos(ph.phibclast2-ph.scphi))");
  varsf->push_back("(ph.ebclast2>0.0)*ph.ebclast2/ph.scrawe");
  varsf->push_back("(ph.ebclast2>0.0)*ph.e3x3bclast2/ph.ebclast2");
  varsf->push_back("(ph.ebclast2>0.0)*ph.e5x5bclast2/ph.ebclast2");
  varsf->push_back("(ph.ebclast2>0.0)*ph.sigietaietabclast2");   
  varsf->push_back("(ph.ebclast2>0.0)*ph.sigiphiphibclast2");   
  varsf->push_back("(ph.ebclast2>0.0)*ph.covietaiphibclast2");

  varsf->push_back("ph.ietaseed");
  varsf->push_back("ph.iphiseed");
  varsf->push_back("ph.ietaseed%5");
  varsf->push_back("ph.iphiseed%2");       
  varsf->push_back("(abs(ph.ietaseed)<=25)*(ph.ietaseed%25) + (abs(ph.ietaseed)>25)*((ph.ietaseed-25*abs(ph.ietaseed)/ph.ietaseed)%20)");
  varsf->push_back("ph.iphiseed%20"); 
  varsf->push_back("ph.etacryseed");
  varsf->push_back("ph.phicryseed");
    
  varsf->push_back("ph.ietabc2");
  varsf->push_back("ph.iphibc2");
  varsf->push_back("ph.ietabc2%5");
  varsf->push_back("ph.iphibc2%2");       
  varsf->push_back("(abs(ph.ietabc2)<=25)*(ph.ietabc2%25) + (abs(ph.ietabc2)>25)*((ph.ietabc2-25*abs(ph.ietabc2)/ph.ietabc2)%20)");
  varsf->push_back("ph.iphibc2%20"); 
  varsf->push_back("ph.etacrybc2");
  varsf->push_back("ph.phicrybc2");     */
  
  for (int i=0; i<varsf->size(); ++i) {
    train->AddInputVar(varsf->at(i));
  }
   


   
  
  ROOT::Cintex::Cintex::Enable();   
  const GBRForest2D *forest = train->TrainForest(100);

  
  TFile *fout = new TFile("fgbrtraintest2d.root","RECREATE");    
  fout->WriteObject(forest,"gbrtrain");
  fout->WriteObject(varsf, "varlist");
  
  
  fout->Close();
  
}