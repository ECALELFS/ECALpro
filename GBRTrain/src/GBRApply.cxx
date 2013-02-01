#include "CalibCode/GBRTrain/interface/GBRApply.h"
#include "CalibCode/GBRTrain/interface/GBREvent.h"
#include "CalibCode/EgammaObjects/interface/GBRForest.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include <assert.h>
#include <malloc.h>

//_______________________________________________________________________
GBRApply::GBRApply()
{

}

//_______________________________________________________________________
GBRApply::~GBRApply() 
{

}

//_______________________________________________________________________
float GBRApply::ApplyAsFriend_calib( const GBRForest *forest, std::vector<float> vars ) const
{
  Float_t *vals = new Float_t[vars.size()];
  for(unsigned int i=0; i<vars.size(); i++ ){
     vals[i] = vars[i];
  }
  
  float Correction = forest->GetResponse(vals);
  return Correction;
}
//_______________________________________________________________________
TTree *GBRApply::ApplyAsFriend(TTree *intree, const GBRForest *forest, const std::vector<std::string> &vars, std::string targetname) const
{
  
  int nvars = vars.size();
  
  //initialize TTreeFormulas to read variables from TTree
  std::vector<TTreeFormula*> inputforms;
  for (std::vector<std::string>::const_iterator it = vars.begin(); 
      it != vars.end(); ++it) {
    inputforms.push_back(new TTreeFormula(it->c_str(),it->c_str(),intree));
  }
  
  Float_t target = 0.;
  Float_t *vals = new Float_t[nvars];
  
  //initialize new friend tree
  TTree *friendtree = new TTree;
  friendtree->Branch(targetname.c_str(),&target,TString::Format("%s/F",targetname.c_str()));
  
  for (Long64_t iev=0; iev<intree->GetEntries(); ++iev) {
    if (iev%100000==0) printf("%i\n",int(iev));
    intree->LoadTree(iev);
        
    for (int i=0; i<nvars; ++i) {
      vals[i] = inputforms[i]->EvalInstance();
    }
    
    target = forest->GetResponse(vals);
    
    friendtree->Fill();

  }
  
  
  //clear TTreeFormulas
  for (std::vector<TTreeFormula*>::const_iterator it = inputforms.begin(); 
        it != inputforms.end(); ++it) {
      delete *it;
  }
  
  delete[] vals;
  
  intree->AddFriend(friendtree);

  //the branch addresses are set to local variables in this function
  //these local variables go out of scope after this function finishes
  //so we need to reset the branch addresses before returning
  friendtree->ResetBranchAddresses();
  return friendtree;
  
}
