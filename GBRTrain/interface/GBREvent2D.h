
#ifndef EGAMMAOBJECTS_GBREvent2D
#define EGAMMAOBJECTS_GBREvent2D

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// GBREvent2D                                                             //
//                                                                      //
// A fast minimal implementation of Gradient-Boosted Regression Trees   //
// which has been especially optimized for size on disk and in memory.  //                                                                  
//                                                                      //
// This is a helper class for GBRTrainer to store  needed information   //
// in memory and facilitate sorting according to target or input        //
// variables.                                                           //
//                                                                      //
//  Josh Bendavid - MIT                                                 //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <stdio.h>
#include <algorithm>
#include <TMath.h>
  
  class GBREvent2D {

    public:

       GBREvent2D(int nvars);       
       ~GBREvent2D();
       
       float Var(int i) const { return fVars[i]; }
       unsigned short Quantile(int i) const { return fQuantiles[i]; }
       float TargetX() const   { return fTargetX;  }
       float TargetY() const   { return fTargetY;  }
       float TargetMag() const   { return fTargetMag;  }
       float TransTargetX() const { return fTransTargetX;  }       
       float TransTargetY() const { return fTransTargetY;  }       
       float Weight() const   { return fWeight;  }
       float WeightedTransTargetX() const { return fWeightedTransTargetX; }
       float WeightedTransTargetY() const { return fWeightedTransTargetY; }
       float WeightedTransTarget2() const { return fWeightedTransTarget2; }
       
       void SetVar(int i, float x) { fVars[i] = x; }
       void SetQuantile(int i, int q) { fQuantiles[i] = q; }
       //cache computed qunatities needed for sorting
       void SetTarget(float x,float y) 
         { fTargetX = x; fTargetY = y;
           fTargetMag = sqrt(fTargetX*fTargetX + fTargetY*fTargetY);
         }
       
       //cache computed qunatities needed for split-search
       void SetTransTarget(float x, float y)
         { fTransTargetX = x; fTransTargetY = y;
           fWeightedTransTargetX = fWeight*fTransTargetX;
           fWeightedTransTargetY = fWeight*fTransTargetY;
           fWeightedTransTarget2 = fWeight*(fTransTargetX*fTransTargetX + fTransTargetY + fTransTargetY);
         }
       void SetWeight(float x)     { fWeight = x; }
       
       
    protected:
      float                    *fVars;
      int                      *fQuantiles;
      float                     fTargetX;
      float                     fTargetY;
      float                     fTargetMag;
      float                     fTransTargetX;
      float                     fTransTargetY;
      float                     fWeight;
      float                     fWeightedTransTargetX;
      float                     fWeightedTransTargetY;
      float                     fWeightedTransTarget2;
  };
  
    
  class GBRAbsTargetCMP : public std::binary_function<GBREvent2D*, GBREvent2D*, bool> {
    public:
      GBRAbsTargetCMP() {}
      bool operator() (const GBREvent2D *ev1, const GBREvent2D *ev2) const { return ev1->TargetMag()<ev2->TargetMag() ? true : false; }
  };  
  
  class GBRVarCMP : public std::binary_function<GBREvent2D*, GBREvent2D*, bool> {
    public:
      GBRVarCMP() {}
      GBRVarCMP(int idx) : fVarIdx(idx) {}      
      bool operator() (const GBREvent2D *ev1, const GBREvent2D *ev2) const { return ev1->Var(fVarIdx)<ev2->Var(fVarIdx) ? true : false; }
      
    protected:
      int fVarIdx;
  };  
  
#endif
