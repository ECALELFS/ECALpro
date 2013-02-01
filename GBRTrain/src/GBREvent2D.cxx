#include "CalibCode/GBRTrain/interface/GBREvent2D.h"

//_______________________________________________________________________
GBREvent2D::GBREvent2D(int nvars) : 
  fVars(new float[nvars]),
  fQuantiles(new int[nvars]),
  fTargetX(0.0),
  fTargetY(0.0),
  fTargetMag(0.0),
  fTransTargetX(0.0),
  fTransTargetY(0.0),
  fWeight(1.0),
  fWeightedTransTargetX(0.),
  fWeightedTransTargetY(0.),
  fWeightedTransTarget2(0.)
{

}

//_______________________________________________________________________
GBREvent2D::~GBREvent2D() 
{
  delete[] fVars;
  delete [] fQuantiles;
  
}
