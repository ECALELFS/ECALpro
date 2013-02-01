#ifndef GlobalFunctions_H
#define GlobalFunctions_H

#include "Rtypes.h"

void convxtalid(Int_t &nphi,Int_t &neta)
{
  // Barrel only
  // Output nphi 0...359; neta 0...84; nside=+1 (for eta>0), or 0 (for eta<0).
  // neta will be [-85,-1] , or [0,84], the minus sign indicates the z<0 side.
  
  if(neta > 0) neta -= 1;
  if(nphi > 359) nphi=nphi-360;
  
} //end of convxtalid

int diff_neta_s(Int_t neta1, Int_t neta2){
  Int_t mdiff;
  mdiff=(neta1-neta2);
  return mdiff;
}

// Calculate the distance in xtals taking into account the periodicity of the Barrel
int diff_nphi_s(Int_t nphi1,Int_t nphi2) {
  Int_t mdiff;
  if(abs(nphi1-nphi2) < (360-abs(nphi1-nphi2))) {
    mdiff=nphi1-nphi2;
  }
  else {
    mdiff=360-abs(nphi1-nphi2);
    if(nphi1>nphi2) mdiff=-mdiff;
  }
  return mdiff;
}


#endif
