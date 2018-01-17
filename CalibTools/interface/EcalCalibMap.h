#ifndef EcalCalibMap_H
#define EcalCalibMap_H

#include <iostream>
#include "unistd.h"
#include <string>

#include "TFile.h"
#include "TH2F.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "CalibCode/CalibTools/interface/ECALGeometry.h"
// #include "CalibCode/FillEpsilonPlot/interface/EndcapTools.h"

//EcalCalibType :: { Xtal, EtaRing, TrigTower };

class EcalCalibMapBase {
    public:
        virtual uint32_t nCoeff() =0;
        virtual float& coeff(const DetId& id) =0;
        virtual float& operator[](const DetId& id) =0;
        virtual const float& operator[](const DetId& id) const =0;
        virtual void loadCalibMapFromFile(const char* cfile = "", const bool useGenPhoton2forEoverEtrue = false) =0;
        virtual int getNRegionsEB() =0;
        virtual int getNRegionsEE() =0;
};


template<typename Type>
class EcalCalibMap : public EcalCalibMapBase {
  public:
    EcalCalibMap() {
     std::cout << "EcalCalibMap::EcalCalibMap" << std::endl;
     for(uint32_t i=0; i<Type::nRegions; ++i) { mapEB[i] = 1.; }
     for(uint32_t i=0; i<Type::nRegionsEE; ++i) { mapEE[i]=1.; }
     bad_coeff = -9999.9999;
    };

    float& coeff(const DetId& id)
    { 
       if(id.subdetId() == EcalBarrel)
          return mapEB[Type::iRegion(id)]; 
       else if(id.subdetId() == EcalEndcap)
          return mapEE[Type::iRegionEE(id)];
       return bad_coeff;
    }

    uint32_t nCoeff() { return nRegionsEB+nRegionsEE; }

    float& operator[](const DetId& id) 
    { 
       if(id.subdetId() == EcalBarrel)
          return mapEB[Type::iRegion(id)]; 
       else if(id.subdetId() == EcalEndcap)
          return mapEE[Type::iRegionEE(id)];
       return bad_coeff;
    }

    const float& operator[](const DetId& id) const 
    { 
       if(id.subdetId() == EcalBarrel)
          return mapEB[Type::iRegion(id)]; 
       else if(id.subdetId() == EcalEndcap)
          return mapEE[Type::iRegionEE(id)];
       return bad_coeff;
    }

    int getNRegionsEB() {return nRegionsEB;}
    int getNRegionsEE() {return nRegionsEE;}

    void loadCalibMapFromFile(const char* cfile = "", const bool useGenPhoton2forEoverEtrue = false); 
    // void setCaloGeometry(ECALGeometry *geom);

  private:
    /* static float mapEB[Type::nRegions]; */
    /* static float mapEE[Type::nRegionsEE]; */
    float mapEB[Type::nRegions];
    float mapEE[Type::nRegionsEE];
    static float bad_coeff;
    static const uint32_t nRegionsEB = Type::nRegions;
    static const uint32_t nRegionsEE = Type::nRegionsEE;
    // EndcapTools eeTool_;

};


template<typename Type> 
void EcalCalibMap<Type>::loadCalibMapFromFile(const char* cfile, const bool useGenPhoton2forEoverEtrue) 
{

  // we are also implementing E/Etrue with MC, in which case we have maps for two different photons with E/Etrue instead of the pi0 mass
  // in this case, useGenPhoton2forEoverEtrue = true implies loading the map for the second gen photon (the first one is stored in the map with the same 
  // name it would have if we were working with pi0 mass

  std::cout << "FIT_EPSILON: [EcalCalibMap] :: photon " << (useGenPhoton2forEoverEtrue ? 2 : 1) << ", loadCalibMapFromFile(" << std::string(cfile) << ") called" << std::endl; 
    TFile* f = TFile::Open(cfile);
    /// keep trying in case of network I/O problems
    for(int iTrial=0; iTrial<10 && !f; iTrial++)
    {
        std::cout << "FIT_EPSILON : [EcalCalibMap] :: Warning: could not open calibMap.root (trial #" << iTrial << ")" << std::endl;
        std::cout << "[EcalCalibMap] ::   '-- Could be a network issue. Trying again in 30s..." << std::endl;
        sleep(30);

        f = TFile::Open(cfile);
    }
    if(!f) throw cms::Exception("LoadCalibMap") << "[EcalCalibMap] :: cannot open calibMap.root after 10 trials in 300s\n";
    else std::cout << "FIT_EPSILON: [EcalCalibMap] :: " << std::string(cfile) << " successfully opened" << std::endl;

    // TH2F* hmap = (TH2F*) f->Get("calibMap");
    TH2F* hmap = nullptr;
    TH2F* hmap_EEp = nullptr;
    TH2F* hmap_EEm = nullptr;
    if (useGenPhoton2forEoverEtrue) {
      hmap = (TH2F*) f->Get("calibMap_EB_g2");
      hmap_EEp = (TH2F*) f->Get("calibMap_EEp_g2");
      hmap_EEm = (TH2F*) f->Get("calibMap_EEm_g2");
      if(!hmap)     throw cms::Exception("LoadCalibMap") << "cannot find TH2F::calibMap_EB_g2 in the file provided\n";
      if(!hmap_EEp) throw cms::Exception("LoadCalibMap") << "cannot find TH2F::calibMap_EEp_g2 in the file provided\n";
      if(!hmap_EEm) throw cms::Exception("LoadCalibMap") << "cannot find TH2F::calibMap_EEm_g2 in the file provided\n";
    } else {
      hmap = (TH2F*) f->Get("calibMap_EB");
      hmap_EEp = (TH2F*) f->Get("calibMap_EEp");
      hmap_EEm = (TH2F*) f->Get("calibMap_EEm");
      if(!hmap)     throw cms::Exception("LoadCalibMap") << "cannot find TH2F::calibMap_EB in the file provided\n";
      if(!hmap_EEp) throw cms::Exception("LoadCalibMap") << "cannot find TH2F::calibMap_EEp in the file provided\n";
      if(!hmap_EEm) throw cms::Exception("LoadCalibMap") << "cannot find TH2F::calibMap_EEm in the file provided\n";
    }

    std::cout << "FIT_EPSILON: loading constants from TH2F::calibMap in <" << cfile << "> ..." << std::endl;

    for(int ix=1; ix<= hmap->GetXaxis()->GetNbins(); ++ix) {
       if(ix==86) continue;
       for(int iphi=1; iphi<= hmap->GetYaxis()->GetNbins(); ++iphi) {
          EBDetId id(ix-1-85,iphi);
          //cout << "id: " << id << " calib coeff: " << hmap->GetBinContent(ix,iphi) << endl;
          this->coeff(id) =  hmap->GetBinContent(ix,iphi);
       }
    }

    std::cout << "FIT_EPSILON: loading constants from TH2F::calibMapEE in <" << cfile << "> ..." << std::endl;

    for(int jR=0; jR<EEDetId::kSizeForDenseIndexing; jR++)
    {
      EEDetId eeid( Type::EEDetIdFromRegion(jR));

      if(eeid.positiveZ())
          this->coeff(eeid) = hmap_EEp->GetBinContent(eeid.ix(), eeid.iy()); 
      else
          this->coeff(eeid) = hmap_EEm->GetBinContent(eeid.ix(), eeid.iy()); 
    }

    std::cout <<  "done." << std::endl;
}


//  template<typename Type> 
//  void EcalCalibMap<Type>::setCaloGeometry(ECALGeometry *geom)
//  {
//      if(!geom) throw cms::Exception("EcalCalibMap::SetCaloGeometry") << "Invalid geometry pointer\n";
//      eeTool_.setCaloGeometry(geom);
//  }


/* template<typename Type> float EcalCalibMap<Type>::mapEB[Type::nRegions]; */
/* template<typename Type> float EcalCalibMap<Type>::mapEE[Type::nRegionsEE]; */
template<typename Type> float EcalCalibMap<Type>::bad_coeff;

#endif
