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
        virtual void loadCalibMapFromFile(const char* cfile) =0;
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

    void loadCalibMapFromFile(const char* cfile); 
    // void setCaloGeometry(ECALGeometry *geom);

  private:
    static float mapEB[Type::nRegions];
    static float mapEE[Type::nRegionsEE];
    static float bad_coeff;
    static const uint32_t nRegionsEB = Type::nRegions;
    static const uint32_t nRegionsEE = Type::nRegionsEE;
    // EndcapTools eeTool_;

};


template<typename Type> 
void EcalCalibMap<Type>::loadCalibMapFromFile(const char* cfile) 
{

    std::cout << "[EcalCalibMap] :: loadCalibMapFromFile(" << std::string(cfile) << ") called" << std::endl; 
    TFile* f = TFile::Open(cfile);
    /// keep trying in case of network I/O problems
    for(int iTrial=0; iTrial<10 && !f; iTrial++)
    {
        std::cout << "[EcalCalibMap] :: Warning: could not open calibMap.root (trial #" << iTrial << ")" << std::endl;
        std::cout << "[EcalCalibMap] ::   '-- Could be a network issue. Trying again in 30s..." << std::endl;
        sleep(30);

        f = TFile::Open(cfile);
    }
    if(!f) throw cms::Exception("LoadCalibMap") << "[EcalCalibMap] :: cannot open calibMap.root after 10 trials in 300s\n";
    else std::cout << "[EcalCalibMap] :: " << std::string(cfile) << " successfully opened" << std::endl;

    // TH2F* hmap = (TH2F*) f->Get("calibMap");
    TH2F* hmap = (TH2F*) f->Get("calibMap_EB");
    TH2F* hmap_EEp = (TH2F*) f->Get("calibMap_EEp");
    TH2F* hmap_EEm = (TH2F*) f->Get("calibMap_EEm");
    if(!hmap)     throw cms::Exception("LoadCalibMap") << "cannot find TH2F::calibMap in the file provided\n";
    if(!hmap_EEp) throw cms::Exception("LoadCalibMap") << "cannot find TH2F::hmap_EEp in the file provided\n";
    if(!hmap_EEm) throw cms::Exception("LoadCalibMap") << "cannot find TH2F::hmap_EEm in the file provided\n";

    std::cout << "loading constants from TH2F::calibMap in <" << cfile << "> ..." << std::endl;

    for(int ix=1; ix<= hmap->GetXaxis()->GetNbins(); ++ix) {
       if(ix==86) continue;
       for(int iphi=1; iphi<= hmap->GetYaxis()->GetNbins(); ++iphi) {
          EBDetId id(ix-1-85,iphi);
          //cout << "id: " << id << " calib coeff: " << hmap->GetBinContent(ix,iphi) << endl;
          this->coeff(id) =  hmap->GetBinContent(ix,iphi);
       }
    }

    std::cout << "loading constants from TH2F::calibMapEE in <" << cfile << "> ..." << std::endl;

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


  template<typename Type> float EcalCalibMap<Type>::mapEB[Type::nRegions];
  template<typename Type> float EcalCalibMap<Type>::mapEE[Type::nRegionsEE];
  template<typename Type> float EcalCalibMap<Type>::bad_coeff;

#endif
