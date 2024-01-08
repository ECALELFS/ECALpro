// $Id: DumpCaloGeometry.h,v 1.1 2012/01/31 11:41:34 mgrassi Exp $
// Shahram Rahatlou, Sapienza Universita` di Roma & INFN
// May 2010
//
#ifndef DumpCaloGeometry_H
#define DumpCaloGeometry_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

class TFile;
class TTree;

class DumpCaloGeometry : public edm::one::EDAnalyzer <>{ 

public:

  DumpCaloGeometry( const edm::ParameterSet& );
  ~DumpCaloGeometry();

protected:
   
  void beginRun(const edm::Run& r, const edm::EventSetup& c);

  void analyze(const edm::Event& e, const edm::EventSetup& c) ;

  void beginLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
                            const edm::EventSetup& context) ;

  void endLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
                          const edm::EventSetup& c);

  void endRun(const edm::Run& r, const edm::EventSetup& c);

private:

  float xtalAxis[3];
  float xtalPos[3];
  uint32_t   id;

  /// file and ttree objects : Pi0 EB
  TFile* m_file;
  TTree* m_tree;
  std::string m_outfilename;
  
  // cbasile [CMSSW_13_3_X]: declare the token here to access calo geometry via esConsumes() in the costructor
  const edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geoToken_;
};

#endif

