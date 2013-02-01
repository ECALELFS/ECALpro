#include "CalibCode/CalibTools/interface/ECALGeometry.h"

#include "TFile.h"
#include "TTree.h"
#include<iostream>
using namespace std;

ECALGeometry* ECALGeometry::instance = 0;


ECALGeometry* ECALGeometry::getGeometry(TFile* f) {

  if(instance != 0) return instance;

  instance = new ECALGeometry(f);
  return instance;

}

ECALGeometry::ECALGeometry(TFile* f) {

  if(f==0) {
     cout << "null pointer. Not valid root file to read geometry. exiting..." << endl;
     throw exception();
  }

  instance = new ECALGeometry();

  TTree* tree = (TTree*) f->Get("Geometry");
  uint32_t id;
  float xtalPos[3];
  float xtalAxis[3];

  tree->SetBranchAddress("id",&id);
  tree->SetBranchAddress("xXtal",&xtalPos[0]);
  tree->SetBranchAddress("yXtal",&xtalPos[1]);
  tree->SetBranchAddress("zXtal",&xtalPos[2]);
  tree->SetBranchAddress("xAxisXtal",&xtalAxis[0]);
  tree->SetBranchAddress("yAxisXtal",&xtalAxis[1]);
  tree->SetBranchAddress("zAxisXtal",&xtalAxis[2]);

  for(int i=0; i<tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    posMap[DetId(id)] = GlobalPoint(xtalPos[0], xtalPos[1], xtalPos[2]);
    axisMap[DetId(id)] = GlobalVector(xtalAxis[0], xtalAxis[1], xtalAxis[2]);
  }
  //cout << "finished loading "<< i << "xtals into geometry." << endl;

}

GlobalPoint ECALGeometry::getPosition(DetId id, float depth) const {
  std::map<DetId,GlobalPoint>::const_iterator it = posMap.find(id);
  std::map<DetId,GlobalVector>::const_iterator iax = axisMap.find(id);
  if( it != posMap.end() )
    return it->second + depth*getAxis(id).unit();
  else throw( std::exception() );
}

GlobalVector ECALGeometry::getAxis(DetId id) const {
  std::map<DetId,GlobalVector>::const_iterator iax = axisMap.find(id);
  if(iax != axisMap.end())  return iax->second;
  else throw (std::exception() );
}

void ECALGeometry::print() const {
   
  for(std::map<DetId,GlobalPoint>::const_iterator it = posMap.begin();
                                                  it != posMap.end();
                                                  ++it) {
     DetId i = it->first;
     cout << "Id: " << i << " position(r,eta,phi): " << it->second.mag() << " " << it->second.eta() << " " << it->second.phi()
                    //<< " xtal axis(theta,phi): " << axisMap[i].theta() << " " << axisMap[i].phi()
                    << endl;
  }
}

