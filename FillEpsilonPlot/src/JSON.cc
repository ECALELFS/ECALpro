#include "CalibCode/FillEpsilonPlot/interface/JSON.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
using std::cout;
using std::endl;
using std::vector;
using std::map;
using std::pair;
using std::stringstream;
using std::runtime_error;

#include "CalibCode/json/reader.h"
#include "CalibCode/json/elements.h"

JSON::JSON(const char* json) {

  goodLS_ = LSRange();
  goodLSCache_ = goodLS_.end();


  std::ifstream jsonFileStream;
  jsonFileStream.open(json);
  if (!jsonFileStream.is_open())
  {
    std::cout << "Unable to open file " << json << std::endl;
    return;
  }

  json::Object elemRootFile;
  json::Reader::Read(elemRootFile, jsonFileStream);

  for (json::Object::const_iterator itRun=elemRootFile.Begin();itRun!=elemRootFile.End();++itRun)
  {
    const json::Array& lsSegment = (*itRun).element;
    GoodLSVector thisRunSegments;
    for (json::Array::const_iterator lsIterator=lsSegment.Begin();lsIterator!=lsSegment.End();++lsIterator)
    {
	json::Array lsSegment=(*lsIterator);
	json::Number lsStart=lsSegment[0];
	json::Number lsEnd=lsSegment[1];
	aLSSegment thisSegment=std::make_pair<int,int>(lsStart.Value(),lsEnd.Value());
	thisRunSegments.push_back(thisSegment);
    }
    goodLS_[atoi((*itRun).name.c_str())]=thisRunSegments;
  }


  std::cout << "[GoodRunLSMap]::Good Run LS map filled with " << goodLS_.size() << " runs" << std::endl;
  for (LSRange::const_iterator itR=goodLS_.begin(); itR!=goodLS_.end(); ++itR)
  {
    std::cout << "[GoodRunLSMap]::Run " << (*itR).first <<  " LS ranges are: ";
    for (GoodLSVector::const_iterator iSeg=(*itR).second.begin();iSeg!=(*itR).second.end();++iSeg)
	std::cout << "[" << (*iSeg).first << "," << (*iSeg).second << "] ";
    std::cout << std::endl;
  }  
}

//========================
bool JSON::isGoodLS(int run, int lumi) {
  //========================
  if( oldRun != run ) {
    oldRun = run;
    goodLSCache_ = goodLS_.find( run );
  }
  // check whether this run is part of the good runs. else retrun false
  if( goodLSCache_ != goodLS_.end() ) {
    // get list of LS intervals
    const GoodLSVector& lsvector =   goodLSCache_->second; 
    // loop over good LS intervals and return as soon as one interval contains this event
    for(GoodLSVector::const_iterator iLS = lsvector.begin(); iLS != lsvector.end(); iLS++) {
	if(lumi >= iLS->first && lumi <= iLS->second ) {
	  //cout << "Accepting run: " << Run << " LS: " << LumiSection << endl; 
	  return true;
	} // check current LS being in the interval
    } // loop over good LS for this run
  }
  return false;
}
