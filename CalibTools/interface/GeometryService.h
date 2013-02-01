#ifndef GeometryService_h
#define GeometryService_h

#include <string>
#include "FWCore/Utilities/interface/Exception.h"
#include "CalibCode/CalibTools/interface/ECALGeometry.h"

using std::string;

class GeometryService
{
    public:
        GeometryService(){}

        static string getGeometryName() {
            if(isNameSet_) return geometryFileName_;
            else throw cms::Exception("GeometryService") << "Geometry file name has not been set\n"; 
        }

        static ECALGeometry* getGeometryPtr() {
            if(isPtrSet_) return geometryPtr_;
            else throw cms::Exception("GeometryService") << "Geometry ptr has not been set\n"; 
        }

        static void setGeometryName(string & inputName) {
            geometryFileName_ = inputName;
            isNameSet_ = true;
        }

        static void setGeometryName(char *inputName) {
            geometryFileName_ = string(inputName);
            isNameSet_ = true;
        }

        static void setGeometryPtr(ECALGeometry *inputPtr) {
            geometryPtr_ = inputPtr;
            isPtrSet_ = true;
        }

    private:
        static string geometryFileName_;
        static ECALGeometry *geometryPtr_; 
        static bool isNameSet_;
        static bool isPtrSet_;
};

#endif
