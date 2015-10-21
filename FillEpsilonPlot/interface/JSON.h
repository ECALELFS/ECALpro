#ifndef JSON_H
#define JSON_H

#include <vector>
#include <map>

typedef std::pair< int, int> aLSSegment;
typedef std::vector< aLSSegment > GoodLSVector;
typedef std::map< int, GoodLSVector  >    LSRange ;
typedef std::pair < int, GoodLSVector > LSRangeElement;

class JSON {
 public:
   JSON(const char* file);

   void ReadJSONFile(const char* json);
   bool isGoodLS(int run, int lumi);

 private:
   int oldRun;   
   LSRange goodLS_;
   LSRange::const_iterator goodLSCache_; // ptr to list of good LS for last run

};
#endif
