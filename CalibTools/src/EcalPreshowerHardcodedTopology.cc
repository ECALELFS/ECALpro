#include "CalibCode/CalibTools/interface/EcalPreshowerHardcodedTopology.h"


ESDetId EcalPreshowerHardcodedTopology::incrementIy(const ESDetId& id) const {

  //Strips orientend along x direction for plane 2
  if (id.plane() == 2)
  {
    if (id.strip() < 32 )
    {
      //Incrementing just strip number
      if (ESDetId::validDetId(id.strip()+1,id.six(),id.siy(),id.plane(),id.zside())) 
         return ESDetId(id.strip()+1,id.six(),id.siy(),id.plane(),id.zside());
      else 
         return ESDetId(0); 
    }
    else
    {
      //Changing wafer
      if (ESDetId::validDetId(1,id.six(),id.siy()+1,id.plane(),id.zside())) 
         return ESDetId(1,id.six(),id.siy()+1,id.plane(),id.zside());
      else 
         return ESDetId(0);
    }
  }
  //Strips orientend along y direction for plane 1
  else if (id.plane() == 1)
  {
     //Changing wafer
     if (ESDetId::validDetId(id.strip(),id.six(),id.siy()+1,id.plane(),id.zside()))
        return ESDetId(id.strip(),id.six(),id.siy()+1,id.plane(),id.zside());
     else
        return ESDetId(0);
  }
  else
    return ESDetId(0);
} 



ESDetId EcalPreshowerHardcodedTopology::decrementIy(const ESDetId& id) const {
  //Strips orientend along x direction for plane 2
  if (id.plane() == 2)
  {
     if (id.strip() >1 )
     {
       //Decrementing just strip number
       if (ESDetId::validDetId(id.strip()-1,id.six(),id.siy(),id.plane(),id.zside()))
         return ESDetId(id.strip()-1,id.six(),id.siy(),id.plane(),id.zside());
       else
         return ESDetId(0);
     }
     else
     {
       //Changing wafer
       if (ESDetId::validDetId(32,id.six(),id.siy()-1,id.plane(),id.zside()))
         return ESDetId(32,id.six(),id.siy()-1,id.plane(),id.zside());
       else
         return ESDetId(0);
     }
  }
  //Strips orientend along y direction for plane 1
  else if (id.plane() == 1)
  {
     //Changing wafer
     if (ESDetId::validDetId(id.strip(),id.six(),id.siy()-1,id.plane(),id.zside()))
         return ESDetId(id.strip(),id.six(),id.siy()-1,id.plane(),id.zside());
     else
        return ESDetId(0);
  }
  else
    return ESDetId(0);
} 



ESDetId EcalPreshowerHardcodedTopology::incrementIx(const ESDetId& id) const {
  //Strips orientend along x direction for plane 2
  if (id.plane() == 2)
  {
     //Changing wafer
     if (ESDetId::validDetId(id.strip(),id.six()+1,id.siy(),id.plane(),id.zside()))
        return ESDetId(id.strip(),id.six()+1,id.siy(),id.plane(),id.zside());
     else
        return ESDetId(0);
  }
  //Strips orientend along y direction for plane 1
  else if (id.plane() == 1)
  {
     if (id.strip() < 32 )
     {
        //Incrementing just strip number
        if (ESDetId::validDetId(id.strip()+1,id.six(),id.siy(),id.plane(),id.zside())) 
           return ESDetId(id.strip()+1,id.six(),id.siy(),id.plane(),id.zside());
        else
           return ESDetId(0);
     }
     else
     {
        //Changing wafer
        if (ESDetId::validDetId(1,id.six()+1,id.siy(),id.plane(),id.zside()))
           return ESDetId(1,id.six()+1,id.siy(),id.plane(),id.zside());
        else
           return ESDetId(0);
     }
  }
  else
    return ESDetId(0);
} 


ESDetId EcalPreshowerHardcodedTopology::decrementIx(const ESDetId& id) const {

  //Strips orientend along x direction for plane 2
  if (id.plane() == 2)
  {
      //Changing wafer
      if (ESDetId::validDetId(id.strip(),id.six()-1,id.siy(),id.plane(),id.zside()))
         return ESDetId(id.strip(),id.six()-1,id.siy(),id.plane(),id.zside());
      else
         return ESDetId(0);
    }
  //Strips orientend along y direction for plane 1
  else if (id.plane() == 1)
  {
     if (id.strip() > 1 )
     {
       //Decrementing just strip number
       if (ESDetId::validDetId(id.strip()-1,id.six(),id.siy(),id.plane(),id.zside()))
         return ESDetId(id.strip()-1,id.six(),id.siy(),id.plane(),id.zside());
       else
         return ESDetId(0);
     }
     else
     {
       //Changing wafer
       if (ESDetId::validDetId(32,id.six()-1,id.siy(),id.plane(),id.zside()))
          return ESDetId(32,id.six()-1,id.siy(),id.plane(),id.zside());
       else
          return ESDetId(0);
     }
  }
  else
    return ESDetId(0);
} 



ESDetId EcalPreshowerHardcodedTopology::incrementIz(const ESDetId& id) const {

  if (ESDetId::validDetId(id.strip(),id.six(),id.siy(),id.plane()+1,id.zside()))
     return ESDetId(id.strip(),id.six(),id.siy(),id.plane()+1,id.zside());
  else
     return ESDetId(0);
} 


ESDetId EcalPreshowerHardcodedTopology::decrementIz(const ESDetId& id) const {

  if (ESDetId::validDetId(id.strip(),id.six(),id.siy(),id.plane()-1,id.zside()))
     return ESDetId(id.strip(),id.six(),id.siy(),id.plane()-1,id.zside());
  else
    return ESDetId(0);
  
} 
