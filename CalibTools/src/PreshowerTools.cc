//#include "Analysis/Modules/interface/PreshowerCluster.h"
#include <utility>
#include <vector>
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "CalibCode/CalibTools/interface/PreshowerTools.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"


const double PreshowerTools::mip_ = 8.108e-05;
const double PreshowerTools::gamma_ = 0.024;
const double PreshowerTools::calib_planeX_ = 1.0;
const double PreshowerTools::calib_planeY_ = 0.7;
const int    PreshowerTools::clusterwindowsize_ = 15;

//PreshowerTools::PreshowerTools(CaloSubdetectorGeometry* extGeom, CaloSubdetectorTopology* topology_p,  edm::Handle< ESRecHitCollection > & esHandle)
PreshowerTools::PreshowerTools(const CaloGeometry* extGeom, CaloSubdetectorTopology* topology_p,  edm::Handle< ESRecHitCollection > & esHandle) : geom_(extGeom)
{
    //geom_ = extGeom; 
    estopology_ = topology_p; 

    for (ESRecHitCollection::const_iterator it = esHandle->begin(); it != esHandle->end(); it++) {
        //Make the map of DetID, EcalRecHit pairs
        rechits_map.insert(std::make_pair(it->id(), *it));
    }
}



PreshowerCluster PreshowerTools::makeOnePreshowerCluster(int stripwindow, ESDetId *strip)
{
   //the output class
   PreshowerCluster finalcluster;


   // The set of used DetID's
   //HitsID *used_s;

   esroad_2d.clear();


   //used_s = used_strips;
   used_strips.clear();

  
   int plane = strip->plane();

  
  // Collection of cluster strips
   EcalRecHitCollection clusterRecHits;
   // Map of strips for position calculation
   RecHitsMap recHits_pos;

   //Make a navigator, and set it to the strip cell.
   EcalPreshowerNavigator navigator(*strip, estopology_);
   navigator.setHome(*strip);
   //search for neighbours in the central road
   findESRoad(stripwindow,*strip,navigator,plane);

   if ( plane == 1 ) {
      ESDetId strip_north = navigator.north();
      findESRoad(stripwindow,strip_north,navigator,plane);
      navigator.home();
      ESDetId strip_south = navigator.south();
      findESRoad(stripwindow,strip_south,navigator,plane);
      navigator.home();
   }
   if ( plane == 2 ) {
      ESDetId strip_east = navigator.east();
      findESRoad(stripwindow,strip_east,navigator,plane);
      navigator.home();
      ESDetId strip_west = navigator.west();
      findESRoad(stripwindow,strip_west,navigator,plane);
      navigator.home();
   }

   // Start clustering from strip with max Energy in the road
   float E_max = 0.;
   bool found = false;
   RecHitsMap::iterator max_it;
   // Loop over strips:
   std::vector<ESDetId>::iterator itID;
   for (itID = esroad_2d.begin(); itID != esroad_2d.end(); itID++) {
     RecHitsMap::iterator strip_it = rechits_map.find(*itID);   
     if(!goodStrip(strip_it)) continue;

     DetId nonblindstripid (itID->rawId());
     //GlobalPoint position = geom_->getPosition(nonblindstripid);
     //if (position.z() > 0.)eventCont->hist("nonBlindstripPositionZpos")->Fill(position.x(),position.y());
     //if (position.z() < 0.)eventCont->hist("nonBlindstripPositionZneg")->Fill(position.x(),position.y());

     float E = strip_it->second.energy();
     if ( E > E_max) {
        E_max = E;
        found = true;
        max_it = strip_it;
     }
   }
  
	     if ( !found ) {//cout<<"WARNING: HOTSTRIP NOT FOUND!!!"<<endl;

           for (itID = esroad_2d.begin(); itID != esroad_2d.end(); itID++) {
               DetId blindstripid (itID->rawId());
               //GlobalPoint position = geom_->getPosition(blindstripid);
               //cout << "X Position: " << position.x() << "Y Position: " << position.y() << "Z Position: " << position.z() << endl;
               //if (position.z() > 0.)eventCont->hist("BlindstripPositionZpos")->Fill(position.x(),position.y());
               //if (position.z() < 0.)eventCont->hist("BlindstripPositionZneg")->Fill(position.x(),position.y());    
           }

           return finalcluster;}

   // First, save the hottest strip
   clusterRecHits.push_back(max_it->second);  
   recHits_pos.insert(std::make_pair(max_it->first, max_it->second));
   used_strips.insert(max_it->first);

   // Find positions of adjacent strips:
   ESDetId next, strip_1, strip_2;
   navigator.setHome(max_it->first);
   ESDetId startES = max_it->first;
  
   if (plane == 1) {
     // Save two neighbouring strips to the east
     int nadjacents_east = 0;
     while ( (next=navigator.east()) != ESDetId(0) && next != startES && nadjacents_east < 2 ) {
       ++nadjacents_east;
       RecHitsMap::iterator strip_it = rechits_map.find(next);
      
		   if(!goodStrip(strip_it)) continue;
       // Save strip for clustering if it exists, not already in use, and satisfies an energy threshold
        clusterRecHits.push_back(strip_it->second);       
        // save strip for position calculation
        if ( nadjacents_east==1 ) strip_1 = next;
        used_strips.insert(strip_it->first);             
     }
     // Save two neighbouring strips to the west
     navigator.home();
     int nadjacents_west = 0;
     while ( (next=navigator.west()) != ESDetId(0) && next != startES && nadjacents_west < 2 ) {
        ++nadjacents_west;
        RecHitsMap::iterator strip_it = rechits_map.find(next);
        if(!goodStrip(strip_it)) continue;
        clusterRecHits.push_back(strip_it->second);
        if ( nadjacents_west==1 ) strip_2 = next;
        used_strips.insert(strip_it->first);       
     }
   }
  else if (plane == 2) {

   // Save two neighbouring strips to the north
     int nadjacents_north = 0;
     while ( (next=navigator.north()) != ESDetId(0) && next != startES && nadjacents_north < 2 ) {
        ++nadjacents_north; 
        RecHitsMap::iterator strip_it = rechits_map.find(next); 
        if(!goodStrip(strip_it)) continue;      
        clusterRecHits.push_back(strip_it->second);
        if ( nadjacents_north==1 ) strip_1 = next;
        used_strips.insert(strip_it->first);    
     }
     // Save two neighbouring strips to the south
     navigator.home();
     int nadjacents_south = 0;
     while ( (next=navigator.south()) != ESDetId(0) && next != startES && nadjacents_south < 2 ) {
        ++nadjacents_south;   
        RecHitsMap::iterator strip_it = rechits_map.find(next);   
        if(!goodStrip(strip_it)) continue;      
        clusterRecHits.push_back(strip_it->second);
        if ( nadjacents_south==1 ) strip_2 = next;
        used_strips.insert(strip_it->first);    
     }
   }
   else {
     std::cout << " Wrong plane number" << plane <<", null cluster will be returned! " << std::endl;
     return finalcluster;
   } // end of if

   // strips for position calculation
   RecHitsMap::iterator strip_it1, strip_it2;
   if ( strip_1 != ESDetId(0)) {
     strip_it1 = rechits_map.find(strip_1);
     recHits_pos.insert(std::make_pair(strip_it1->first, strip_it1->second));  
   }
   if ( strip_2 != ESDetId(0) ) {
     strip_it2 = rechits_map.find(strip_2);
     recHits_pos.insert(std::make_pair(strip_it2->first, strip_it2->second));  
   }
   
   RecHitsMap::iterator cp;
   double energy_pos = 0;
   double x_pos = 0;
   double y_pos = 0;
   double z_pos = 0;
   for (cp = recHits_pos.begin(); cp!=recHits_pos.end(); cp++ ) {
      double E = cp->second.energy();
      energy_pos += E; 
      GlobalPoint position = geom_->getPosition(cp->first);
      x_pos += E * position.x();
      y_pos += E * position.y();
      z_pos += E * position.z();     
   }
  if(energy_pos>0.) {
     x_pos /= energy_pos;
     y_pos /= energy_pos;
     z_pos /= energy_pos;
  }

  EcalRecHitCollection::iterator it;
  double Eclust = 0;
  int stripscounter = 0;

  for (it=clusterRecHits.begin(); it != clusterRecHits.end(); it++) {
     Eclust += it->energy();
     stripscounter++;
  }

  //Filling PreshowerCluster

  //finalcluster.set_x(x_pos);
  //finalcluster.set_y(y_pos);
  //finalcluster.set_z(z_pos);

  //finalcluster.set_energy(Eclust);
  //finalcluster.set_plane(plane);

  //finalcluster.set_goodcluster(true);
 
  std::vector< std::pair<DetId, float> > usedHits;
  PreshowerCluster output(Eclust,   math::XYZPoint(x_pos,y_pos, z_pos) , usedHits , plane);
  //finalcluster = PreshowerCluster(Eclust,   math::XYZPoint(x_pos,y_pos, z_pos) , std::vector< std::pair<DetId, float> > usedHits, plane);

  //used for debugging purposes 
/*
  cout << "//-------------------------------------------//"<<endl;
  cout << " ES Cluster is created with " << endl;
  cout << " energy = " << finalcluster.get_energy() << endl;
  cout << " plane = " << finalcluster.get_plane() << endl;
  cout << " stripscounter = "<<stripscounter<< endl;
  cout << " (x,y,z) = " << "(" << finalcluster.get_x() <<", "<< finalcluster.get_y() <<", "<< finalcluster.get_z()<<")"<< std::endl; 
  cout << "//-------------------------------------------//"<<endl;
*/

  //return finalcluster;
  return output;

}



void PreshowerTools::findESRoad(int stripwindow, ESDetId strip, EcalPreshowerNavigator theESNav, int plane) {
  
   if ( strip == ESDetId(0) ) return;

    ESDetId next;
    theESNav.setHome(strip);
 // First, add a central strip to the road 
    esroad_2d.push_back(strip);   
   
    if (plane == 1) {
      // east road
      int n_east= 0;
      while ( ((next=theESNav.east()) != ESDetId(0) && next != strip) ) {
         esroad_2d.push_back(next);   
         ++n_east;  
         if (n_east == stripwindow) break; 
      }
      // west road

      int n_west= 0;
      theESNav.home();
      while ( ((next=theESNav.west()) != ESDetId(0) && next != strip )) {
         esroad_2d.push_back(next);   
         ++n_west;  
         if (n_west == stripwindow) break; 
      }
   } 
   else if (plane == 2) {
     // north road
     int n_north= 0;
     while ( ((next=theESNav.north()) != ESDetId(0) && next != strip) ) {       
        esroad_2d.push_back(next);   
        ++n_north;  
        if (n_north == stripwindow) break; 
     }
     // south road
     int n_south= 0;
     theESNav.home();
     while ( ((next=theESNav.south()) != ESDetId(0) && next != strip) ) {
        esroad_2d.push_back(next);   
        ++n_south;  
        if (n_south == stripwindow) break; 
     }
   } 

   theESNav.home();
 }



 // returns true if the candidate strip fulfills the requirements to be added to the cluster:
 //=====================================================================================================
 bool PreshowerTools::goodStrip(RecHitsMap::iterator candidate_it)
 //======================================================================================================
 {
   // crystal should not be included...

      if ( (used_strips.find(candidate_it->first) != used_strips.end())  ||        //...if it already belongs to a cluster
        (candidate_it == rechits_map.end() )                    ||        //...if it corresponds to a hit
        (candidate_it->second.energy() <= 0. ) )   // ...if it has a negative or zero energy
     {
     return false;
     }
      //Used for debug
      /*   if (used_strips.find(candidate_it->first) != used_strips.end()){cout<<"The Strip Already belongs to a cluster"<<endl;
       return false;} else if (candidate_it == rechits_map.end()){cout<<"The Strip correspond to a hit"<<endl;
       return false;} else if (candidate_it->second.energy() <= 0. ) {cout<<"the strip has a negative or zero energy"<<endl;
       return false;}*/
     
   return true;
 }


