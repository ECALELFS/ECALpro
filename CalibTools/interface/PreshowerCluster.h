#ifndef PreshowerCluster_H
#define PreshowerCluster_H

class PreshowerCluster {
 public: 
  
  PreshowerCluster(){
    x = 0.;
    y = 0.;
    z = 0.;
    energy = 0.;
    plane = 0.;
    goodcluster = false;
  };
  
  ~PreshowerCluster(){};

  double get_x(){return x;}
  double get_y(){return y;}
  double get_z(){return z;}
  double get_energy(){return energy;}
  int get_plane(){return plane;}
  bool get_goodcluster(){return goodcluster;}

  void set_x(double a) { x = a;}
  void set_y(double a) { y = a;}
  void set_z(double a) { z = a;}
  void set_energy(double a) { energy = a;}
  void set_plane(int a){plane = a;}
  void set_goodcluster(bool a){goodcluster = a;}

  

 private:

  double x;
  double y;
  double z;

  bool goodcluster;
  double energy;
  int plane;
  
  };
#endif
