#ifndef EcalRecHitCompare_H
#define EcalRecHitCompare_H

class ecalRecHitPtrLess : public std::binary_function<EcalRecHit*, EcalRecHit*, bool>
{
public:
  bool operator()(EcalRecHit* x, EcalRecHit* y)
  {
    return (x->energy() > y->energy());
  }
};

class ecalRecHitLess : public std::binary_function<EcalRecHit, EcalRecHit, bool>
{
public:
  bool operator()(EcalRecHit x, EcalRecHit y)
  {
    return (x.energy() > y.energy());
  }
};


#endif
