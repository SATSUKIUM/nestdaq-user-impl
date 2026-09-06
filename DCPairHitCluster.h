/*
  DCPairHitCluster.hh

  2018/12  K.Shirotori
*/

#ifndef DCPairHitCluster_h
#define DCPairHitCluster_h 1

#include "DCLTrackHit.h"

class DCLTrackHit;

class DCPairHitCluster
{
public:
  DCPairHitCluster( DCLTrackHit *hitA, DCLTrackHit *hitB=0 );
  ~DCPairHitCluster();

private:
  DCLTrackHit *hitA_, *hitB_;
  int nhits_;

public:
  int NumberOfHits( void ) const { return nhits_; }
  DCLTrackHit *GetHit( int i ) const
  {
    if(i==0)      return hitA_;
    else if(i==1) return hitB_;
    else          return 0;
  }

};

#endif // DCPairHitCluster_h
